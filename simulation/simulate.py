"""
simulate.py
Main simulation loop for Arctic PFOA transport.

Two-layer horizontal advection driven by TOPAZ4 monthly velocity fields.
  Layer 1: 0-200m   (depth indices 0-19), depth-weighted mean velocity
  Layer 2: 250-4000m (depth indices 20-39), depth-weighted mean velocity

Timestep: 6 hours (~122 steps per monthly forcing interval)

Sink term: first-order decay with 5-year half-life, applied each timestep.
  k = ln(2) / (5 * 365.25 * 86400) = 4.40e-9 s^-1
  Rationale: lumped representation of sediment burial, degradation, and
  sea ice partitioning. Keeps simulation bounded while delegating spatial
  and temporal structure of removal to PINN bias correction. Consistent
  with Ice-BCNet (Yuan et al. 2024) approach of omitting detailed sinks
  from the forward model.

Output (Option B):
  - Monthly layer-mean concentration fields (nlat x nlon)
  - Monthly concentration at observation station locations
"""

import os
import sys
import numpy as np
import netCDF4
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from grid_utils import load_grid
from advection import upwind_step
from sources import (load_riverine, load_atmospheric,
                     get_riverine_source, get_atmospheric_source,
                     get_boundary_1d)

# ------------------------------------------------------------------
# paths (relative to repo root)
# ------------------------------------------------------------------
TOPAZ4_PATH  = 'data/topaz4_arctic_velocity_2004_2025.nc'
RIVERINE_CSV = 'data/riverine_pfoa_flux.csv'
ATMOS_CSV    = 'data/atmospheric_pfoa_deposition.csv'
OBS_CSV      = 'data/PFAS In-Situ Observations - Sheet1.csv'
OUTPUT_PATH  = 'output/simulation_pfoa.nc'

# ------------------------------------------------------------------
# simulation parameters
# ------------------------------------------------------------------
DT          = 1.0 * 3600.0    # 1-hour timestep in seconds
N_MONTHS    = 264              # 2004-01 through 2025-12
KAPPA_V     = 1e-5             # m2/s, background diapycnal diffusivity
H1          = 200.0            # m, layer 1 thickness (fixed)

# first-order decay sink: 5-year half-life
# k = ln(2) / half_life_seconds
HALF_LIFE_YR  = 94.0
DECAY_K       = np.log(2.0) / (HALF_LIFE_YR * 365.25 * 86400.0)  # s^-1
DECAY_FACTOR  = np.exp(-DECAY_K * DT)   # dimensionless, applied each 6hr step
# DECAY_FACTOR ~ 0.9999997 per step, ~0.9986 per month, ~0.870 per year


def _dz(depth):
    """Layer thicknesses from TOPAZ4 depth level centres."""
    depth = np.asarray(depth, dtype=float)
    mids  = 0.5 * (depth[:-1] + depth[1:])
    dz        = np.empty(len(depth))
    dz[0]     = mids[0]
    dz[1:-1]  = mids[1:] - mids[:-1]
    dz[-1]    = depth[-1] - mids[-1]
    return dz


def layer_mean_velocity(vx_all, vy_all, layer_idxs, dz):
    """
    Depth-weighted mean velocity over a layer.

    Parameters
    ----------
    vx_all, vy_all : (ndepth, nlat, nlon) float, NaN over land
    layer_idxs     : slice
    dz             : (ndepth,) float, level thicknesses in metres

    Returns
    -------
    vx_mean, vy_mean : (nlat, nlon), NaN over land
    """
    vx = vx_all[layer_idxs]
    vy = vy_all[layer_idxs]
    w  = dz[layer_idxs][:, np.newaxis, np.newaxis]

    all_land = np.all(np.isnan(vx), axis=0)
    vx_mean  = np.sum(np.where(np.isnan(vx), 0.0, vx) * w, axis=0) / w.sum()
    vy_mean  = np.sum(np.where(np.isnan(vy), 0.0, vy) * w, axis=0) / w.sum()

    vx_mean[all_land] = np.nan
    vy_mean[all_land] = np.nan

    return vx_mean, vy_mean


def apply_diffusion(C1, C2, grid, dt):
    """
    Diffusive exchange between layer 1 (0-200m) and layer 2 (200m+).

    Physics:
        dC1/dt = -kappa_v * (C1 - C2) / (dz_grad * H1)
        dC2/dt = +kappa_v * (C1 - C2) / (dz_grad * H2)
    where dz_grad = H1/2 + H2/2 (gradient length scale between layer centres).
    kappa_v = 1e-5 m2/s (background diapycnal diffusivity).
    """
    H2      = np.maximum(grid['model_depth'] - H1, 0.0)
    active  = (~grid['land_mask']) & (H2 > 0.0)

    H2_safe  = np.where(H2 > 0.0, H2, 1.0)
    dz_grad  = np.where(active, H1 / 2.0 + H2_safe / 2.0, 1.0)
    flux     = np.where(active, KAPPA_V * (C1 - C2) / dz_grad, 0.0)

    C1_new = C1 - flux / H1 * dt
    C2_new = C2 + np.where(active, flux / H2_safe * dt, 0.0)

    np.maximum(C1_new, 0.0, out=C1_new)
    np.maximum(C2_new, 0.0, out=C2_new)
    C1_new[grid['land_mask']] = 0.0
    C2_new[grid['land_mask']] = 0.0

    return C1_new, C2_new


def load_stations(obs_csv, grid):
    df = pd.read_csv(obs_csv)
    df = df[(df['compound'] == 'PFOA') & (df['lat'] >= 60.0)].copy()
    st = (df[['lat', 'lon']]
          .drop_duplicates()
          .reset_index(drop=True)
          .rename(columns={'lat': 'obs_lat', 'lon': 'obs_lon'}))
    st['station_id'] = st.index

    # snap to nearest OCEAN cell (not land)
    ocean_i, ocean_j = np.where(~grid['land_mask'])

    def snap_to_ocean(obs_lat, obs_lon):
        dlat = grid['lat'][ocean_i] - obs_lat
        dlon = grid['lon'][ocean_j] - obs_lon
        # handle longitude wrapping
        dlon = np.where(dlon > 180, dlon - 360, dlon)
        dlon = np.where(dlon < -180, dlon + 360, dlon)
        dist = dlat**2 + dlon**2
        k = np.argmin(dist)
        return int(ocean_i[k]), int(ocean_j[k])

    ij = st.apply(lambda row: snap_to_ocean(row['obs_lat'], row['obs_lon']), axis=1)
    st['i'] = [x[0] for x in ij]
    st['j'] = [x[1] for x in ij]
    return st

def create_output(path, grid, stations):
    """Initialise output NetCDF4 file. Returns open Dataset."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ds = netCDF4.Dataset(path, 'w', format='NETCDF4')

    ds.createDimension('time',    N_MONTHS)
    ds.createDimension('lat',     grid['nlat'])
    ds.createDimension('lon',     grid['nlon'])
    ds.createDimension('station', len(stations))

    t = ds.createVariable('time', 'i4', ('time',))
    t.units    = 'months since 2004-01-01'
    t[:]       = np.arange(N_MONTHS)

    la = ds.createVariable('lat', 'f4', ('lat',))
    la.units = 'degrees_north'; la[:] = grid['lat']

    lo = ds.createVariable('lon', 'f4', ('lon',))
    lo.units = 'degrees_east';  lo[:] = grid['lon']

    sl = ds.createVariable('station_lat', 'f4', ('station',))
    sl[:] = stations['obs_lat'].values
    so = ds.createVariable('station_lon', 'f4', ('station',))
    so[:] = stations['obs_lon'].values

    for name, desc in [
        ('C_layer1',     'PFOA concentration layer 1 mean 0-200m ng/L'),
        ('C_layer2',     'PFOA concentration layer 2 mean 200-4000m ng/L'),
    ]:
        v = ds.createVariable(name, 'f4', ('time', 'lat', 'lon'),
                              fill_value=np.nan, zlib=True, complevel=4)
        v.units = 'ng/L'; v.long_name = desc

    for name, desc in [
        ('C_station_L1', 'PFOA at obs stations layer 1 ng/L'),
        ('C_station_L2', 'PFOA at obs stations layer 2 ng/L'),
    ]:
        v = ds.createVariable(name, 'f4', ('time', 'station'), fill_value=np.nan)
        v.units = 'ng/L'; v.long_name = desc

    ds.description  = 'Arctic PFOA 2-layer upwind advection simulation, TOPAZ4 forcing'
    ds.kappa_v      = f'{KAPPA_V} m2/s diapycnal diffusivity'
    ds.dt_seconds   = str(DT)
    ds.decay_k      = f'{DECAY_K:.4e} s^-1 (5-yr half-life lumped sink)'
    ds.layer1       = '0-200m depth-weighted mean velocity'
    ds.layer2       = '250-4000m depth-weighted mean velocity'

    return ds


def run():
    print('loading grid...')
    grid = load_grid(TOPAZ4_PATH)
    dz   = _dz(grid['depth'])

    print(f'decay factor per step: {DECAY_FACTOR:.8f}  '
          f'(half-life = {HALF_LIFE_YR} yr, k = {DECAY_K:.4e} s^-1)')

    print('loading sources...')
    riv_df = load_riverine(RIVERINE_CSV)
    atm_df = load_atmospheric(ATMOS_CSV)

    print('loading observation stations...')
    stations = load_stations(OBS_CSV, grid)
    print(f'  {len(stations)} unique PFOA stations above 60N')

    # static boundary concentrations
    c_bnd_L1 = get_boundary_1d(grid, 100.0)
    c_bnd_L2 = get_boundary_1d(grid, 250.0)

    # initialise concentration fields to zero
    C1 = np.zeros((grid['nlat'], grid['nlon']))
    C2 = np.zeros((grid['nlat'], grid['nlon']))

    steps_per_month = int(round((365.25 / 12.0 * 86400.0) / DT))
    print(f'steps per month: {steps_per_month}  ({steps_per_month * DT / 86400:.1f} days)')

    print(f'creating output: {OUTPUT_PATH}')
    out = create_output(OUTPUT_PATH, grid, stations)

    topaz = netCDF4.Dataset(TOPAZ4_PATH, 'r')
    lat_s = grid['lat_idx_min']
    lat_e = grid['lat_idx_max']

    print('running simulation...')
    for t_idx in range(N_MONTHS):
        year  = 2004 + t_idx // 12
        month = 1    + t_idx  % 12

        if month == 1:
            print(f'  {year}  (C1 max={C1.max():.4f} ng/L  mean={C1[~grid["land_mask"]].mean():.4f})')

        # load velocity for this month
        vxo = np.array(topaz.variables['vxo'][t_idx, :, lat_s:lat_e, :])
        vyo = np.array(topaz.variables['vyo'][t_idx, :, lat_s:lat_e, :])

        # depth-weighted layer mean velocities
        vx1, vy1 = layer_mean_velocity(vxo, vyo, grid['layer1_idxs'], dz)
        vx2, vy2 = layer_mean_velocity(vxo, vyo, grid['layer2_idxs'], dz)

        S = get_riverine_source(grid, year, month, riv_df, DT)
        # atmospheric deposition omitted -- net air-sea flux small relative to
        # oceanic inflow; MacInnis 2017 Devon Ice Cap record not representative
        # of open ocean deposition. Contribution subsumed into PINN correction.
        # + get_atmospheric_source(grid, year, atm_df, DT)
        # time integration
        for _ in range(steps_per_month):
            C1 = upwind_step(C1, vx1, vy1, grid['dx'], grid['dy'],
                             grid['land_mask'], DT, c_south_bnd=c_bnd_L1)
            C2 = upwind_step(C2, vx2, vy2, grid['dx'], grid['dy'],
                             grid['land_mask'], DT, c_south_bnd=c_bnd_L2)
            C1 = C1 + S
            C1[grid['land_mask']] = 0.0

            # first-order decay sink (lumped: sediment burial, degradation,
            # sea ice partitioning). PINN will learn residual correction.
            C1 *= DECAY_FACTOR
            C2 *= DECAY_FACTOR

            C1, C2 = apply_diffusion(C1, C2, grid, DT)

        # write monthly snapshot
        out.variables['C_layer1'][t_idx] = C1
        out.variables['C_layer2'][t_idx] = C2

        si = stations['i'].values
        sj = stations['j'].values
        out.variables['C_station_L1'][t_idx] = C1[si, sj]
        out.variables['C_station_L2'][t_idx] = C2[si, sj]
        out.sync()

    topaz.close()
    out.close()
    print(f'simulation complete. output: {OUTPUT_PATH}')


if __name__ == '__main__':
    run()
