"""
sources.py
Source term computation for Arctic PFOA simulation.
Provides riverine injection, atmospheric deposition, and 60N boundary
concentration fields. All concentrations in ng/L.
"""

import numpy as np
import pandas as pd

# ------------------------------------------------------------------
# river mouth coordinates
# nearest ocean cell to approximate delta outlet, verified against
# TOPAZ4 model_depth mask (all coordinates confirmed as ocean cells)
# ------------------------------------------------------------------
RIVER_MOUTHS = {
    'Lena':      (74.000,  126.500),   # Laptev Sea outlet
    'Yenisei':   (73.750,   82.000),   # Kara Sea outlet
    'Ob':        (72.500,   73.625),   # Kara Sea outlet
    'Mackenzie': (70.000, -134.000),   # Beaufort Sea outlet
}

# 60N boundary concentrations (ng/L), applied to 0-200m inflow cells
# Atlantic sector mean: Benskin 2012 (Labrador Sea) + Joerss 2020 (Norwegian coast)
# Pacific sector: Li 2018 Bering Sea regional mean (n=1 regional estimate)
BOUNDARY_ATLANTIC_NGPL = 0.032
BOUNDARY_PACIFIC_NGPL  = 0.025


def _compute_dz(depth):
    """
    Layer thicknesses in metres from TOPAZ4 depth level centres.
    thickness[i] = distance between midpoints of adjacent levels.
    """
    depth = np.asarray(depth, dtype=float)
    mids  = 0.5 * (depth[:-1] + depth[1:])
    dz        = np.empty(len(depth))
    dz[0]     = mids[0]
    dz[1:-1]  = mids[1:] - mids[:-1]
    dz[-1]    = depth[-1] - mids[-1]
    return dz


def _mouth_indices(grid):
    """Return dict of river -> (i, j) indices in the 60-90N domain."""
    out = {}
    for river, (rlat, rlon) in RIVER_MOUTHS.items():
        i = int(np.argmin(np.abs(grid['lat'] - rlat)))
        j = int(np.argmin(np.abs(grid['lon'] - rlon)))
        out[river] = (i, j)
    return out


def load_riverine(csv_path):
    """Load riverine_pfoa_flux.csv. Returns DataFrame."""
    return pd.read_csv(csv_path)


def load_atmospheric(csv_path):
    """Load atmospheric_pfoa_deposition.csv. Returns DataFrame."""
    return pd.read_csv(csv_path)


def get_riverine_source(grid, year, month, riverine_df, dt):
    """
    Riverine source field for one timestep.

    Flux is distributed evenly across ocean cells in a 5x5 box centred on
    the nearest ocean cell to each river delta outlet. Distributing over a
    box rather than a single cell prevents unrealistic concentration spikes
    at the mouth; physically represents rapid delta mixing. Layer 1 depth
    (200m) used as injection volume so source is diluted over the full
    surface layer rather than the ~1m surface skin.

    Parameters
    ----------
    grid        : dict from grid_utils.load_grid()
    year, month : int
    riverine_df : DataFrame from load_riverine()
    dt          : float, timestep in seconds

    Returns
    -------
    S : (nlat, nlon) ng/L to add to layer 1 concentration field.
    """
    S  = np.zeros((grid['nlat'], grid['nlon']))
    dz_layer1  = 200.0                         # layer 1 thickness (m)
    BOX_HALF   = 16                            # 33x33 injection box half-width

    seconds_per_month = (365.25 / 12.0) * 86400.0
    mouth_idx = _mouth_indices(grid)

    rows = riverine_df[
        (riverine_df['year'] == year) &
        (riverine_df['month'] == month)
    ]

    for _, row in rows.iterrows():
        river = row['river']
        if river not in mouth_idx:
            continue

        i, j = mouth_idx[river]
        flux_kg_month = row['pfoa_flux_kg_month']
        if pd.isna(flux_kg_month):
            continue

        # mass injected this timestep (ng)
        flux_ng_step = flux_kg_month * 1e12 * (dt / seconds_per_month)

        # collect ocean cells in box
        box_cells = []
        for di in range(-BOX_HALF, BOX_HALF + 1):
            for dj in range(-BOX_HALF, BOX_HALF + 1):
                ii = i + di
                jj = (j + dj) % grid['nlon']   # periodic longitude
                if 0 <= ii < grid['nlat'] and not grid['land_mask'][ii, jj]:
                    box_cells.append((ii, jj))

        if len(box_cells) == 0:
            continue

        # distribute flux equally across ocean cells in box
        for ii, jj in box_cells:
            cell_vol_L = grid['dx'][ii, jj] * grid['dy'] * dz_layer1 * 1000.0
            S[ii, jj] += flux_ng_step / cell_vol_L / len(box_cells)

    return S


def get_atmospheric_source(grid, year, atmos_df, dt):
    """
    Atmospheric deposition source field for one timestep.

    Applied uniformly to all ocean surface cells. Flux from MacInnis et al.
    2017 (1993-2007); post-2007 held constant at 2004-2007 mean per Yeung
    et al. 2017.

    Parameters
    ----------
    grid      : dict from grid_utils.load_grid()
    year      : int
    atmos_df  : DataFrame from load_atmospheric()
    dt        : float, timestep in seconds

    Returns
    -------
    S : (nlat, nlon) ng/L to add to concentration at depth index 0.
    """
    S = np.zeros((grid['nlat'], grid['nlon']))
    dz = _compute_dz(grid['depth'])
    dz_surface = dz[0]

    row = atmos_df[atmos_df['year'] == year]
    if len(row) == 0:
        return S
    flux_ng_m2_yr = float(row.iloc[0]['pfoa_deposition_ng_m2_yr'])

    seconds_per_year  = 365.25 * 86400.0
    flux_ng_m2_step   = flux_ng_m2_yr * (dt / seconds_per_year)

    # ng/m² -> ng/L: divide by surface layer depth (m) and unit conversion
    delta_C = flux_ng_m2_step / (dz_surface * 1000.0)

    S[~grid['land_mask']] = delta_C
    return S


def get_boundary_1d(grid, depth_m):
    """
    Prescribed 60N southern boundary concentration for one depth level.

    Atlantic sector (-80 to 40E): 0.093 ng/L
    Pacific sector  (>120E or <-120W): 0.083 ng/L
    Transition zones: 0.0 ng/L (conservative)
    Below 200m: 0.0 ng/L (Yeung et al. 2017)

    Parameters
    ----------
    grid    : dict from grid_utils.load_grid()
    depth_m : float, depth of this level in metres

    Returns
    -------
    c_bnd : (nlon,) ng/L
    """
    lon   = grid['lon']
    c_bnd = np.zeros(len(lon))

    if depth_m > 200.0:
        return c_bnd

    atlantic = (lon >= -80.0) & (lon <= 40.0)
    pacific  = (lon >  120.0) | (lon < -120.0)

    c_bnd[atlantic] = BOUNDARY_ATLANTIC_NGPL
    c_bnd[pacific]  = BOUNDARY_PACIFIC_NGPL

    return c_bnd


if __name__ == '__main__':
    import os
    from grid_utils import load_grid

    data_dir  = os.path.join(os.path.dirname(__file__), '..', 'data')
    topaz4    = os.path.join(data_dir, 'topaz4_arctic_velocity_2004_2025.nc')
    riv_csv   = os.path.join(data_dir, 'riverine_pfoa_flux.csv')
    atm_csv   = os.path.join(data_dir, 'atmospheric_pfoa_deposition.csv')

    grid = load_grid(topaz4)
    riv_df = load_riverine(riv_csv)
    atm_df = load_atmospheric(atm_csv)

    dt = 6 * 3600.0   # 6-hour timestep

    # --- riverine ---
    S_riv = get_riverine_source(grid, 2004, 7, riv_df, dt)
    mouth_idx = _mouth_indices(grid)
    print('Riverine source (ng/L per 6hr step) at mouth cells, July 2004:')
    for river, (i, j) in mouth_idx.items():
        print(f'  {river}: {S_riv[i, j]:.4f}  lat={grid["lat"][i]:.3f} lon={grid["lon"][j]:.3f}')

    # --- atmospheric ---
    S_atm = get_atmospheric_source(grid, 2004, atm_df, dt)
    ocean_vals = S_atm[~grid['land_mask']]
    print(f'\nAtmospheric source 2004 (ng/L per 6hr step):')
    print(f'  uniform ocean value: {ocean_vals[0]:.6e}')
    print(f'  land cells: {(S_atm[grid["land_mask"]] == 0).all()} (all zero)')

    # --- boundary ---
    print('\n60N boundary concentration (ng/L):')
    for depth_m in [0.0, 100.0, 200.0, 250.0]:
        c = get_boundary_1d(grid, depth_m)
        atl = c[(grid['lon'] >= -80) & (grid['lon'] <= 40)]
        pac = c[(grid['lon'] > 120) | (grid['lon'] < -120)]
        print(f'  depth={depth_m}m  atlantic={atl[atl>0].mean() if len(atl[atl>0]) else 0:.3f}'
              f'  pacific={pac[pac>0].mean() if len(pac[pac>0]) else 0:.3f}')