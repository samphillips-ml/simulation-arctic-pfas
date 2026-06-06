"""
sources.py
Source term computation for Arctic PFOA simulation.
Provides riverine injection, atmospheric deposition, Bering Strait
inflow, and 60N Atlantic boundary concentration fields. All
concentrations in ng/L.
"""

import numpy as np
import pandas as pd

# ------------------------------------------------------------------
# river mouth coordinates
# nearest ocean cell to approximate delta/bay outlet, verified against
# TOPAZ4 model_depth mask (all coordinates confirmed as ocean cells).
#
# All ten coordinates checked against the TOPAZ4 land mask via
# check_dvina_mask.py / check_new_river_mouths.py. Where the natural
# delta/bay coordinate fell on land (common at 0.125 deg resolution for
# narrow deltas), the nearest wet cell was used instead. All shifts are
# well within the BOX_HALF=16 (~2 deg) injection box half-width, so the
# injected signal still lands in the correct general shelf area.
#
# Special case -- NorthernDvina: TOPAZ4 masks the entire White Sea
# basin (no wet cells found south of ~70.5N in 32-42E). Flux is injected
# at the southernmost available Barents Sea cell, representing the
# White Sea / Barents Sea exchange zone. Implied transport delay through
# the White Sea is absorbed by bias correction. See DECISIONS.md sec. 9.
# ------------------------------------------------------------------
RIVER_MOUTHS = {
    'Lena':          (74.000,  126.500),   # Laptev Sea outlet
    'Yenisei':       (73.750,   82.000),   # Kara Sea outlet
    'Ob':            (72.500,   73.625),   # Kara Sea outlet
    'Mackenzie':     (70.000, -134.000),   # Beaufort Sea outlet
    'Kolyma':        (70.000,  162.125),   # Kolyma Gulf, East Siberian Sea
    'Pechora':       (69.125,   54.000),   # Pechora Sea, Barents
    'Yana':          (71.750,  136.500),   # Yana Bay, Laptev Sea
    'Indigirka':     (71.750,  151.750),   # East Siberian Sea, off delta
    'Olenek':        (73.500,  122.500),   # Olenek Bay, Laptev Sea
    'NorthernDvina': (70.500,   38.250),   # White Sea masked; inject S. Barents
}

# 60N boundary concentration (ng/L), applied to 0-200m inflow cells.
# Atlantic sector only -- the Pacific sector at 60N is now land
# (Bering Sea masked, see grid_utils.py), so the Pacific branch that
# used to live here has been replaced by the Bering Strait source term
# below. See DECISIONS.md sec. 10.
# Atlantic sector mean: Benskin 2012 (Labrador Sea) + Joerss 2020 (Norwegian coast)
BOUNDARY_ATLANTIC_NGPL = 0.032

# ------------------------------------------------------------------
# Bering Strait inflow concentration
#
# Central value 0.05 ng/L (50 pg/L), converging estimate from two
# independent PFOA-specific surface-seawater datasets for the
# Pacific-inflow / western Arctic water mass:
#   - Cai et al. 2012 (Environ. Sci. Technol. 46(2):661-668,
#     doi:10.1021/es2026278): North Pacific PFOA average 56 pg/L,
#     range <20-100 pg/L.
#   - Yamazaki et al. 2021 (Chemosphere 272:129869,
#     doi:10.1016/j.chemosphere.2020.128803): western Arctic Ocean
#     PFOA 48-87 pg/L.
# Sensitivity bounds span the low end of Cai (<20 pg/L) and the high
# end of Yamazaki (87 pg/L), rounded to 0.02 and 0.09 ng/L.
# See DECISIONS.md sec. 10 for full derivation and caveats (no
# strait-specific measurement exists; both sources are regional
# proxies for the water mass transiting the strait; values are
# 2010-2014-era and may slightly overestimate current PFOA given the
# US EPA PFOA Stewardship Program phase-down).
# ------------------------------------------------------------------
BERING_PFOA_NGPL          = 0.05    # central value, ng/L
BERING_PFOA_NGPL_LOW      = 0.02    # sensitivity bound, ng/L
BERING_PFOA_NGPL_HIGH     = 0.09    # sensitivity bound, ng/L

# Bering Strait throat cells -- narrow lat/lon band spanning TOPAZ4's
# actual resolved Pacific-side opening.
#
# NOTE: the real-world strait (~65.7N, 168.5-169.5W) is NOT open water
# in TOPAZ4 at 0.125 deg -- diagnostic check (diag_bering.py) against
# the raw model_depth mask found zero wet cells anywhere in the
# -172 to -165W band south of 66.375N. The model's coastline geometry
# pushes the effective channel north of the true strait, likely because
# the real ~85 km wide, Diomede-Island-obstructed channel is narrower
# than this grid resolves. This box targets TOPAZ4's actual opening
# (first wet row at 66.375N, widening northward) rather than the
# real-world coordinate. See DECISIONS.md sec. 10 for the full
# discussion and the implied limitation (the source term now
# represents Pacific inflow somewhat north of the literal strait,
# not the strait itself).
# Must sit north of grid_utils.BERING_SEA_LAT_MAX (64.5) so the cells
# survive land-masking.
BERING_STRAIT_LAT_RANGE = (66.375, 67.0)
BERING_STRAIT_LON_RANGE = (-172.0, -166.625)


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


def _bering_strait_cells(grid):
    """
    Return (i_idx, j_idx) arrays of wet-cell indices within the Bering
    Strait throat box. Computed fresh each call rather than cached --
    cheap (small box) and avoids stale-index bugs if grid changes.
    """
    lat_mask = ((grid['lat'] >= BERING_STRAIT_LAT_RANGE[0]) &
                (grid['lat'] <= BERING_STRAIT_LAT_RANGE[1]))
    lon_mask = ((grid['lon'] >= BERING_STRAIT_LON_RANGE[0]) &
                (grid['lon'] <= BERING_STRAIT_LON_RANGE[1]))
    box = lat_mask[:, np.newaxis] & lon_mask[np.newaxis, :]
    wet = box & (~grid['land_mask'])
    i_idx, j_idx = np.where(wet)
    return i_idx, j_idx


def load_riverine(csv_path):
    """Load riverine_pfoa_flux.csv. Returns DataFrame."""
    return pd.read_csv(csv_path)


def load_atmospheric(csv_path):
    """Load atmospheric_pfoa_deposition.csv. Returns DataFrame."""
    return pd.read_csv(csv_path)


def get_riverine_source(grid, year, month, riverine_df, dt):
    """
    Riverine source field for one timestep.

    Flux is distributed evenly across ocean cells in a box centred on
    the nearest ocean cell to each river delta/bay outlet (see
    RIVER_MOUTHS for placement notes). Distributing over a box rather
    than a single cell prevents unrealistic concentration spikes at the
    mouth; physically represents rapid delta mixing. Layer 1 depth
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

    # ng/m^2 -> ng/L: divide by surface layer depth (m) and unit conversion
    delta_C = flux_ng_m2_step / (dz_surface * 1000.0)
    S[~grid['land_mask']] = delta_C
    return S


def get_bering_strait_source(grid, vx, vy, dt, c_strait=BERING_PFOA_NGPL):
    """
    Bering Strait inflow source field for one timestep.

    Replaces the old diffuse 60N Pacific boundary. The Bering Sea south
    of the strait is now land-masked (grid_utils.py), so the strait
    throat cells are the sole Pacific-side opening into the domain.

    Mechanism: injects mass at strait throat cells sized as
        concentration x inflow velocity x cross-sectional area x dt
    using the live TOPAZ4 velocity field for that month, restricted to
    cells with net northward (into-domain) flow. This is a source-term
    approximation rather than a true clamped open boundary (the
    upwind_step southern-boundary mechanism in advection.py only
    operates at the fixed domain edge, lat index 0 / 60N, and cannot
    prescribe inflow concentration at an interior latitude without
    restructuring the advection scheme -- see DECISIONS.md sec. 10 for
    the tradeoff discussion). Outflow cells (vy <= 0, water leaving the
    domain southward through the strait) are left untouched; advection
    handles outflow via the normal upwind scheme.

    Parameters
    ----------
    grid     : dict from grid_utils.load_grid()
    vx, vy   : (nlat, nlon) float, layer 1 depth-weighted mean velocity
               in m/s (NaN over land). Same arrays simulate.py already
               computes for the advection step this timestep.
    dt       : float, timestep in seconds
    c_strait : float, prescribed inflow concentration in ng/L
               (default BERING_PFOA_NGPL; pass BERING_PFOA_NGPL_LOW/HIGH
               for sensitivity runs)

    Returns
    -------
    S : (nlat, nlon) ng/L to add to layer 1 concentration field.
    """
    S = np.zeros((grid['nlat'], grid['nlon']))
    dz_layer1 = 200.0   # layer 1 thickness (m), consistent with rivers

    i_idx, j_idx = _bering_strait_cells(grid)
    if len(i_idx) == 0:
        # strait throat fully masked -- almost certainly a mask bounds
        # bug (see grid_utils.py __main__ self-check), not a real state
        return S

    vy_strait = vy[i_idx, j_idx]
    vy_strait = np.where(np.isnan(vy_strait), 0.0, vy_strait)
    inflow = vy_strait > 0.0   # northward = into the domain

    if not np.any(inflow):
        return S

    ii = i_idx[inflow]
    jj = j_idx[inflow]
    v_in = vy_strait[inflow]   # m/s, positive northward

    # volume flux through each cell's northern face this step (L)
    cell_dx = grid['dx'][ii, jj]
    vol_L_step = v_in * cell_dx * dz_layer1 * dt * 1000.0

    # mass injected (ng), then converted to a concentration bump in
    # that cell's own volume (same convention as get_riverine_source:
    # cell_vol_L is the receiving cell's volume, not the flux volume)
    mass_ng = c_strait * vol_L_step   # ng/L x L = ng

    cell_vol_L = cell_dx * grid['dy'] * dz_layer1 * 1000.0
    S[ii, jj] += mass_ng / cell_vol_L

    return S


def get_boundary_1d(grid, depth_m):
    """
    Prescribed 60N southern boundary concentration for one depth level.

    Atlantic sector only (-80 to 40E): 0.032 ng/L
    Everywhere else at 60N: 0.0 ng/L -- the Pacific sector is now land
    (Bering Sea masked south of the strait), so there is no Pacific
    inflow at the 60N domain edge anymore. Pacific inflow enters via
    get_bering_strait_source() instead. See DECISIONS.md sec. 10.
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
    c_bnd[atlantic] = BOUNDARY_ATLANTIC_NGPL
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

    dt = 6 * 3600.0   # 6-hour timestep, diagnostic only (actual sim dt
                       # is passed in from simulate.py)

    # --- riverine ---
    S_riv = get_riverine_source(grid, 2004, 7, riv_df, dt)
    mouth_idx = _mouth_indices(grid)
    print('Riverine source (ng/L per 6hr step) at mouth cells, July 2004:')
    for river, (i, j) in mouth_idx.items():
        wet = not grid['land_mask'][i, j]
        print(f'  {river:<14} {S_riv[i, j]:.4f}  '
              f'lat={grid["lat"][i]:.3f} lon={grid["lon"][j]:.3f}  '
              f'wet={wet}')

    # --- atmospheric ---
    S_atm = get_atmospheric_source(grid, 2004, atm_df, dt)
    ocean_vals = S_atm[~grid['land_mask']]
    print(f'\nAtmospheric source 2004 (ng/L per 6hr step):')
    print(f'  uniform ocean value: {ocean_vals[0]:.6e}')
    print(f'  land cells: {(S_atm[grid["land_mask"]] == 0).all()} (all zero)')

    # --- Bering Strait ---
    i_idx, j_idx = _bering_strait_cells(grid)
    print(f'\nBering Strait throat: {len(i_idx)} wet cells in box '
          f'lat={BERING_STRAIT_LAT_RANGE} lon={BERING_STRAIT_LON_RANGE}')
    if len(i_idx) > 0:
        print(f'  lat range of wet cells: '
              f'{grid["lat"][i_idx].min():.3f} - {grid["lat"][i_idx].max():.3f}')
        print(f'  lon range of wet cells: '
              f'{grid["lon"][j_idx].min():.3f} - {grid["lon"][j_idx].max():.3f}')
        # diagnostic source call with a plausible northward velocity
        # to confirm the function produces sane magnitudes
        vx_test = np.full((grid['nlat'], grid['nlon']), np.nan)
        vy_test = np.full((grid['nlat'], grid['nlon']), np.nan)
        vy_test[i_idx, j_idx] = 0.3   # m/s, plausible strait inflow speed
        S_bering = get_bering_strait_source(grid, vx_test, vy_test, dt)
        nz = S_bering[i_idx, j_idx]
        print(f'  diagnostic S at 0.3 m/s inflow (ng/L per 6hr step): '
              f'min={nz.min():.6f} max={nz.max():.6f} mean={nz.mean():.6f}')
    else:
        print('  WARNING: no wet cells found -- check BERING_STRAIT_LAT_RANGE '
              'against grid_utils.BERING_SEA_LAT_MAX, the mask may have '
              'sealed the strait shut')

    # --- boundary ---
    print('\n60N boundary concentration (ng/L):')
    for depth_m in [0.0, 100.0, 200.0, 250.0]:
        c = get_boundary_1d(grid, depth_m)
        atl = c[(grid['lon'] >= -80) & (grid['lon'] <= 40)]
        print(f'  depth={depth_m}m  atlantic={atl[atl>0].mean() if len(atl[atl>0]) else 0:.3f}'
              f'  (pacific sector now land -- see Bering Strait source above)')
