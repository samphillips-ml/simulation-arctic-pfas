"""
grid_utils.py
Grid metadata for Arctic PFOA simulation.
Loads lat/lon/depth from TOPAZ4 NetCDF, computes cell dimensions,
builds land mask, and defines two-layer depth split.
"""

import numpy as np
import netCDF4

# ------------------------------------------------------------------
# constants
# ------------------------------------------------------------------
EARTH_RADIUS_M   = 6371000.0
LAT_MIN          = 60.0        # southern boundary of simulation domain
LAT_MAX          = 90.0
LAYER_SPLIT_M    = 200.0       # layer 1: 0-200m, layer 2: 250-4000m

# Bering Sea land-mask box -- everything south of the Bering Strait
# throat is masked to land, leaving the strait itself as the sole
# Pacific-side opening into the domain. Box bounds chosen to fully
# remove the Bering Sea (St. Lawrence Island, Norton Sound, Gulf of
# Anadyr, Anadyr Strait) while leaving the strait throat open.
#
# IMPORTANT: TOPAZ4 at 0.125 deg does not resolve the real-world strait
# location (~65.7N) as an open channel -- diagnostic check
# (diag_bering.py) against the raw model_depth mask found zero wet
# cells in the -172 to -165W band at any latitude south of 66.375N.
# The first open row appears at 66.375N (7 wet cells, -170.125 to
# -169.0). TOPAZ4's coastline geometry pushes the effective Pacific
# opening north of the true strait location, presumably because the
# real channel (~85 km wide, obstructed by the Diomede Islands) is too
# narrow/shallow to resolve at this grid spacing. BERING_SEA_LAT_MAX
# is therefore set well south of the actual TOPAZ4 opening (64.5N,
# comfortable margin below 66.375N) rather than at the real-world
# strait latitude, to avoid masking the only channel TOPAZ4 actually
# provides. See DECISIONS.md sec. 10 and sources.py
# BERING_STRAIT_LAT_RANGE for the corresponding throat-cell box.
# Longitude range expressed in the -180..180 convention used by the
# TOPAZ4 grid; wraps are not needed here since the Bering Sea does not
# cross the antimeridian in this convention (it spans roughly 163E to
# -165W, i.e. positive lon > 163 OR negative lon > -165, see
# BERING_SEA_LON_RANGES below).
BERING_SEA_LAT_MAX        = 64.5   # mask everything south of this
BERING_SEA_LON_RANGES     = [(163.0, 180.0), (-180.0, -165.0)]  # deg E


def _bering_sea_mask(lat, lon):
    """
    Boolean (nlat, nlon) mask, True where a cell falls within the
    Bering Sea box south of the Bering Strait throat.
    """
    lat_mask = lat < BERING_SEA_LAT_MAX                  # (nlat,)

    lon_mask = np.zeros(len(lon), dtype=bool)
    for lo, hi in BERING_SEA_LON_RANGES:
        lon_mask |= (lon >= lo) & (lon <= hi)            # (nlon,)

    return lat_mask[:, np.newaxis] & lon_mask[np.newaxis, :]


def load_grid(topaz4_path):
    """
    Load grid metadata from TOPAZ4 NetCDF file.

    Parameters
    ----------
    topaz4_path : str
        Path to topaz4_arctic_velocity_2004_2025.nc

    Returns
    -------
    dict with keys:
        lat           : (nlat,)      degrees N, 60-90N
        lon           : (nlon,)      degrees E, -180 to 179.875
        depth         : (40,)        depth level centres in metres
        dx            : (nlat, nlon) zonal cell width in metres
        dy            : float        meridional cell height in metres (constant)
        land_mask     : (nlat, nlon) bool, True = land/no-data cell
        model_depth   : (nlat, nlon) float, seafloor depth in metres (NaN = land)
        lat_idx_min   : int          index into full TOPAZ4 lat array for LAT_MIN
        lat_idx_max   : int          index into full TOPAZ4 lat array for LAT_MAX + 1
        nlat          : int          241
        nlon          : int          2880
        layer_split_idx : int        first depth index of layer 2 (index 20, 250m)
        layer1_idxs   : slice        depth indices for layer 1 (0:20, 0-200m)
        layer2_idxs   : slice        depth indices for layer 2 (20:, 250-4000m)
    """
    ds = netCDF4.Dataset(topaz4_path, 'r')

    # --- coordinates ---
    lat_full = ds.variables['latitude'][:]    # (321,) 50-90N
    lon      = ds.variables['longitude'][:]   # (2880,) -180 to 179.875
    depth    = ds.variables['depth'][:]       # (40,)

    # slice to simulation domain (60-90N)
    lat_idx_min = int(np.searchsorted(lat_full, LAT_MIN))
    lat_idx_max = int(np.searchsorted(lat_full, LAT_MAX)) + 1
    lat = lat_full[lat_idx_min:lat_idx_max]   # (241,)

    nlat = len(lat)
    nlon = len(lon)

    # --- cell dimensions in metres ---
    # regular grid so dy is constant; dx shrinks with cos(lat)
    dlat_rad = np.radians(float(lat[1] - lat[0]))
    dlon_rad = np.radians(float(lon[1] - lon[0]))

    dy = EARTH_RADIUS_M * dlat_rad            # ~13900 m, scalar

    dx_1d = EARTH_RADIUS_M * np.cos(np.radians(lat)) * dlon_rad  # (nlat,)
    dx    = np.broadcast_to(dx_1d[:, np.newaxis], (nlat, nlon)).copy()  # (nlat, nlon)
    # note: dx -> 0 at 90N (cos 90 = 0). the north pole is land in TOPAZ4
    # and will be masked, but advection.py should guard against dx == 0.

    # --- land mask ---
    md_full    = np.array(ds.variables['model_depth'][:])          # (321, 2880)
    model_depth = md_full[lat_idx_min:lat_idx_max, :]              # (241, 2880)
    land_mask  = np.isnan(model_depth)                             # True = land

    # mask polar cap above 89N -- tiny dx cells cause numerical instability
    # on regular lat/lon grid; no observations exist above 89N in training data
    POLAR_CAP_MASK_LAT = 89.0
    land_mask[lat[:, np.newaxis].repeat(nlon, axis=1) >= POLAR_CAP_MASK_LAT] = True

    # mask Bering Sea south of the Bering Strait throat -- the strait
    # becomes the sole Pacific-side opening into the domain, replacing
    # the diffuse 60N Pacific boundary. See DECISIONS.md sec. 10.
    land_mask[_bering_sea_mask(lat, lon)] = True

    # --- layer split ---
    # depth[19] = 200.0 m exactly; layer 1 includes 200m, layer 2 starts at 250m
    layer_split_idx = int(np.searchsorted(depth, LAYER_SPLIT_M))
    assert float(depth[layer_split_idx]) == LAYER_SPLIT_M, (
        f"Layer split mismatch: depth[{layer_split_idx}] = {depth[layer_split_idx]}, "
        f"expected {LAYER_SPLIT_M}"
    )
    layer1_idxs = slice(0, layer_split_idx + 1)    # 0:20  depths 0-200m
    layer2_idxs = slice(layer_split_idx + 1, None) # 20:   depths 250-4000m

    ds.close()

    return {
        'lat':             lat,
        'lon':             lon,
        'depth':           depth,
        'dx':              dx,
        'dy':              dy,
        'land_mask':       land_mask,
        'model_depth':     model_depth,
        'lat_idx_min':     lat_idx_min,
        'lat_idx_max':     lat_idx_max,
        'nlat':            nlat,
        'nlon':            nlon,
        'layer_split_idx': layer_split_idx,
        'layer1_idxs':     layer1_idxs,
        'layer2_idxs':     layer2_idxs,
    }


if __name__ == '__main__':
    import sys
    import os

    if len(sys.argv) >= 2:
        topaz4_path = sys.argv[1]
    else:
        topaz4_path = os.path.join(
            os.path.dirname(__file__), '..', 'data',
            'topaz4_arctic_velocity_2004_2025.nc'
        )

    print(f'loading grid from {topaz4_path}')
    g = load_grid(topaz4_path)

    print(f'lat:         {g["lat"][0]:.3f} to {g["lat"][-1]:.3f} N  ({g["nlat"]} points)')
    print(f'lon:         {g["lon"][0]:.3f} to {g["lon"][-1]:.3f} E  ({g["nlon"]} points)')
    print(f'dy:          {g["dy"]:.1f} m')
    print(f'dx at 60N:   {g["dx"][0, 0]:.1f} m')
    print(f'dx at 75N:   {g["dx"][120, 0]:.1f} m')
    print(f'dx at 89N:   {g["dx"][-2, 0]:.1f} m')
    print(f'land cells:  {g["land_mask"].sum()} of {g["land_mask"].size} '
          f'({100 * g["land_mask"].mean():.1f}%)')
    print(f'layer 1:     depth indices 0:{g["layer_split_idx"] + 1}  '
          f'(0 - {g["depth"][g["layer_split_idx"]]:.0f} m)')
    print(f'layer 2:     depth indices {g["layer_split_idx"] + 1}:   '
          f'({g["depth"][g["layer_split_idx"] + 1]:.0f} - '
          f'{g["depth"][-1]:.0f} m)')

    # sanity check: confirm strait throat remains open. NOTE: TOPAZ4
    # does not resolve the real-world strait latitude (~65.7N) as open
    # water -- diag_bering.py found the first open row at 66.375N. This
    # check tests the actual TOPAZ4 opening, not the real-world location.
    strait_lat_idx = np.argmin(np.abs(g['lat'] - 66.5))
    strait_lon_mask = (g['lon'] >= -171.5) & (g['lon'] <= -167.0)
    strait_wet = (~g['land_mask'][strait_lat_idx, strait_lon_mask]).sum()
    print(f'\nBering Strait throat check (lat idx {strait_lat_idx}, '
          f'lat={g["lat"][strait_lat_idx]:.3f}N):')
    print(f'  wet cells in 171.5-167.0W band: {strait_wet}  '
          f'(expect > 0 -- if 0, BERING_SEA_LAT_MAX is masking the strait itself)')