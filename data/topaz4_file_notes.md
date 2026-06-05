# TOPAZ4 NetCDF File Notes
# File: topaz4_arctic_velocity_2004_2025.nc
# Inspected: 2026-06-05

## File identity
- Product: ARCTIC_MULTIYEAR_PHY_002_003
- Dataset: cmems_mod_arc_phy_my_topaz4_P1M (monthly means)
- Source: NERSC-HYCOM model fields
- Institution: NERSC, Jahnebakken 3, N-5007 Bergen, Norway
- DOI: https://doi.org/10.48670/moi-00007
- Format: NETCDF4 / HDF5
- Title: Arctic Ocean Physics Reanalysis, 12.5km monthly mean
- field_date: 2024-12
- copernicusmarine_version: 2.3.0

## Dimensions
- depth:     40
- latitude:  321
- longitude: 2880
- time:      264

## Coordinates
- lat:   50.000 to 90.000 N, step 0.1250 deg (~13.9 km)
- lon:  -180.000 to 179.875 E, step 0.1250 deg (~3.5 km at 75N)
- time: 2004-01-01 to 2025-12-01, monthly, hours since 1950-01-01
- time calendar: proleptic_gregorian
- depth levels (m):
    [0, 2, 4, 6, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100,
     125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900,
     1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000]

## Variables
All 4D variables have shape (time, depth, latitude, longitude).
All use _FillValue: nan. No missing_value attribute set.
All have regrid_method: bilinear -- regridded from native polar stereographic
to regular lat/lon grid. Velocity components are therefore in geographic
east/north directions; no rotation required.

| Name        | Dims                        | Units      | Description                        |
|-------------|-----------------------------|------------|------------------------------------|
| vxo         | (time, depth, lat, lon)     | m/s        | Eastward sea water velocity        |
| vyo         | (time, depth, lat, lon)     | m/s        | Northward sea water velocity       |
| thetao      | (time, depth, lat, lon)     | degrees_C  | Sea water potential temperature    |
| so          | (time, depth, lat, lon)     | 1e-3       | Sea water salinity                 |
| mlotst      | (time, lat, lon)            | m          | Mixed layer depth (sigma-theta)    |
| zos         | (time, lat, lon)            | m          | Sea surface height above geoid     |
| model_depth | (lat, lon)                  | m          | Sea floor depth below sea level    |

## Land mask
- model_depth is NaN over land (466725 of 924480 cells = ~50.5%)
- model_depth min/max over ocean: 10.47 to 4968.67 m
- No zero-depth ocean cells -- land identified purely by NaN mask
- vxo/vyo are also NaN over land, consistent with model_depth mask

## Velocity statistics (t=0, z=0, surface January 2004)
- vxo min/max/mean: -0.4502 / 0.3335 / -0.005992 m/s
- vyo min/max/mean: -0.4560 / 0.3147 / +0.000679 m/s
- Max speed ~0.45 m/s. At 0.125 deg (~3.5 km at 75N), 12-hour timestep
  gives CFL ~ 0.45 * 43200 / 3500 ~ 5.6 at surface max.
  NOTE: need to verify typical vs peak velocities -- 0.45 m/s may be
  a localized extreme (Fram Strait, coastlines). Domain-wide typical
  CFL at 12hr step is well under 1.

## Two-layer scheme
- Layer 1 (surface):  depth indices 0-19,  0 to 200 m
- Layer 2 (deep):     depth indices 20-39, 250 to 4000 m
- Split at depth index 19 = exactly 200 m

## Simulation domain slice
- Use latitude indices 80:321 (60.0 to 90.0 N, 241 points)
- Full longitude range (2880 points)
- Lat index 80 = 60.0 N: (50.0 + 80 * 0.125 = 60.0) confirmed

## Notes
- thetao, so, zos not needed for simulation but present for PINN feature use later
- mlotst useful for inter-layer diffusion parameterization
- File is a single monolithic NetCDF4; load variables as memory-mapped slices,
  do not load full 4D arrays at once (~39 GB per 4D variable)