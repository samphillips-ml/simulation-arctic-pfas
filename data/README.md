On the Orion Cluster, the .nc file is called topaz4_arctic_velocity_2004_2025.nc and lives in this directory.

Info about topaz4_arctic_velocity_2004_2025.nc
## topaz4_arctic_velocity_2004_2025.nc

CMEMS Arctic Ocean Physics Reanalysis (TOPAZ4b), 12.5km monthly mean.
Dataset ID: cmems_mod_arc_phy_my_topaz4_P1M
DOI: https://doi.org/10.48670/moi-00007

Time: 2004-01-01 to 2025-12-01, 264 monthly steps
Lat: 50-90N (321 points), Lon: -180 to 180 (2880 points), Depth: 0-4000m (40 levels)

Variables:
- vxo, vyo       eastward/northward velocity (m/s), 4D, 39GB each
- thetao         potential temperature (degrees_C), 4D, 39GB
- so             salinity (PSU), 4D, 39GB
- mlotst         mixed layer depth (m), 3D, 976MB
- zos            sea surface height (m), 3D, 976MB
- model_depth    bathymetry (m), 2D, 4MB

Velocity is already rotated to geographic east/north by CMEMS at download -- no further rotation needed.
Always open lazy and slice to 65N+, target depth range before loading anything into memory.

## \[River\] Discharge
These files are gotten from ArcticGRO and contain the historical and current discharge data. they are used to calculate the flux values.