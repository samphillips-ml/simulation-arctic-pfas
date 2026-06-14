"""
visualize_simulation.py

Animated GIF of Arctic PFOA layer 1 concentration in North Polar
Stereographic projection. Colorscale spans the full range of
positive concentrations (no floor -- atmospheric deposition is not
represented in this simulation, so there is no "background" value
to clip against).

Usage:
    python simulation/visualize_simulation.py
    python simulation/visualize_simulation.py output/simulation_pfoa.nc output/pfoa_animation.gif
"""

import sys
import numpy as np
import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.animation import FuncAnimation, PillowWriter
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# -- paths ---------------------------------------------------------------------
OUTPUT_NC = sys.argv[1] if len(sys.argv) > 1 else 'output/simulation_pfoa.nc'
GIF_PATH  = sys.argv[2] if len(sys.argv) > 2 else 'output/pfoa_layer1.gif'
FPS       = 4   # slower for polar stereo -- rendering is heavier

# -- load ----------------------------------------------------------------------
print(f'loading {OUTPUT_NC}...')
ds     = netCDF4.Dataset(OUTPUT_NC, 'r')
C1_raw = np.array(ds.variables['C_layer1'][:])   # (n_months, nlat, nlon)
lat    = np.array(ds.variables['lat'][:])
lon    = np.array(ds.variables['lon'][:])
ds.close()

n_months, nlat, nlon = C1_raw.shape
print(f'  {n_months} months, {nlat} lat x {nlon} lon')

# -- prepare -------------------------------------------------------------------
C1 = C1_raw.astype(float)
C1[~np.isfinite(C1)] = np.nan
C1[C1 <= 0.0] = np.nan

# colorscale spans the full positive range of the data (no floor):
# vmin = smallest positive concentration, vmax = 99th percentile
# (99th percentile still used on the upper end to avoid blowup cells
# at river-mouth injection points dominating the scale)
finite   = C1[np.isfinite(C1)]
vmin_lin = float(np.nanmin(finite)) if len(finite) > 0 else 1e-6
vmax_lin = float(np.nanpercentile(finite, 99)) if len(finite) > 0 else 1.0

# use log10 colorscale
vmin = np.log10(vmin_lin)
vmax = np.log10(vmax_lin)
print(f'  colorscale: {vmin_lin:.3e} to {vmax_lin:.3f} ng/L  (log10: {vmin:.1f} to {vmax:.1f})')

C1_log = np.log10(C1)

# 2D meshgrid for pcolormesh
dlat = float(lat[1] - lat[0])
dlon = float(lon[1] - lon[0])
lat_edges = np.append(lat - dlat/2, lat[-1] + dlat/2)
lon_edges = np.append(lon - dlon/2, lon[-1] + dlon/2)
lon_e2d, lat_e2d = np.meshgrid(lon_edges, lat_edges)

# domain stats for right panel
domain_max  = np.nanmax(C1,     axis=(1, 2))
domain_mean = np.nanmean(C1,    axis=(1, 2))

# -- figure setup --------------------------------------------------------------
data_crs = ccrs.PlateCarree()
proj     = ccrs.NorthPolarStereo()

fig = plt.figure(figsize=(14, 7))
ax_map  = fig.add_subplot(1, 2, 1, projection=proj)
ax_stat = fig.add_subplot(1, 2, 2)

ax_map.set_extent([-180, 180, 55, 90], crs=data_crs)

# static land from cartopy (faster than re-drawing TOPAZ mask each frame)
ax_map.add_feature(cfeature.LAND, facecolor='#cccccc', zorder=2)
ax_map.add_feature(cfeature.COASTLINE, linewidth=0.4, zorder=3)
ax_map.gridlines(linewidth=0.3, color='gray', alpha=0.5)

# initial pcolormesh -- will update data each frame
mesh = ax_map.pcolormesh(
    lon_e2d, lat_e2d, C1_log[0],
    transform=data_crs,
    cmap='plasma',
    vmin=vmin, vmax=vmax,
    shading='flat',
    zorder=1,
)

# colorbar
tick_vals = np.arange(np.ceil(vmin), np.floor(vmax) + 1)
cbar = fig.colorbar(mesh, ax=ax_map, fraction=0.04, pad=0.04,
                    shrink=0.7, ticks=tick_vals)
cbar.set_label('PFOA (ng/L)', fontsize=10)
cbar.ax.set_yticklabels([f'$10^{{{int(v)}}}$' for v in tick_vals])

title = ax_map.set_title('', fontsize=11, fontweight='bold', pad=8)

# right panel: time series
months_x = np.arange(n_months)
ax_stat.plot(months_x, domain_max,  color='crimson',  lw=1.5, label='max')
ax_stat.plot(months_x, domain_mean, color='steelblue', lw=1.5, label='mean')
ax_stat.axhline(0.05, color='green', lw=1.0, linestyle='--',
                label='obs median (0.05)')
ax_stat.axhline(0.17, color='orange', lw=1.0, linestyle='--',
                label='obs max (0.17)')
ax_stat.set_xlabel('Month since 2004-01')
ax_stat.set_ylabel('PFOA (ng/L)')
ax_stat.set_title('Domain statistics', fontsize=10)
ax_stat.legend(fontsize=7)
ax_stat.set_xlim(0, n_months - 1)
ax_stat.grid(True, alpha=0.3)
ax_stat.set_yscale('log')
# year tick marks
year_ticks = [i * 12 for i in range(n_months // 12 + 1)]
ax_stat.set_xticks(year_ticks)
ax_stat.set_xticklabels([str(2004 + i) for i in range(len(year_ticks))],
                         rotation=45, fontsize=7)

vline = ax_stat.axvline(x=0, color='black', lw=1.5, zorder=5)

fig.tight_layout()

# -- animation -----------------------------------------------------------------
def update(frame):
    yr = 2004 + frame // 12
    mo = 1    + frame  % 12
    mesh.set_array(C1_log[frame].ravel())
    title.set_text(f'PFOA Layer 1 (0-200m)   {yr}-{mo:02d}')
    vline.set_xdata([frame, frame])
    return mesh, title, vline

print(f'rendering {n_months} frames at {FPS} fps...')
print('  (polar stereo is slow -- expect several minutes)')
ani = FuncAnimation(fig, update, frames=n_months,
                    interval=1000 // FPS, blit=True)
ani.save(GIF_PATH, writer=PillowWriter(fps=FPS), dpi=100)
plt.close()
print(f'saved: {GIF_PATH}')