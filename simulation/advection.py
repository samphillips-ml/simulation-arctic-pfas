"""
advection.py
Upwind finite-difference advection step for Arctic PFOA simulation.
Operates on a single 2D horizontal concentration slice (one depth level).
Caller loops over depth levels and timesteps.
"""

import numpy as np


def upwind_step(C, vx, vy, dx, dy, land_mask, dt, c_south_bnd=None):
    """
    One upwind advection timestep for a 2D concentration field.

    Upwind scheme: for each cell, the upstream gradient is used based on
    the sign of the velocity. Periodic boundary in longitude. Fixed
    boundaries at 60N (south) and 90N (north).

    Parameters
    ----------
    C           : (nlat, nlon) float
                  Concentration field in ng/L.
    vx          : (nlat, nlon) float
                  Eastward velocity in m/s. NaN over land.
    vy          : (nlat, nlon) float
                  Northward velocity in m/s. NaN over land.
    dx          : (nlat, nlon) float
                  Zonal cell width in metres.
    dy          : float
                  Meridional cell height in metres (constant).
    land_mask   : (nlat, nlon) bool
                  True = land cell. Land cells are zeroed after each step.
    dt          : float
                  Timestep in seconds.
    c_south_bnd : (nlon,) float or None
                  Prescribed concentration at the 60N southern boundary
                  in ng/L. Applied to inflow cells (vy > 0 at i=0).
                  If None, zero-gradient is assumed (no inflow from south).

    Returns
    -------
    C_new : (nlat, nlon) float
            Updated concentration field in ng/L.
    """
    # --- zero velocity over land (no flux through land) ---
    vx = np.where(np.isnan(vx), 0.0, vx)
    vy = np.where(np.isnan(vy), 0.0, vy)

    # --- guard against dx = 0 at the pole ---
    dx_safe = np.maximum(dx, 1.0)
    # --- cap velocity to CFL <= 1 in every cell ---
    # near 90N dx is tiny (30m at 89.875N), giving extreme CFL at any
    # non-trivial velocity; also caps Bering Strait and other high-speed cells
    vx = np.sign(vx) * np.minimum(np.abs(vx), dx_safe / dt)
    vy = np.sign(vy) * np.minimum(np.abs(vy), dy     / dt)

    # --- shifted concentration arrays ---

    # longitude: periodic (wraps at -180/180)
    C_west = np.roll(C, 1, axis=1)     # C_west[i,j]  = C[i, j-1]
    C_east = np.roll(C, -1, axis=1)    # C_east[i,j]  = C[i, j+1]

    # latitude: non-periodic
    # interior: C_south[i,j] = C[i-1,j],  C_north[i,j] = C[i+1,j]
    C_south = np.empty_like(C)
    C_north = np.empty_like(C)

    C_south[1:, :] = C[:-1, :]
    C_north[:-1, :] = C[1:, :]

    # southern boundary (i=0, 60N)
    # inflow (vy > 0, northward): use prescribed boundary concentration
    # outflow (vy <= 0, southward): zero-gradient
    if c_south_bnd is not None:
        inflow = vy[0, :] > 0.0
        C_south[0, :] = np.where(inflow, c_south_bnd, C[0, :])
    else:
        C_south[0, :] = C[0, :]        # zero-gradient

    # northern boundary (i=-1, 90N): zero-gradient
    C_north[-1, :] = C[-1, :]

    # --- upwind differences ---
    # vx >= 0 (eastward): tracer arrives from west, backward difference
    # vx <  0 (westward): tracer arrives from east, forward difference
    dCdx = np.where(vx >= 0.0,
                    (C - C_west) / dx_safe,
                    (C_east - C) / dx_safe)

    # vy >= 0 (northward): tracer arrives from south, backward difference
    # vy <  0 (southward): tracer arrives from north, forward difference
    dCdy = np.where(vy >= 0.0,
                    (C - C_south) / dy,
                    (C_north - C) / dy)

    # --- forward Euler update ---
    C_new = C - dt * (vx * dCdx + vy * dCdy)

    # --- enforce non-negativity ---
    # CFL > 1 in rare coastal cells can produce small negative values;
    # clip rather than let them accumulate
    np.maximum(C_new, 0.0, out=C_new)

    # --- zero land cells ---
    C_new[land_mask] = 0.0

    return C_new
