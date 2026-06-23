"""
Phase B (lives in validation/ - this is the real test: is the simulation
right, not just stable). Collapses the master observation CSV to PFOA-only
physical sampling events, matches each to the corresponding simulation
output (exact month, no smoothing window - per decision), and computes
residual = obs - sim.

Run from repo root:

    python validation/compute_residuals.py

Decisions baked into this script (document in DECISIONS.md):
  - Compound filter: PFOA only. Events with no PFOA measurement among
    their compound-level rows are dropped entirely (not kept as
    NaN/zero), since "no PFOA measured" is not the same as "PFOA absent."
  - BDL handling: PFOA rows with below_detection=TRUE are EXCLUDED from
    the residual calc. Detection-limit-as-value is appropriate for some
    purposes but would bias residuals toward whatever the lab's MDL
    happened to be, not the true (unknown) concentration.
  - Comparison convention: EXACT matching month only. An obs dated
    2009-08 matches sim step index = (2009-2004)*12 + (8-1) = 67.
    No averaging window.
  - Depth -> layer mapping: depth_m < 200 -> C_layer1, depth_m >= 200
    -> C_layer2. Events exactly at 200m go to layer 2 (boundary
    convention, adjust if you want the opposite).
  - Held-out Sept 2025 validation samples: NOT YET IN THE MASTER CSV
    as of this writing. HOLDOUT_FILTER below is a stub - once those
    samples are added, fill in the actual exclusion logic (by date,
    source_doi, or a dedicated flag column) so they can never
    accidentally leak into this residual computation.
"""

import numpy as np
import pandas as pd
import xarray as xr

# --- config -----------------------------------------------------------
OBS_CSV_PATH = "data/master_pfas_observations.csv"  # ADJUST to actual path
NC_PATH = "output/FINAL_simulation_pfoa-20_06_2026.nc"
OUT_PATH = "validation/obs_vs_sim_residuals.csv"

DATE_COL = "data"  # yes, really - known bug, kept intentionally (refactor later)
DEPTH_BOUNDARY_M = 200  # depth_m >= this -> layer 2

SIM_START_YEAR = 2004  # step 0 = Jan 2004


def sim_step_index(date_str):
    """'YYYY-MM' or 'YYYY-MM-DD' -> integer step index into the 264-step
    sim output. Per project data conventions, both formats appear in the
    master CSV; day-of-month (if present) is ignored since the simulation
    is monthly-resolution."""
    parts = str(date_str).split("-")
    year, month = parts[0], parts[1]
    return (int(year) - SIM_START_YEAR) * 12 + (int(month) - 1)


def holdout_filter(df):
    """
    STUB. The 21 Sept 2025 held-out validation samples are not yet in
    the master CSV. Once they're added, replace this with the actual
    exclusion logic - e.g.:

        df = df[df["source_doi"] != "<2025 expedition doi/identifier>"]

    or by date:

        df = df[~df[DATE_COL].astype(str).str.startswith("2025-09")]

    Until then this is a no-op, but it's called explicitly below so
    that filling it in later is a one-line change, not a retrofit.
    """
    return df


def main():
    # --- load and filter observations ---
    obs = pd.read_csv(OBS_CSV_PATH)
    print(f"Loaded {len(obs)} compound-level rows from {OBS_CSV_PATH}")

    obs = holdout_filter(obs)

    pfoa = obs[obs["compound"] == "PFOA"].copy()
    print(f"Rows with compound == 'PFOA': {len(pfoa)}")

    n_bdl = int(pfoa["below_detection"].astype(str).str.upper().eq("TRUE").sum())
    pfoa = pfoa[~pfoa["below_detection"].astype(str).str.upper().eq("TRUE")]
    print(f"Excluded {n_bdl} below-detection PFOA rows; {len(pfoa)} remain")

    # Collapse to physical events. Since we already filtered to PFOA only,
    # each (lat, lon, depth_m, date) group should now be a single row in
    # almost all cases - but group explicitly anyway in case of true
    # duplicate PFOA measurements at the same event (e.g. replicate
    # samples), in which case we average them.
    group_cols = ["lat", "lon", "depth_m", DATE_COL]
    before = len(pfoa)
    events = (
        pfoa.groupby(group_cols, as_index=False)
        .agg(
            PFAS_concentration_ngL=("PFAS_concentration_ngL", "mean"),
            n_replicates=("PFAS_concentration_ngL", "count"),
            region=("region", "first"),
            source_doi=("source_doi", "first"),
        )
    )
    n_dupes = (events["n_replicates"] > 1).sum()
    print(f"Collapsed {before} PFOA rows into {len(events)} physical events "
          f"({n_dupes} events had >1 replicate, averaged)")

    # --- domain filter ---
    # Hard cutoff at the simulation grid's actual bounds (60-90N, confirmed
    # via check_grid_coords.py). Without this, nearest-neighbor matching
    # below will silently match far-outside-domain stations (e.g. East
    # China Sea, Narragansett Bay, Bay of Biscay - all present in the raw
    # multi-paper master CSV) to the nearest Arctic grid cell, producing
    # huge spurious residuals that reflect "this station isn't in the
    # simulated domain," not "the model is wrong." This is NOT optional.
    GRID_LAT_MIN = 60.0
    n_before_domain = len(events)
    events = events[events["lat"] >= GRID_LAT_MIN].copy()
    n_dropped_domain = n_before_domain - len(events)
    print(f"Domain filter (lat >= {GRID_LAT_MIN}N): dropped {n_dropped_domain} "
          f"out-of-domain events, {len(events)} remain")

    # --- compute sim step index per event, drop out-of-range ---
    n_dashes = events[DATE_COL].astype(str).str.count("-")
    print(f"Date format check: {(n_dashes == 1).sum()} events as YYYY-MM, "
          f"{(n_dashes == 2).sum()} events as YYYY-MM-DD")

    events["sim_step"] = events[DATE_COL].apply(sim_step_index)
    in_range = (events["sim_step"] >= 0) & (events["sim_step"] < 264)
    n_out = (~in_range).sum()
    if n_out > 0:
        print(f"WARNING: {n_out} events fall outside the 2004-2025 sim window - dropping")
    events = events[in_range].copy()

    # --- depth -> layer mapping ---
    events["layer"] = np.where(events["depth_m"] < DEPTH_BOUNDARY_M, 1, 2)

    # --- open sim output and load grid coords (needed for both the
    # land-mask filter below and the sim-matching step after it) ---
    ds = xr.open_dataset(NC_PATH, decode_times=False)
    lat_vals = ds["lat"].values
    lon_vals = ds["lon"].values

    def nearest_idx(arr, val):
        return int(np.abs(arr - val).argmin())

    # --- land-mask filter ---
    # Confirmed via check_greenland_cells.py: 2 events landed on TOPAZ4
    # land cells (NaN model_depth) near Greenland's coast - same category
    # of grid-resolution limit as the Bering Strait throat (see grid_utils.py
    # docstring, DECISIONS.md sec. 10). Not a simulation bug; these obs
    # events simply can't be represented at this grid resolution. Filtering
    # by actual land_mask (rather than hardcoding those 2 coordinates)
    # generalizes to any future obs event landing on a land cell.
    import sys as _sys
    _sys.path.insert(0, "simulation")
    import grid_utils as _grid_utils

    _grid = _grid_utils.load_grid("data/topaz4_arctic_velocity_2004_2025.nc")
    _land_mask = _grid["land_mask"]

    def _is_land(lat_val, lon_val):
        i_lat = nearest_idx(lat_vals, lat_val)
        i_lon = nearest_idx(lon_vals, lon_val)
        return bool(_land_mask[i_lat, i_lon])

    n_before_land = len(events)
    events["matched_cell_is_land"] = events.apply(
        lambda r: _is_land(r["lat"], r["lon"]), axis=1
    )
    events_before_land_filter = events.copy()
    events = events[~events["matched_cell_is_land"]].copy()
    n_dropped_land = n_before_land - len(events)
    print(f"Land-mask filter: dropped {n_dropped_land} events matching a TOPAZ4 "
          f"land cell (coastline/fjord resolution limit), {len(events)} remain")
    if n_dropped_land > 0:
        dropped = events_before_land_filter[events_before_land_filter["matched_cell_is_land"]]
        print("\nDropped events (lat, lon, region, source_doi):")
        print(dropped[["lat", "lon", "region", "source_doi"]].to_string(index=False))
        print()

    # --- pull matching sim values ---
    sim_concs = []
    match_dist_deg = []
    for _, row in events.iterrows():
        i_lat = nearest_idx(lat_vals, row["lat"])
        i_lon = nearest_idx(lon_vals, row["lon"])
        var = "C_layer1" if row["layer"] == 1 else "C_layer2"
        val = ds[var].isel(time=int(row["sim_step"]), lat=i_lat, lon=i_lon).values.item()
        sim_concs.append(val)
        match_dist_deg.append(
            np.hypot(lat_vals[i_lat] - row["lat"], lon_vals[i_lon] - row["lon"])
        )

    events["sim_concentration"] = sim_concs
    events["match_dist_deg"] = match_dist_deg
    events["residual"] = events["PFAS_concentration_ngL"] - events["sim_concentration"]

    # Flag (don't auto-drop) matches that are suspiciously far from the
    # actual observation - grid spacing is 0.125deg, so anything much
    # larger than that suggests a coastline/land-mask mismatch worth a
    # manual look, not necessarily a reason to exclude.
    DIST_FLAG_THRESHOLD_DEG = 0.5
    n_flagged = (events["match_dist_deg"] > DIST_FLAG_THRESHOLD_DEG).sum()
    if n_flagged > 0:
        print(f"NOTE: {n_flagged} events matched a grid cell >{DIST_FLAG_THRESHOLD_DEG} deg "
              f"away from the actual obs coordinate - worth a manual look (coastline/land-mask?)")

    ds.close()

    # --- save ---
    events.to_csv(OUT_PATH, index=False)
    print(f"\nSaved {len(events)} obs-vs-sim residual rows to {OUT_PATH}")

    # --- quick look ---
    print("\nResidual summary:")
    print(events["residual"].describe())
    print(f"\nMean residual (obs - sim): {events['residual'].mean():.6f}")
    print(f"Median residual:           {events['residual'].median():.6f}")
    print("(positive = simulation underestimates; matches the documented")
    print(" ~2-2.5x underestimate pattern if mean/median are positive)")


if __name__ == "__main__":
    main()