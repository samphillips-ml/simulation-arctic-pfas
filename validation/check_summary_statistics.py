"""
Phase A sanity check: per-timestep summary stats for the frozen simulation
output. Run from repo root:

    python validation/check_summary_stats.py

Checks, per the Phase A plan:
  1. No-blowup / no-NaN: min, max, mean, NaN count per timestep
  2. Magnitude/trend: basin-mean concentration trend across the run

Does NOT touch TOPAZ4 input. Reads only the frozen output .nc.
"""

import numpy as np
import xarray as xr

# --- config -----------------------------------------------------------
NC_PATH = "output/FINAL_simulation_pfoa-20_06_2026.nc"

# Update this if the variable name in your output differs.
CONC_VAR = "pfoa_concentration"

# Expected basin-mean trend bounds (ng/L), loose sanity bounds based on
# prior runs: ~0.005 (2005) to ~0.040 (early 2020s), allowing some slack.
EARLY_EXPECTED_RANGE = (0.001, 0.015)   # first few timesteps
LATE_EXPECTED_RANGE = (0.015, 0.080)    # last few timesteps
# ------------------------------------------------------------------------


def main():
    ds = xr.open_dataset(NC_PATH)

    if CONC_VAR not in ds:
        print(f"ERROR: variable '{CONC_VAR}' not found in dataset.")
        print(f"Available variables: {list(ds.data_vars)}")
        return

    da = ds[CONC_VAR]

    # Figure out the time dimension name (commonly 'time' or 'month')
    time_dim = None
    for cand in ("time", "month", "t"):
        if cand in da.dims:
            time_dim = cand
            break
    if time_dim is None:
        print(f"ERROR: could not identify time dimension. Dims: {da.dims}")
        return

    n_steps = da.sizes[time_dim]
    print(f"Found {n_steps} timesteps along dim '{time_dim}'.\n")

    # --- per-timestep stats ---
    print(f"{'step':>5} {'min':>12} {'max':>12} {'mean':>12} {'nan_count':>10}")
    print("-" * 56)

    means = []
    any_nan = False
    any_negative = False

    for i in range(n_steps):
        frame = da.isel({time_dim: i}).values
        nan_count = int(np.isnan(frame).sum())
        valid = frame[~np.isnan(frame)]

        if nan_count > 0:
            any_nan = True

        if valid.size == 0:
            print(f"{i:>5} {'ALL NAN':>12} {'-':>12} {'-':>12} {nan_count:>10}")
            means.append(np.nan)
            continue

        vmin, vmax, vmean = valid.min(), valid.max(), valid.mean()
        if vmin < 0:
            any_negative = True
        means.append(vmean)

        # Only print every step if small; otherwise print first/last 5 + flags
        flag = ""
        if nan_count > 0:
            flag += " <-- NaNs present"
        if vmin < 0:
            flag += " <-- NEGATIVE VALUE"

        print(f"{i:>5} {vmin:>12.6f} {vmax:>12.6f} {vmean:>12.6f} {nan_count:>10}{flag}")

    means = np.array(means)

    # --- summary ---
    print("\n" + "=" * 56)
    print("SUMMARY")
    print("=" * 56)
    print(f"Total timesteps:        {n_steps}")
    print(f"Any NaNs present:       {any_nan}")
    print(f"Any negative values:    {any_negative}  (should be False - upwind scheme guarantees non-negativity)")

    early = means[:5]
    late = means[-5:]
    print(f"\nBasin mean, first 5 steps: {early}")
    print(f"Basin mean, last 5 steps:  {late}")

    early_mean = np.nanmean(early)
    late_mean = np.nanmean(late)
    print(f"\nEarly-period mean: {early_mean:.6f} ng/L  (expected range: {EARLY_EXPECTED_RANGE})")
    print(f"Late-period mean:  {late_mean:.6f} ng/L  (expected range: {LATE_EXPECTED_RANGE})")

    early_ok = EARLY_EXPECTED_RANGE[0] <= early_mean <= EARLY_EXPECTED_RANGE[1]
    late_ok = LATE_EXPECTED_RANGE[0] <= late_mean <= LATE_EXPECTED_RANGE[1]

    print(f"\nEarly period within expected bounds: {early_ok}")
    print(f"Late period within expected bounds:  {late_ok}")

    if not (early_ok and late_ok):
        print("\n*** Trend is outside expected bounds - investigate before freezing this run. ***")

    if any_nan or any_negative:
        print("\n*** NaNs or negative values present - investigate before freezing this run. ***")

    ds.close()


if __name__ == "__main__":
    main()