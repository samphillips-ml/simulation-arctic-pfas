"""
visualize_observations.py

Diagnostics and visuals for PFOA in-situ observations above 60N.
  1. Histogram of PFOA concentrations (detected only)
  2. Concentration vs depth (all detected stations)
  3. Concentration vs latitude
  4. Time series by year (box plots)
  5. Print summary statistics

Run from repo root:
  python visualize_observations.py

Output: output/observation_diagnostics.png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

OBS_PATH = "data/PFAS In-Situ Observations - Sheet1.csv"
OUT_PATH = "output/observation_diagnostics.png"
os.makedirs("output", exist_ok=True)

# -- load and filter -----------------------------------------------------------
obs  = pd.read_csv(OBS_PATH)
pfoa = obs[
    (obs["compound"] == "PFOA") &
    (obs["below_detection"].astype(str).str.upper() != "TRUE") &
    (obs["lat"] >= 60.0)
].copy()

pfoa["PFAS_concentration_ngL"] = pd.to_numeric(
    pfoa["PFAS_concentration_ngL"], errors="coerce")
pfoa = pfoa.dropna(subset=["PFAS_concentration_ngL"])

# parse year from date column
pfoa["year"] = pd.to_numeric(
    pfoa["data"].astype(str).str[:4], errors="coerce")

print(f"PFOA detected observations above 60N: {len(pfoa)}")
print(f"Unique stations: {pfoa[['lat','lon']].drop_duplicates().shape[0]}")
print(f"Year range: {pfoa['year'].min():.0f} - {pfoa['year'].max():.0f}")
print(f"\nConcentration summary (ng/L):")
print(pfoa["PFAS_concentration_ngL"].describe().round(4))
print(f"\nPercentiles:")
for p in [5, 10, 25, 50, 75, 90, 95, 99]:
    v = np.percentile(pfoa["PFAS_concentration_ngL"], p)
    print(f"  p{p:2d}: {v:.4f} ng/L")

# -- plot ----------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("PFOA In-Situ Observations -- Detected, Above 60N",
             fontsize=13, fontweight="bold")

# panel 1: histogram
ax = axes[0, 0]
vals = pfoa["PFAS_concentration_ngL"]
ax.hist(vals, bins=40, color="steelblue", edgecolor="white", linewidth=0.3)
ax.axvline(vals.median(), color="red", lw=1.2, linestyle="--",
           label=f"median={vals.median():.3f}")
ax.axvline(vals.mean(), color="orange", lw=1.2, linestyle="--",
           label=f"mean={vals.mean():.3f}")
ax.set_xlabel("PFOA concentration (ng/L)")
ax.set_ylabel("Count")
ax.set_title("Concentration distribution")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# panel 2: concentration vs depth
ax = axes[0, 1]
depth = pd.to_numeric(pfoa["depth_m"], errors="coerce")
ax.scatter(vals, depth, alpha=0.4, s=15, color="steelblue")
ax.invert_yaxis()
ax.set_xlabel("PFOA concentration (ng/L)")
ax.set_ylabel("Depth (m)")
ax.set_title("Concentration vs depth")
ax.grid(True, alpha=0.3)

# panel 3: concentration vs latitude
ax = axes[1, 0]
ax.scatter(pfoa["lat"], vals, alpha=0.4, s=15, color="darkorange")
ax.set_xlabel("Latitude (deg N)")
ax.set_ylabel("PFOA concentration (ng/L)")
ax.set_title("Concentration vs latitude")
ax.grid(True, alpha=0.3)

# panel 4: box plot by year
ax = axes[1, 1]
years = sorted(pfoa["year"].dropna().unique())
data_by_year = [pfoa[pfoa["year"] == y]["PFAS_concentration_ngL"].values
                for y in years]
bp = ax.boxplot(data_by_year, patch_artist=True,
                medianprops=dict(color="red", lw=1.5))
for patch in bp["boxes"]:
    patch.set_facecolor("steelblue")
    patch.set_alpha(0.6)
ax.set_xticks(range(1, len(years)+1))
ax.set_xticklabels([str(int(y)) for y in years], rotation=45, fontsize=8)
ax.set_xlabel("Year")
ax.set_ylabel("PFOA concentration (ng/L)")
ax.set_title("Concentration by year")
ax.grid(True, alpha=0.3, axis="y")

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
print(f"\nSaved: {OUT_PATH}")