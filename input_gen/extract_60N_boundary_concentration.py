import pandas as pd
import numpy as np
import os

# Paths relative to repo root
DATA_PATH = os.path.join("experiment", "PFAS In-Situ Observations - Sheet1.csv")
OUTPUT_PATH = os.path.join("experiment", "boundary_60N_pfoa_surface.csv")

# Known regional-mean DOIs — these are valid data but not independent station obs.
# Li 2018 (10.1016/j.envpol.2018.03.018): regional min/max/mean entered per region,
# not per station. Per-station data not published. Treat as n=1 regional estimate
# per region, not as independent observations.
REGIONAL_MEAN_DOIS = {
    "10.1016/j.envpol.2018.03.018": "Li 2018 — regional summary stats, not per-station"
}

# Load
df = pd.read_csv(DATA_PATH)

# Filter: PFOA only
df = df[df["compound"].str.upper() == "PFOA"]

# Filter: above 200m
df = df[df["depth_m"] <= 200]

# Filter: exclude below detection
df = df[df["below_detection"] != True]
df = df[df["below_detection"] != "TRUE"]

# Flag regional mean rows
def classify_doi(doi):
    for key in REGIONAL_MEAN_DOIS:
        if key in str(doi):
            return "regional_mean"
    return "station_obs"

df["data_type"] = df["source_doi"].apply(classify_doi)

# Filter: lat band 55-65N
df_band = df[(df["lat"] >= 55) & (df["lat"] <= 65)].copy()

print(f"Total PFOA obs above 200m, non-BDL: {len(df)}")
print(f"Obs in 55-65N band: {len(df_band)}")
print()

# Split by data type
df_station = df_band[df_band["data_type"] == "station_obs"]
df_regional = df_band[df_band["data_type"] == "regional_mean"]

# --- Station observations ---
print("=== STATION OBSERVATIONS (55-65N) ===")
if len(df_station) > 0:
    print(df_station[["lat", "lon", "depth_m", "PFAS_concentration_ngL", "source_doi"]].to_string())
    print()

    atl = df_station[(df_station["lon"] >= -80) & (df_station["lon"] <= 40)]
    pac = df_station[(df_station["lon"] > 120) | (df_station["lon"] < -120)]

    if len(atl) > 0:
        print(f"Atlantic (lon -80 to 40):    n={len(atl)}, mean={atl['PFAS_concentration_ngL'].mean():.4f} ng/L, median={atl['PFAS_concentration_ngL'].median():.4f} ng/L")
    else:
        print("Atlantic: no station obs")

    if len(pac) > 0:
        print(f"Pacific  (lon >120 or <-120): n={len(pac)}, mean={pac['PFAS_concentration_ngL'].mean():.4f} ng/L, median={pac['PFAS_concentration_ngL'].median():.4f} ng/L")
    else:
        print("Pacific: no station obs")
else:
    print("No station observations in 55-65N band.")

print()

# --- Regional mean entries ---
print("=== REGIONAL MEAN ENTRIES (55-65N) — not independent observations ===")
if len(df_regional) > 0:
    for doi, note in REGIONAL_MEAN_DOIS.items():
        subset = df_regional[df_regional["source_doi"].str.contains(doi, na=False)]
        if len(subset) > 0:
            print(f"  {note}")
            print(f"  Rows: {len(subset)}, values: {subset['PFAS_concentration_ngL'].unique().tolist()} ng/L")
            print(f"  Treat as single regional estimate per distinct value, not n={len(subset)}")
            print()
else:
    print("None.")

print()

# --- Recommended boundary values ---
print("=== RECOMMENDED BOUNDARY VALUES ===")
if len(atl) > 0:
    print(f"Atlantic 60N surface layer: {atl['PFAS_concentration_ngL'].mean():.4f} ng/L (mean of {len(atl)} station obs)")
else:
    print("Atlantic 60N surface layer: no station obs — source manually")

# For Pacific, use Li 2018 BS regional mean (0.0826 ng/L) as single estimate
li2018_pac = df_regional[(df_regional["lon"] > 120) | (df_regional["lon"] < -120)]
if len(li2018_pac) > 0:
    pac_val = li2018_pac["PFAS_concentration_ngL"].unique()
    print(f"Pacific 60N surface layer:  {pac_val[0]:.4f} ng/L (Li 2018 Bering Sea regional mean, n=1 regional estimate)")

# Save full band output with data_type column
df_band.to_csv(OUTPUT_PATH, index=False)
print(f"\nSaved to {OUTPUT_PATH}")