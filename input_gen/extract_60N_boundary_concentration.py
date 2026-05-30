import pandas as pd
import numpy as np
import os

# Paths relative to repo root
DATA_PATH = os.path.join("experiment", "PFAS In-Situ Observations - Sheet1.csv")
OUTPUT_PATH = os.path.join("experiment", "boundary_60N_pfoa_surface.csv")

# Load
df = pd.read_csv(DATA_PATH)

# Filter: PFOA only
df = df[df["compound"].str.upper() == "PFOA"]

# Filter: above 200m
df = df[df["depth_m"] <= 200]

# Filter: exclude below detection
df = df[df["below_detection"] != True]
df = df[df["below_detection"] != "TRUE"]

# Filter: lat band 55-65N to capture 60N boundary
df_band = df[(df["lat"] >= 55) & (df["lat"] <= 65)].copy()

print(f"Total PFOA obs above 200m, non-BDL: {len(df)}")
print(f"Obs in 55-65N band: {len(df_band)}")
print()

if len(df_band) > 0:
    print("--- 55-65N band summary ---")
    print(df_band[["lat", "lon", "depth_m", "PFAS_concentration_ngL", "source_doi"]].to_string())
    print()
    print(f"Mean concentration:   {df_band['PFAS_concentration_ngL'].mean():.4f} ng/L")
    print(f"Median concentration: {df_band['PFAS_concentration_ngL'].median():.4f} ng/L")
    print(f"Std:                  {df_band['PFAS_concentration_ngL'].std():.4f} ng/L")
    print(f"N:                    {len(df_band)}")
    print()

    # Split Atlantic vs Pacific by longitude
    atl = df_band[(df_band["lon"] >= -80) & (df_band["lon"] <= 40)]
    pac = df_band[(df_band["lon"] > 120) | (df_band["lon"] < -120)]

    if len(atl) > 0:
        print(f"Atlantic sector (lon -80 to 40):      n={len(atl)}, mean={atl['PFAS_concentration_ngL'].mean():.4f} ng/L, median={atl['PFAS_concentration_ngL'].median():.4f} ng/L")
    else:
        print("Atlantic sector: no obs")

    if len(pac) > 0:
        print(f"Pacific sector (lon >120 or <-120):   n={len(pac)}, mean={pac['PFAS_concentration_ngL'].mean():.4f} ng/L, median={pac['PFAS_concentration_ngL'].median():.4f} ng/L")
    else:
        print("Pacific sector: no obs")

    df_band.to_csv(OUTPUT_PATH, index=False)
    print(f"\nSaved filtered obs to {OUTPUT_PATH}")

else:
    print("No observations in 55-65N band. Widening to 50-70N...")
    df_wide = df[(df["lat"] >= 50) & (df["lat"] <= 70)].copy()
    print(f"Obs in 50-70N band: {len(df_wide)}")

    if len(df_wide) > 0:
        print(df_wide[["lat", "lon", "depth_m", "PFAS_concentration_ngL", "source_doi"]].to_string())
        print()
        print(f"Mean concentration:   {df_wide['PFAS_concentration_ngL'].mean():.4f} ng/L")
        print(f"Median concentration: {df_wide['PFAS_concentration_ngL'].median():.4f} ng/L")
        print(f"Std:                  {df_wide['PFAS_concentration_ngL'].std():.4f} ng/L")
        print(f"N:                    {len(df_wide)}")

        atl = df_wide[(df_wide["lon"] >= -80) & (df_wide["lon"] <= 40)]
        pac = df_wide[(df_wide["lon"] > 120) | (df_wide["lon"] < -120)]

        if len(atl) > 0:
            print(f"\nAtlantic sector: n={len(atl)}, mean={atl['PFAS_concentration_ngL'].mean():.4f} ng/L, median={atl['PFAS_concentration_ngL'].median():.4f} ng/L")
        else:
            print("Atlantic sector: no obs")

        if len(pac) > 0:
            print(f"Pacific sector:  n={len(pac)}, mean={pac['PFAS_concentration_ngL'].mean():.4f} ng/L, median={pac['PFAS_concentration_ngL'].median():.4f} ng/L")
        else:
            print("Pacific sector: no obs")

        df_wide.to_csv(OUTPUT_PATH, index=False)
        print(f"\nSaved filtered obs to {OUTPUT_PATH}")
    else:
        print("No observations found even in 50-70N band. Check data.")
