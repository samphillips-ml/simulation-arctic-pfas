"""
gen_riverine_flux.py
Generate riverine_pfoa_flux.csv for Arctic PFOA simulation.

Ten rivers, all discharge data from:
  The Arctic Great Rivers Observatory. 2026. Discharge Dataset,
  Version 20260210. https://arcticgreatrivers.org/data

Total annual PFOA flux: 14,800 kg/yr (Stemmler & Lammel 2010,
doi:10.5194/acp-10-9965-2010). Partitioned by mean annual discharge
fraction (2004-2025), shaped seasonally assuming chemostatic behavior
(flux proportional to discharge; PFOA concentration uniform year-round).

Northern Dvina injected at 70.5N, 38.25E (southernmost Barents cell)
because TOPAZ4 masks the White Sea entirely. See DECISIONS.md sec. 9.

Gap handling (see DECISIONS.md sec. 9):
  Yana:       2000-2017 gap (204 months). Filled using calendar-month
              climatology (mean discharge for that month across all
              other available years for this river). Linear time-
              interpolation was tested and rejected: Yana has extreme
              seasonality (near-zero in winter, >40,000 m3/s at the
              summer freshet peak), and a straight line across a 17-year
              gap applies a near-constant intermediate value to every
              month including winter, roughly doubling the true annual
              mean (1826 vs ~1071 m3/s from the raw daily record).
  Indigirka:  no data 2004-2015. Zero-filled (not climatology-filled) --
              this is a genuine absence of injection during this period,
              not a data gap to be patched. Injection begins 2016 when
              real discharge data resumes.
  Olenek:     same as Indigirka.
  All other rivers with minor gaps: filled using the same calendar-month
  climatology method for consistency.

Required input files in data/:
  Lena-Discharge.xlsx                                  ArcticGRO, Lena at Kyusyur
  Yenisei-Discharge.xlsx                               ArcticGRO, Yenisei at Igarka
  Ob-Discharge.xlsx                                    ArcticGRO, Ob at Salekhard
  Mackenzie-Discharge.xlsx                             ArcticGRO, Mackenzie at Arctic Red River
  Kolyma_Kolymskoe_Version_20260210.xlsx               ArcticGRO, Kolyma at Kolymskoe
  Pechora_UstTsilma_Version_20260210.xlsx              ArcticGRO, Pechora at Ust'-Tsil'ma
                                                        *** NOT Pechora_Oksino -- that gauge
                                                        only covers 1980-1993, no overlap
                                                        with 2004-2025 sim period ***
  Yana_Ubileynaya_Version_20260210.xlsx                ArcticGRO, Yana at Ubileynaya
  Indigirka_Indigirskiy_Version_20260210.xlsx          ArcticGRO, Indigirka at Indigirskiy
  Olenek_7.5_km_down_of_Buurs_mouth_Version_20260210.xlsx
                                                        ArcticGRO, Olenek at Buur's mouth
  NorthernDvina_UstPenega_Version_20260210.xlsx        ArcticGRO, Northern Dvina at Ust'-Penega

All files expected to have columns: date (parseable), discharge (m3/s).
"""

import pandas as pd
import numpy as np

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
STEMMLER_PFOA_KG_YR = 14800.0
SIM_YEAR_START      = 2004
SIM_YEAR_END        = 2025

# Rivers with no real discharge data before this year; injection is
# zero-filled (not interpolated) prior to this year.
ZERO_FILL_RIVERS = {
    'Indigirka': 2016,
    'Olenek':    2016,
}

OUTPUT_CSV = 'experiment/riverine_pfoa_flux.csv'

RIVERS = {
    'Lena':          'data/Lena-Discharge.xlsx',
    'Yenisei':       'data/Yenisei-Discharge.xlsx',
    'Ob':            'data/Ob-Discharge.xlsx',
    'Mackenzie':     'data/Mackenzie-Discharge.xlsx',
    'Kolyma':        'data/Kolyma_Kolymskoe_Version_20260210.xlsx',
    # NOTE: Pechora_Oksino_Version_20260210.xlsx is the WRONG gauge --
    # Oksino only covers 1980-1993, no overlap with 2004-2025 sim period.
    # Must use Ust'-Tsil'ma gauge instead (re-download from ArcticGRO).
    'Pechora':       'data/Pechora_UstTsilma_Version_20260210.xlsx',
    'Yana':          'data/Yana_Ubileynaya_Version_20260210.xlsx',
    'Indigirka':     'data/Indigirka_Indigirskiy_Version_20260210.xlsx',
    'Olenek':        'data/Olenek_7.5_km_down_of_Buur.s_mouth_Version_20260210.xlsx',
    'NorthernDvina': 'data/NorthernDvina_UstPenega_Version_20260210.xlsx',
}

# ------------------------------------------------------------------
# Load and preprocess each river
# ------------------------------------------------------------------
dfs = []
for river, path in RIVERS.items():
    df = pd.read_excel(path)
    df['date'] = pd.to_datetime(df['date'])
    df = df[(df['date'].dt.year >= SIM_YEAR_START) &
            (df['date'].dt.year <= SIM_YEAR_END)].copy()
    df['river'] = river
    dfs.append(df)

data = pd.concat(dfs, ignore_index=True)
data['year']  = data['date'].dt.year
data['month'] = data['date'].dt.month

# ------------------------------------------------------------------
# Monthly mean discharge per river/year/month
# ------------------------------------------------------------------
monthly = (data.groupby(['river', 'year', 'month'])['discharge']
               .mean()
               .reset_index()
               .rename(columns={'discharge': 'discharge_m3s'}))

# ------------------------------------------------------------------
# Build complete river x year x month index so gaps show up as NaN
# ------------------------------------------------------------------
complete_index = pd.MultiIndex.from_product(
    [list(RIVERS.keys()),
     range(SIM_YEAR_START, SIM_YEAR_END + 1),
     range(1, 13)],
    names=['river', 'year', 'month']
)
monthly = (monthly
           .set_index(['river', 'year', 'month'])
           .reindex(complete_index)
           .reset_index())

n_gaps_before = monthly['discharge_m3s'].isna().sum()
print(f'Total missing river/month values before gap handling: {n_gaps_before}')

# ------------------------------------------------------------------
# Step 1: zero-fill rivers with a known pre-data-availability gap.
# This must happen BEFORE interpolation so the interpolator never
# bridges across this gap with a fabricated trend.
# ------------------------------------------------------------------
for river, start_year in ZERO_FILL_RIVERS.items():
    mask = (monthly['river'] == river) & (monthly['year'] < start_year)
    n_zeroed = monthly.loc[mask, 'discharge_m3s'].isna().sum()
    monthly.loc[mask, 'discharge_m3s'] = 0.0
    print(f'  {river}: zero-filled {n_zeroed} months before {start_year} '
          f'(no injection prior to {start_year})')

# ------------------------------------------------------------------
# Step 2: fill remaining gaps (e.g. Yana 2000-2017) using calendar-month
# climatology -- for each river, fill a missing (year, month) with the
# mean discharge for that calendar month across all other available
# years for that river.
#
# Linear time-interpolation is NOT used here: for rivers with strong
# seasonality (e.g. Yana, near-zero in winter and >40,000 m3/s at the
# summer freshet peak), a straight line across a multi-year gap applies
# a near-constant intermediate value to every month including winter,
# severely inflating the annual mean. Climatological fill preserves the
# seasonal cycle through the gap.
# ------------------------------------------------------------------
n_gaps_remaining = monthly['discharge_m3s'].isna().sum()
print(f'Remaining gaps to fill via calendar-month climatology: {n_gaps_remaining}')

monthly = monthly.sort_values(['river', 'year', 'month']).reset_index(drop=True)

climatology = (monthly.groupby(['river', 'month'])['discharge_m3s']
                       .transform('mean'))

n_before = monthly['discharge_m3s'].isna().sum()
monthly['discharge_m3s'] = monthly['discharge_m3s'].fillna(climatology)
n_filled = n_before - monthly['discharge_m3s'].isna().sum()
print(f'  Filled {n_filled} gaps using calendar-month climatology')

# Fallback for any (river, month) combination with zero data anywhere
# (should not occur given current inputs, but guard against it).
remaining = monthly['discharge_m3s'].isna().sum()
if remaining > 0:
    print(f'  WARNING: {remaining} gaps remain after climatology fill -- '
          f'falling back to ffill/bfill')
    monthly['discharge_m3s'] = (
        monthly.groupby('river')['discharge_m3s']
               .transform(lambda x: x.ffill().bfill())
    )

assert monthly['discharge_m3s'].isna().sum() == 0, \
    'NaNs remain after gap handling -- check input data coverage'

# ------------------------------------------------------------------
# Discharge fractions from 2004-2025 mean annual discharge.
# Note: for Indigirka/Olenek, the zero-filled 2004-2015 period pulls
# down their mean annual discharge and therefore their fraction.
# This is intentional -- their fraction reflects the fact that they
# contribute flux for only part of the simulation period.
# ------------------------------------------------------------------
mean_annual     = monthly.groupby('river')['discharge_m3s'].mean()
river_fractions = mean_annual / mean_annual.sum()

print('\nMean annual discharge (m3/s, 2004-2025 incl. zero-fill) and PFOA fractions:')
print(f'  {"River":<16} {"Q_mean m3/s":>12} {"Fraction":>10} {"kg/yr":>10}')
print(f'  {"-"*52}')
for river in RIVERS:
    q  = mean_annual[river]
    f  = river_fractions[river]
    kg = STEMMLER_PFOA_KG_YR * f
    print(f'  {river:<16} {q:>12.1f} {f:>10.4f} {kg:>10.1f}')
print(f'  {"-"*52}')
print(f'  {"TOTAL":<16} {mean_annual.sum():>12.1f} {river_fractions.sum():>10.4f} '
      f'{STEMMLER_PFOA_KG_YR:>10.1f}')

# ------------------------------------------------------------------
# Monthly PFOA flux (kg/month) -- chemostatic seasonal shaping
# flux = annual_total_kg_yr * fraction * (Q_month / Q_annual_mean) / 12
#
# For zero-filled years (Indigirka/Olenek pre-2016), Q_annual_mean
# for that year is 0, which would divide by zero. Guard against this:
# zero-filled years get pfoa_flux_kg_month = 0 directly.
# ------------------------------------------------------------------
annual_mean = (monthly.groupby(['river', 'year'])['discharge_m3s']
                      .mean()
                      .reset_index()
                      .rename(columns={'discharge_m3s': 'annual_mean_m3s'}))

monthly = monthly.merge(annual_mean, on=['river', 'year'])

with np.errstate(divide='ignore', invalid='ignore'):
    seasonal_ratio = np.where(
        monthly['annual_mean_m3s'] > 0,
        monthly['discharge_m3s'] / monthly['annual_mean_m3s'],
        0.0
    )

monthly['pfoa_flux_kg_month'] = (
    STEMMLER_PFOA_KG_YR
    * monthly['river'].map(river_fractions)
    * seasonal_ratio
    / 12.0
)

monthly = monthly.drop(columns=['annual_mean_m3s'])

# ------------------------------------------------------------------
# Diagnostics: confirm zero-fill rivers actually inject zero pre-2016
# and nonzero post-2016
# ------------------------------------------------------------------
print('\nZero-fill river injection check:')
for river, start_year in ZERO_FILL_RIVERS.items():
    pre  = monthly[(monthly['river'] == river) &
                    (monthly['year'] <  start_year)]['pfoa_flux_kg_month'].sum()
    post = monthly[(monthly['river'] == river) &
                    (monthly['year'] >= start_year)]['pfoa_flux_kg_month'].sum()
    n_years_post = SIM_YEAR_END - start_year + 1
    print(f'  {river}: pre-{start_year} total={pre:.4f} kg '
          f'(expect 0), post-{start_year} total={post:.2f} kg '
          f'over {n_years_post} yr ({post/n_years_post:.2f} kg/yr)')

# ------------------------------------------------------------------
# Mass conservation check.
# NOTE: total injected mass will be LESS than STEMMLER_PFOA_KG_YR
# during 2004-2015 because Indigirka/Olenek contribute zero, but
# their fraction of the 14,800 kg/yr was still subtracted from the
# other rivers (since fractions are computed over the full 2004-2025
# period and sum to 1.0). This means 2004-2015 years are slightly
# UNDER 14,800 kg/yr total, and 2016-2025 years are at 14,800 kg/yr.
# This is a real, documented effect of the zero-fill approach -- see
# DECISIONS.md sec. 9.
# ------------------------------------------------------------------
total_per_year = (monthly.groupby('year')['pfoa_flux_kg_month']
                         .sum()
                         .reset_index())

print(f'\nMass conservation check (target {STEMMLER_PFOA_KG_YR:.0f} kg/yr '
      f'for 2016-2025; less for 2004-2015 due to Indigirka/Olenek zero-fill):')
pre_2016  = total_per_year[total_per_year['year'] <  2016]
post_2016 = total_per_year[total_per_year['year'] >= 2016]
print(f'  2004-2015: min={pre_2016["pfoa_flux_kg_month"].min():.1f}  '
      f'max={pre_2016["pfoa_flux_kg_month"].max():.1f}  '
      f'mean={pre_2016["pfoa_flux_kg_month"].mean():.1f} kg/yr')
print(f'  2016-2025: min={post_2016["pfoa_flux_kg_month"].min():.1f}  '
      f'max={post_2016["pfoa_flux_kg_month"].max():.1f}  '
      f'mean={post_2016["pfoa_flux_kg_month"].mean():.1f} kg/yr')

# ------------------------------------------------------------------
# Write output
# ------------------------------------------------------------------
monthly.to_csv(OUTPUT_CSV, index=False)
print(f'\nSaved to {OUTPUT_CSV}')
print(f'Rows: {len(monthly)} ({len(RIVERS)} rivers x '
      f'{SIM_YEAR_END - SIM_YEAR_START + 1} years x 12 months = '
      f'{len(RIVERS) * (SIM_YEAR_END - SIM_YEAR_START + 1) * 12} expected)')