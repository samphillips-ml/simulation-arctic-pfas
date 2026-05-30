import pandas as pd
import numpy as np

STEMMLER_PFOA_KG_YR = 14800.0

rivers = {
    'Lena':      'data/Lena-Discharge.xlsx',
    'Mackenzie': 'data/Mackenzie-Discharge.xlsx',
    'Ob':        'data/Ob-Discharge.xlsx',
    'Yenisei':   'data/Yenisei-Discharge.xlsx',
}

dfs = []
for river, path in rivers.items():
    df = pd.read_excel(path)
    df['date'] = pd.to_datetime(df['date'])
    df = df[(df['date'].dt.year >= 2004) & (df['date'].dt.year <= 2025)]
    df['river'] = river
    dfs.append(df)

data = pd.concat(dfs, ignore_index=True)

data['year'] = data['date'].dt.year
data['month'] = data['date'].dt.month
monthly = (data.groupby(['river', 'year', 'month'])['discharge']
               .mean()
               .reset_index()
               .rename(columns={'discharge': 'discharge_m3s'}))

mean_annual = monthly.groupby('river')['discharge_m3s'].mean()
river_fractions = mean_annual / mean_annual.sum()
print("River discharge fractions:")
print(river_fractions.round(3))

annual_mean = (monthly.groupby(['river', 'year'])['discharge_m3s']
                      .mean()
                      .reset_index()
                      .rename(columns={'discharge_m3s': 'annual_mean_m3s'}))
monthly = monthly.merge(annual_mean, on=['river', 'year'])

monthly['pfoa_flux_kg_month'] = (
    STEMMLER_PFOA_KG_YR
    * monthly['river'].map(river_fractions)
    * (monthly['discharge_m3s'] / monthly['annual_mean_m3s'])
    / 12.0
)

monthly = monthly.drop(columns=['annual_mean_m3s'])
monthly.to_csv('experiment/riverine_pfoa_flux.csv', index=False)
print("\nSaved to experiment/riverine_pfoa_flux.csv")
print(monthly.head(12)) # debug, TODO: Remove