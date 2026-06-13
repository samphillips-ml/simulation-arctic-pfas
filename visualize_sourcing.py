"""
visualize_riverine.py

Simple time series of PFOA flux (kg/month) per river, 2004-2025.

Run from repo root:
  python visualize_riverine.py

Output: output/riverine_flux.png
"""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

RIVERINE_CSV = "data/riverine_pfoa_flux.csv"
OUT_PATH     = "output/riverine_flux.png"
os.makedirs("output", exist_ok=True)

riv = pd.read_csv(RIVERINE_CSV)

# build a single datetime-like x axis
riv['t'] = (riv['year'] - 2004) * 12 + (riv['month'] - 1)

colors = {
    'Lena':      'steelblue',
    'Yenisei':   'darkorange',
    'Ob':        'green',
    'Mackenzie': 'red',
}

fig, ax = plt.subplots(figsize=(16, 5))

for river, group in riv.groupby('river'):
    group = group.sort_values('t')
    ax.plot(group['t'], group['pfoa_flux_kg_month'],
            color=colors.get(river, 'gray'),
            lw=1.2, label=river)

# year tick marks
year_ticks = [i * 12 for i in range(22)]
ax.set_xticks(year_ticks)
ax.set_xticklabels([str(2004 + i) for i in range(22)], rotation=45, fontsize=8)

ax.set_xlabel("Year")
ax.set_ylabel("PFOA flux (kg/month)")
ax.set_title("Riverine PFOA flux by river, 2004-2025", fontweight="bold")
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
print(f"Saved: {OUT_PATH}")