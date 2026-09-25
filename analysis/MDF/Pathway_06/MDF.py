import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from equilibrator_api import ComponentContribution
from equilibrator_pathway import ThermodynamicModel

plt.rcParams.update({
    'font.size': 18,
    'axes.labelsize': 20,
    'axes.titlesize': 26,
    'axes.labelweight': 'bold',
    'axes.titleweight': 'bold',
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'axes.edgecolor': 'black',
    'axes.linewidth': 2,
    'figure.figsize': [12, 7],
    'legend.fontsize': 16,
    'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans', 'sans-serif'],
    'axes.grid': False,
    'savefig.dpi': 300
})
sns.set_palette("Set2")

# Load pathway
comp_contrib = ComponentContribution()
pp = ThermodynamicModel.from_sbtab("pathways.tsv", comp_contrib=comp_contrib)
mdf_result = pp.mdf_analysis()
df = mdf_result.reaction_df.copy()

# Convert units for plotting and saving
df['Driving Force (kJ/mol)'] = df['optimized_dg_prime'].apply(lambda x: x.to('kilojoule / mole').m)
df['Standard dg_prime (kJ/mol)'] = df['standard_dg_prime'].apply(lambda x: x.to('kilojoule / mole').m)

# ----------------- Excel Export -----------------
# Save reaction IDs and driving forces to Excel, update column names
df[['reaction_id', 'Driving Force (kJ/mol)', 'Standard dg_prime (kJ/mol)']].rename(
    columns={
        'reaction_id': 'Reaction ID',
        'Driving Force (kJ/mol)': 'MDF Driving Force (kJ/mol)',
        'Standard dg_prime (kJ/mol)': 'Standard Driving Force (kJ/mol)'
    }
).to_excel("reactions.xlsx", index=False)
print("Saved reaction driving forces to reactions.xlsx")
# ------------------------------------------------

# Cumulative sums
df['Cumulative MDF-optimized'] = np.cumsum(df['Driving Force (kJ/mol)'])
df['Cumulative Phys'] = np.cumsum(df['Standard dg_prime (kJ/mol)'])
reaction_ids = list(df['reaction_id'])
bottleneck_idx = df['Driving Force (kJ/mol)'].idxmin()

# ----- Plot -----
fig, ax = plt.subplots(figsize=(12,7))
# Dotted green for standard ΔG' (proxy for physiological)
ax.plot(
    reaction_ids,
    df['Cumulative Phys'],
    linestyle='--', color='#66c2a5', linewidth=3,
    label='Physiological concentrations (standard $\Delta G\'^\circ$)'
)
# Solid grey for MDF-optimized
ax.plot(
    reaction_ids,
    df['Cumulative MDF-optimized'],
    linestyle='-', color='dimgray', linewidth=4,
    label='MDF-optimized concentrations'
)
# Bottleneck in red
if bottleneck_idx > 0:
    ax.plot(
        [reaction_ids[bottleneck_idx-1], reaction_ids[bottleneck_idx]],
        [df['Cumulative MDF-optimized'][bottleneck_idx-1], df['Cumulative MDF-optimized'][bottleneck_idx]],
        linestyle='-', color='red', linewidth=6,
        label='Bottleneck reactions'
    )
ax.set_xlabel("Reaction Step", labelpad=15, fontweight='bold')
ax.set_ylabel(r"Cumulative $\Delta G^\circ$ (kJ/mol)", labelpad=15, fontweight='bold')
plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
ax.legend(frameon=True, fontsize=16, loc='best')
sns.despine(ax=ax)
plt.tight_layout()
plt.show()