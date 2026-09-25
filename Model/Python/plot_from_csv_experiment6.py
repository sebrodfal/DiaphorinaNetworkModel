"""Summary figure for Experiment 6, read straight from the saved CSVs."""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import results_path

data = pd.read_csv(results_path("experiment6_dose_frequency.csv"))
delta_grid = [10, 20, 30, 40, 50]
T_grid = [5, 10, 15, 20, 30]
colors = plt.get_cmap("tab10").colors

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 13), gridspec_kw={"hspace": 0.35})

# Panel 1: CostEfficiency vs delta, one line per T, topology=Complete, Uniform
s1 = data[(data.Topology == "Complete") & (data.Strategy == "Uniform")]
for k, T in enumerate(T_grid):
    ys = [float(s1[(s1["Delta"] == d) & (s1["T"] == T)].iloc[0]["CostEfficiency"]) for d in delta_grid]
    ax1.plot(delta_grid, ys, "-o", color=colors[k], label=f"T={T}")
ax1.set_xlabel("\u03b4 (release amount)")
ax1.set_ylabel("Cost efficiency (AUC avoided / parasitoid released)")
ax1.set_title("Complete topology, Uniform release: cost efficiency vs dose, by period T")
ax1.grid(True, linestyle="--", alpha=0.6)
ax1.legend(title="Period", loc="upper right")

# Panel 2: PestAUC by topology, best (Uniform) vs best (TargetedAtSource)
best = pd.read_csv(results_path("experiment6_best_dose_frequency.csv"))
topologies = ["Path", "Star", "Complete", "Random"]
uni = [float(best[(best.Topology == t) & (best.Strategy == "Uniform")].iloc[0]["PestAUC"]) for t in topologies]
tar = [float(best[(best.Topology == t) & (best.Strategy == "TargetedAtSource")].iloc[0]["PestAUC"])
       for t in topologies]
x = np.arange(len(topologies))
w = 0.38
ax2.bar(x - w / 2, uni, w, color=colors[0], label="Uniform")
ax2.bar(x + w / 2, tar, w, color=colors[1], label="TargetedAtSource")
ax2.set_xticks(x, topologies)
ax2.set_xlabel("Topology")
ax2.set_ylabel("Pest AUC")
ax2.set_title("Best-cost-efficiency PestAUC by topology and strategy")
ax2.grid(True, linestyle="--", alpha=0.6)
ax2.legend()

png_path = results_path("experiment6_summary.png")
fig.savefig(png_path, dpi=120, bbox_inches="tight")
print("Saved:", png_path)
