"""Re-draws the Experiment 3 summary figure straight from the saved CSV
(decoupled from the heavy simulation run)."""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import results_path

summary = pd.read_csv(results_path("experiment3_replicates_summary.csv"))
all_topologies = ["Path", "Cycle", "Star", "Complete", "Random", "SmallWorld", "ScaleFree"]

sl = summary[(abs(summary.Dd - 0.03) < 1e-9) & (abs(summary.Dt - 0.03) < 1e-9)
             & (summary.Strategy == "Uniform")].set_index("Topology")
mean_vals = [float(sl.loc[t, "PestAUC_Mean"]) for t in all_topologies]
std_vals = [float(sl.loc[t, "PestAUC_Std"]) for t in all_topologies]
x = np.arange(1, len(all_topologies) + 1)

fig, ax = plt.subplots(figsize=(10, 6))
ax.errorbar(x, mean_vals, yerr=std_vals, fmt="none", ecolor=(0.85, 0.5, 0.1), elinewidth=2.5, capsize=6)
ax.plot(x, mean_vals, "o", color="black", ms=8)
ax.set_xticks(x, all_topologies)
ax.set_xlabel("Topology")
ax.set_ylabel("Pest AUC (network-wide, Uniform release, Dd=Dt=0.03)")
ax.set_title("Mean \u00b1 std of Pest AUC (10 seeds for stochastic topologies)")
ax.grid(True, linestyle="--", alpha=0.6)

png_path = results_path("experiment3_summary_errorbars.png")
fig.savefig(png_path, dpi=120, bbox_inches="tight")
print("Saved:", png_path)
print("Means:", mean_vals)
print("Stds:", std_vals)
