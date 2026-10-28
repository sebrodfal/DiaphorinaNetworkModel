"""Experiment 5: bifurcation map of the SINGLE-NODE impulsive model over
(delta, T). Classifies each combination as "Regular" (settles to a period-T
or low-order period-kT cycle) or "Irregular" (wide dispersion) using the
dispersion of the last 30 post-jump states out of a 120-cycle run, and saves
a classification CSV + heatmap PNG."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import PARAMS, results_path, write_csv
from core import release_times, single_node_rhs, solve_impulsive
from functions import TEMPERATURE_MAP

params = PARAMS
n_cycles = 120
n_tail = 30


def classify_regime(delta_amt, T):
    tt = n_cycles * T

    def jump(y):
        y = y.copy()
        y[2] += delta_amt
        return y
    sol = solve_impulsive(single_node_rhs(params), [50, 10, 15, 30], (0.0, float(tt)),
                          release_times(T, tt), jump)
    tail = sol(np.array([k * T + 1e-6 for k in range(n_cycles - n_tail, n_cycles)])).T
    centre = tail.mean(axis=0)
    dispersion = float(np.mean(np.linalg.norm(tail - centre, axis=1)))
    return {"Delta": delta_amt, "T": T, "Dispersion": dispersion,
            "Regime": "Regular" if dispersion < 1.0 else "Irregular", "s0": tail[-1]}


delta_grid = [5, 10, 15, 20, 25, 30, 35, 40]
T_grid = [5, 8, 11, 15, 20, 25, 30, 40]
results = [classify_regime(d, tp) for d in delta_grid for tp in T_grid]
print("Completed", len(results), "classifications.")

csv_path = results_path("experiment5_bifurcation_map.csv")
write_csv(csv_path, [{k: v for k, v in r.items() if k != "s0"} for r in results])
print("Saved:", csv_path)

# Heatmap of log10(dispersion), T on x-axis, delta on y-axis (delta grows upwards)
disp = np.array([[next(r for r in results if r["Delta"] == d and r["T"] == tp)["Dispersion"]
                  for tp in T_grid] for d in delta_grid])
log_disp = np.log10(disp + 1e-3)
fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(log_disp[::-1], cmap=TEMPERATURE_MAP, aspect="auto")
ax.set_xticks(range(len(T_grid)), T_grid)
ax.set_yticks(range(len(delta_grid)), delta_grid[::-1])
ax.set_xlabel("T (release period)")
ax.set_ylabel("\u03b4 (release amount)")
ax.set_title("log10(dispersion of last 30 post-jump states) \u2014 low = periodic, high = irregular/chaotic",
             fontsize=10)
fig.colorbar(im, ax=ax)
png_path = results_path("experiment5_bifurcation_heatmap.png")
fig.savefig(png_path, dpi=120, bbox_inches="tight")
print("Saved:", png_path)

print("Regular (delta,T) points:", [(r["Delta"], r["T"]) for r in results if r["Regime"] == "Regular"])
