"""Experiment 7: one-at-a-time local sensitivity analysis of the 18 local
parameters on the point-source outbreak scenario (tt=100). Each parameter is
perturbed +/-10% (topology=Complete, Dd=Dt=0.03, delta=35, T=15, Uniform
release) and the % change in PestAUC and PestLastCycleMean is recorded."""

import copy

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import PARAMS, base_outbreak_ics, results_path, write_csv
from network_experiments import experiment_metrics, make_topology, solve_network_data

base_params = PARAMS

# (block, index, name) with 1-based block/index as in Mathematica
param_names = [
    (1, 1, "r1"), (1, 2, "b12"), (1, 3, "b13"), (1, 4, "a1"), (1, 5, "c1"),
    (2, 1, "r2"), (2, 2, "b21"), (2, 3, "b24"), (2, 4, "a2"), (2, 5, "c2"),
    (3, 1, "r3"), (3, 2, "b31"), (3, 3, "a3"), (3, 4, "c3"),
    (4, 1, "r4"), (4, 2, "b42"), (4, 3, "a4"), (4, 4, "c4"),
]

n = 8
tt = 100
Dd_val, Dt_val = 0.03, 0.03
delta_amt, T_val = 35, 15
ics = base_outbreak_ics(n, tt)
graph = make_topology("Complete", n)


def run_with_params(p):
    return experiment_metrics(solve_network_data(graph, Dd_val, Dt_val, delta_amt, T_val, p, ics, "All"))


base_metrics = run_with_params(base_params)
print("Baseline PestAUC =", base_metrics["PestAUC"], "   PestLastCycleMean =", base_metrics["PestLastCycleMean"])


def perturb(block, idx, factor):
    p = copy.deepcopy(base_params)
    p[block - 1][idx - 1] = base_params[block - 1][idx - 1] * factor
    return p


results = []
for block, idx, name in param_names:
    for factor in (0.9, 1.1):
        m = run_with_params(perturb(block, idx, factor))
        results.append({
            "Parameter": name, "Factor": factor,
            "PerturbedValue": base_params[block - 1][idx - 1] * factor,
            "PestAUC": m["PestAUC"],
            "PctChangeAUC": 100.0 * (m["PestAUC"] - base_metrics["PestAUC"]) / base_metrics["PestAUC"],
            "PestLastCycleMean": m["PestLastCycleMean"],
            "PctChangeLastCycleMean": 100.0 * (m["PestLastCycleMean"] - base_metrics["PestLastCycleMean"])
            / base_metrics["PestLastCycleMean"],
        })

print("Completed", len(results), "perturbation runs.")
csv_path = results_path("experiment7_sensitivity.csv")
write_csv(csv_path, results)
print("Saved:", csv_path)

# Sensitivity score per parameter: max |%change in AUC| across the two perturbations
sensitivity = []
for _, _, name in param_names:
    vals = [abs(r["PctChangeAUC"]) for r in results if r["Parameter"] == name]
    sensitivity.append({"Parameter": name, "MaxAbsPctChangeAUC": max(vals)})
sensitivity_sorted = sorted(sensitivity, key=lambda s: s["MaxAbsPctChangeAUC"], reverse=True)

rank_path = results_path("experiment7_sensitivity_ranking.csv")
write_csv(rank_path, sensitivity_sorted)
print("Saved:", rank_path)

print("Top 5 most sensitive parameters: ")
for s in sensitivity_sorted[:5]:
    print(f"  {s['Parameter']}: {s['MaxAbsPctChangeAUC']}% max change in PestAUC")

# Tornado chart of top 10
top10 = sensitivity_sorted[:10]
names = [t["Parameter"] for t in top10]
low = [next(r for r in results if r["Parameter"] == nm and r["Factor"] == 0.9)["PctChangeAUC"] for nm in names]
high = [next(r for r in results if r["Parameter"] == nm and r["Factor"] == 1.1)["PctChangeAUC"] for nm in names]
x = np.arange(len(names))
fig, ax = plt.subplots(figsize=(11, 6))
w = 0.4
ax.bar(x - w / 2, low, w, label="param \u00d7 0.9")
ax.bar(x + w / 2, high, w, label="param \u00d7 1.1")
ax.set_xticks(x, names)
ax.set_title("Sensitivity of PestAUC to \u00b110% parameter perturbations (top 10, sorted)")
ax.set_xlabel("Parameter")
ax.set_ylabel("% change in PestAUC")
ax.grid(True, linestyle="--", alpha=0.6)
ax.legend()
png_path = results_path("experiment7_tornado.png")
fig.savefig(png_path, dpi=120, bbox_inches="tight")
print("Saved:", png_path)
