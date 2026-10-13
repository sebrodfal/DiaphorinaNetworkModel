"""Experiment 3: replicate stochastic topologies to quantify variance.
"Random", "SmallWorld", "ScaleFree" are redrawn n_seeds times (independent
graph instances) and the point-source outbreak scenario of Experiment 2 is
re-run on each, so results can be reported as mean +/- std per
(topology, Dd, Dt, strategy). Deterministic topologies are included once."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import PARAMS, base_outbreak_ics, results_path, write_csv
from network_experiments import experiment_metrics, make_topology, solve_network_data

n = 8
tt = 100
T = 15
delta_amt = 35
ics = base_outbreak_ics(n, tt)

deterministic_topologies = ["Path", "Cycle", "Star", "Complete"]
stochastic_topologies = ["Random", "SmallWorld", "ScaleFree"]
n_seeds = 10
diffusion_grid = [0, 0.01, 0.03, 0.1]
strategies = ["None", "Uniform", "TargetedAtSource"]
STRATEGY = {"None": (0, [], "None"),
            "Uniform": (delta_amt, "All", "Uniform"),
            "TargetedAtSource": (delta_amt, [0], "TargetedAtSource")}

results = []


def run_one(topo, seed, graph):
    for Dd in diffusion_grid:
        for Dt in diffusion_grid:
            for strategy in strategies:
                delta, release_nodes, label = STRATEGY[strategy]
                data = solve_network_data(graph, Dd, Dt, delta, T, PARAMS, ics, release_nodes)
                row = {"Topology": topo, "Seed": seed, "N": n, "Dd": Dd, "Dt": Dt,
                       "Strategy": label, "Delta": delta, "T": T}
                row.update(experiment_metrics(data))
                results.append(row)


for topo in deterministic_topologies:
    run_one(topo, 1, make_topology(topo, n, 1))

for topo in stochastic_topologies:
    for seed in range(1, n_seeds + 1):
        run_one(topo, seed, make_topology(topo, n, seed))

print("Completed", len(results), "runs.")

raw_path = results_path("experiment3_replicates_raw.csv")
write_csv(raw_path, results)
print("Saved:", raw_path)

# Summary: mean/std per (Topology, Dd, Dt, Strategy), collapsing over Seed
metric_keys = ["PestPeak", "PestMin", "PestAUC", "PestLastCycleMean", "PestSpatialCV",
               "PestSync", "ParasitoidLastCycleMean", "ParasitoidSync", "AlgebraicConnectivity",
               "MeanDegree"]
groups = {}
for r in results:
    groups.setdefault((r["Topology"], r["Dd"], r["Dt"], r["Strategy"]), []).append(r)

summary_rows = []
for key, group in groups.items():
    row = {"Topology": key[0], "Dd": key[1], "Dt": key[2], "Strategy": key[3],
           "NumSeeds": len(group)}
    for m in metric_keys:
        vals = np.array([g[m] for g in group], dtype=float)
        row[m + "_Mean"] = float(vals.mean())
        row[m + "_Std"] = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
    summary_rows.append(row)

summary_path = results_path("experiment3_replicates_summary.csv")
write_csv(summary_path, summary_rows)
print("Saved:", summary_path)

# Figure: PestAUC mean +/- std by topology, Strategy="Uniform", Dd=Dt=0.03
all_topologies = deterministic_topologies + stochastic_topologies
by_topo = {r["Topology"]: r for r in summary_rows
           if r["Dd"] == 0.03 and r["Dt"] == 0.03 and r["Strategy"] == "Uniform"}
mean_vals = [by_topo[t]["PestAUC_Mean"] for t in all_topologies]
std_vals = [by_topo[t]["PestAUC_Std"] for t in all_topologies]
x = np.arange(1, len(all_topologies) + 1)

fig, ax = plt.subplots(figsize=(10, 6))
ax.errorbar(x, mean_vals, yerr=std_vals, fmt="o", color="black", ecolor="black", capsize=5)
ax.set_xticks(x, all_topologies)
ax.set_xlabel("Topology")
ax.set_ylabel("Pest AUC (network-wide, Uniform release, Dd=Dt=0.03)")
ax.set_title(f"Mean \u00b1 std of Pest AUC across {n_seeds} seeds (stochastic topologies)")
ax.grid(True, linestyle="--", alpha=0.6)
png_path = results_path("experiment3_summary_errorbars.png")
fig.savefig(png_path, dpi=120, bbox_inches="tight")
print("Saved:", png_path)
