"""Experiment 2: point-source outbreak containment across network topologies.
Node 1 starts as an outbreak focus (high Diaphorina); all other nodes start
pest-free-ish. Sweeps topology, pest/parasitoid diffusion rates (Dd,Dt), and
release strategy (no control / uniform release at every node / targeted
release only at the outbreak node).
Output: a long-format CSV in Results/, one row per run, plus a summary PNG."""

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

topologies = ["Path", "Cycle", "Star", "Complete", "Random", "SmallWorld", "ScaleFree"]
diffusion_grid = [0, 0.01, 0.03, 0.1]
strategies = ["None", "Uniform", "TargetedAtSource"]
STRATEGY = {"None": (0, [], "None"),
            "Uniform": (delta_amt, "All", "Uniform"),
            "TargetedAtSource": (delta_amt, [0], "TargetedAtSource")}

results = []
for topo in topologies:
    graph = make_topology(topo, n)
    for Dd in diffusion_grid:
        for Dt in diffusion_grid:
            for strategy in strategies:
                delta, release_nodes, label = STRATEGY[strategy]
                data = solve_network_data(graph, Dd, Dt, delta, T, PARAMS, ics, release_nodes)
                row = {"Topology": topo, "N": n, "Dd": Dd, "Dt": Dt,
                       "Strategy": label, "Delta": delta, "T": T}
                row.update(experiment_metrics(data))
                results.append(row)

print("Completed", len(results), "runs.")
csv_path = results_path("experiment2_outbreak_containment.csv")
write_csv(csv_path, results)
print("Saved:", csv_path)

# Quick summary figure: mean PestAUC by strategy at Dd=Dt=0.03, averaged over topologies
subset = [r for r in results if r["Dd"] == 0.03 and r["Dt"] == 0.03]
means = [np.mean([r["PestAUC"] for r in subset if r["Strategy"] == s]) for s in strategies]
fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(strategies, means)
ax.set_title("Mean network Pest AUC by release strategy (Dd=Dt=0.03, avg over topologies)")
ax.set_xlabel("Release strategy")
ax.set_ylabel("Total Diaphorina AUC (network-wide)")
ax.grid(True, linestyle="--", alpha=0.6)
png_path = results_path("experiment2_summary_strategy.png")
fig.savefig(png_path, dpi=120, bbox_inches="tight")
print("Saved:", png_path)
