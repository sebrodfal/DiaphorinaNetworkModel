"""Experiment 6: dose-frequency sweep (delta x T) crossed with topology, over the
short horizon used throughout the preprint (tt=100).

Cost-effectiveness: TotalReleased = delta * (number of release events in
[0,tt]) * (number of release nodes). CostEfficiency = AUC avoided per unit of
parasitoid released, relative to the no-control baseline for the same topology."""

import numpy as np

from common import PARAMS, base_outbreak_ics, results_path, write_csv
from network_experiments import experiment_metrics, make_topology, solve_network_data

n = 8
tt = 100
Dd_val = 0.03
Dt_val = 0.03
ics = base_outbreak_ics(n, tt)

topologies = ["Path", "Star", "Complete", "Random"]
delta_grid = [10, 20, 30, 40, 50]
T_grid = [5, 10, 15, 20, 30]
strategies = ["Uniform", "TargetedAtSource"]

results = []
baseline_auc = {}

for topo in topologies:
    graph = make_topology(topo, n)
    # No-control baseline for this topology, to normalize cost-effectiveness
    base_data = solve_network_data(graph, Dd_val, Dt_val, 0, 15, PARAMS, ics, [])
    baseline_auc[topo] = experiment_metrics(base_data)["PestAUC"]
    for delta_amt in delta_grid:
        for T_val in T_grid:
            for strategy in strategies:
                release_nodes = "All" if strategy == "Uniform" else [0]
                n_release_nodes = n if strategy == "Uniform" else 1
                data = solve_network_data(graph, Dd_val, Dt_val, delta_amt, T_val, PARAMS, ics, release_nodes)
                metrics = experiment_metrics(data)
                n_events = int(np.floor(tt / T_val))
                total_released = delta_amt * n_events * n_release_nodes
                auc_avoided = baseline_auc[topo] - metrics["PestAUC"]
                row = {"Topology": topo, "Delta": delta_amt, "T": T_val, "Strategy": strategy,
                       "TotalReleased": total_released, "BaselineAUC": baseline_auc[topo],
                       "AUCAvoided": auc_avoided,
                       "CostEfficiency": auc_avoided / total_released if total_released > 0 else 0}
                row.update(metrics)
                results.append(row)

print("Completed", len(results), "runs.")
csv_path = results_path("experiment6_dose_frequency.csv")
write_csv(csv_path, results)
print("Saved:", csv_path)

# Best (delta,T) per topology x strategy by CostEfficiency
groups = {}
for r in results:
    groups.setdefault((r["Topology"], r["Strategy"]), []).append(r)
best = [max(g, key=lambda r: r["CostEfficiency"]) for g in groups.values()]
best_path = results_path("experiment6_best_dose_frequency.csv")
write_csv(best_path, best)
print("Saved:", best_path)
for b in best:
    print(f"{b['Topology']} / {b['Strategy']}: best delta={b['Delta']} T={b['T']}"
          f"  CostEfficiency={b['CostEfficiency']}  PestAUC={b['PestAUC']}")
