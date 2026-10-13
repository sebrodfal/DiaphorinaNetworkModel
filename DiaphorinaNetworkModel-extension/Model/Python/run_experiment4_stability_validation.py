"""Experiment 4: validates the Master-Stability-without-periodicity machinery
(network_stability.master_lyapunov_exponent / master_lyapunov_spectrum) and
applies it to the real topologies used in Experiments 2/3.

Part A -- mu=0 must reduce to the plain single-node problem, independent of
Dd/Dt. Checked against an independently written single-node Lyapunov
calculation, and by comparing master_lyapunov_exponent at mu=0 for two Dd/Dt.

Part B -- master_lyapunov_spectrum on the Laplacian spectrum of each of the 7
topologies (n=8, delta=35, T=15, Dd=Dt=0.03), cross-referenced against the
empirical PestSync/PestSpatialCV of Experiment 3's summary CSV."""

import os

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from common import PARAMS, results_path, write_csv
from core import (ATOL, METHOD, RTOL, local_jacobian_function, release_times, single_node_rhs,
                  solve_impulsive)
from network_experiments import laplacian, make_topology
from network_stability import master_lyapunov_exponent, master_lyapunov_spectrum

params = PARAMS
delta_base, T_base, Dd_base, Dt_base, n = 35, 15, 0.03, 0.03, 8

# ============================================================================
# PART A -- ground truth check: mu=0 must reduce to the single-node problem
# ============================================================================
print("### PART A: mu=0 reduction check ###")


def raw_single_node_lyapunov(params, delta, T, t_total, burn_frac, seg_len):
    """Independent single-node Lyapunov exponent: its own reference solve and
    its own segment loop (only the local field/Jacobian are shared)."""
    Jf = local_jacobian_function(params)
    jumps = [] if delta == 0 else release_times(T, t_total)

    def jump(y):
        y = y.copy()
        y[2] += delta
        return y
    w = solve_impulsive(single_node_rhs(params), [50, 10, 15, 30], (0.0, float(t_total)), jumps, jump)
    burn = burn_frac * t_total
    n_seg = int(np.floor((t_total - burn) / seg_len))
    t_c, v, log_sum = burn, np.array([1.0, 0.0, 0.0, 0.0]), 0.0
    for _ in range(n_seg):
        r = solve_ivp(lambda t, u: Jf(*w(t)) @ u, (t_c, t_c + seg_len), v,
                      method=METHOD, rtol=RTOL, atol=ATOL)
        v_end = r.y[:, -1]
        nrm = np.linalg.norm(v_end)
        v, log_sum, t_c = v_end / nrm, log_sum + np.log(nrm), t_c + seg_len
    return log_sum / (n_seg * seg_len)


raw_lambda0 = raw_single_node_lyapunov(params, delta_base, T_base, 1200, 0.3, T_base)
via_master_dd0 = master_lyapunov_exponent(params, 0, 0, 0, delta_base, T_base)["LyapunovExponent"]
via_master_ddreal = master_lyapunov_exponent(params, Dd_base, Dt_base, 0, delta_base, T_base)["LyapunovExponent"]

print("Raw single-node lambda (independent calc):     ", raw_lambda0)
print("master_lyapunov_exponent, mu=0, Dd=Dt=0:        ", via_master_dd0)
print(f"master_lyapunov_exponent, mu=0, Dd=Dt={Dd_base}:  ", via_master_ddreal)
diffs = np.abs(np.array([raw_lambda0, via_master_dd0, via_master_ddreal]) - raw_lambda0)
print("Max pairwise diff (should be ~0):", diffs.max())
part_a_pass = diffs.max() < 1e-4
print("PASS: mu=0 reduces correctly to the single-node problem, independent of Dd/Dt." if part_a_pass
      else "FAIL: mismatch -- check reference_trajectory/perturbation_segment.")

# ============================================================================
# PART B -- master_lyapunov_spectrum on the real topologies from Experiments 2/3
# ============================================================================
print("\n### PART B: synchronized-state stability prediction per topology ###")

topologies = ["Path", "Cycle", "Star", "Complete", "Random", "SmallWorld", "ScaleFree"]


def topology_result(topo):
    graph = make_topology(topo, n)
    mu_all = np.sort(np.linalg.eigvalsh(laplacian(graph)))
    mu_transverse = mu_all[1:]          # drop mu_1=0, the synchronized mode itself
    spectrum = master_lyapunov_spectrum(params, Dd_base, Dt_base, mu_transverse, delta_base, T_base)
    max_lambda = max(s["LyapunovExponent"] for s in spectrum)
    return {"Topology": topo, "AlgebraicConnectivity": float(mu_all[1]),
            "MaxTransverseLyapunov": float(max_lambda),
            "PredictedSyncStable": bool(max_lambda < 0),
            "AllReliable": all(s["Reliable"] for s in spectrum)}


results = [topology_result(t) for t in topologies]
csv_path = results_path("experiment4_stability_prediction.csv")
write_csv(csv_path, results)
print("Saved:", csv_path)
for r in results:
    print(f"  {r['Topology']}: max transverse lambda={r['MaxTransverseLyapunov']}"
          f"  predicted stable sync={r['PredictedSyncStable']}  all segments reliable={r['AllReliable']}")

# Cross-reference against Experiment 3 (same Dd=Dt=0.03, Uniform strategy)
summary_path = results_path("experiment3_replicates_summary.csv")
if os.path.exists(summary_path):
    summary = pd.read_csv(summary_path)
    sl = summary[(abs(summary.Dd - Dd_base) < 1e-9) & (abs(summary.Dt - Dt_base) < 1e-9)
                 & (summary.Strategy == "Uniform")].set_index("Topology")
    comparison = []
    for r in results:
        emp = sl.loc[r["Topology"]]
        comparison.append({"Topology": r["Topology"], "PredictedSyncStable": r["PredictedSyncStable"],
                           "MaxTransverseLyapunov": r["MaxTransverseLyapunov"],
                           "EmpiricalPestSync_Mean": float(emp["PestSync_Mean"]),
                           "EmpiricalPestSpatialCV_Mean": float(emp["PestSpatialCV_Mean"])})
    print("\n### Theory vs. simulation (Experiment 3, Dd=Dt=0.03, Uniform) ###")
    for c in comparison:
        print(f"  {c['Topology']}: predicted stable={c['PredictedSyncStable']}"
              f"  (lambda={c['MaxTransverseLyapunov']})   empirical PestSync={c['EmpiricalPestSync_Mean']}"
              f"  PestSpatialCV={c['EmpiricalPestSpatialCV_Mean']}")
    comp_path = results_path("experiment4_theory_vs_simulation.csv")
    write_csv(comp_path, comparison)
    print("Saved:", comp_path)
else:
    print("\nSkipped theory-vs-simulation cross-reference: could not find", summary_path,
          "(run run_experiment3_replicates.py first)")

print("\nDONE.")
