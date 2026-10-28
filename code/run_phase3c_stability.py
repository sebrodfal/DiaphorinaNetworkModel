"""Fase 3c: stability of the synchronous state with transit mortality.

With m > 0 a synchronous state x_i(t) = s(t) exists only when every node has the
same degree k (here: cycle, k = 2, and complete graph, k = 7). On a k-regular
graph L = kI - A, so for a Laplacian eigenvector v with eigenvalue mu,
    (L + m A) v = (m k + (1 - m) mu) v.
The synchronous mode (mu = 0) therefore follows the single-node model with
    r1' = r1 - D_d m_d k,   r3' = r3 - D_t m_t k,
and each transverse mode is the master-stability problem of that modified node
with coupling (1 - m) D instead of D. This script
  1. checks the reduction numerically (complete graph and cycle, uniform release,
     synchronous initial condition, m_d = m_t > 0), and
  2. computes the largest transverse Lyapunov exponent for m_t in
     {0, 0.05, 0.1, 0.2, 0.4} (m_d = 0, D_d = D_t = 0.03, delta = 35, T = 15),
     which for m_t = 0 must coincide with experiment4_stability_prediction.csv.
Output: Results/plan/phase3/phase3c_stability.csv, phase3c_reduction_check.json
"""

import copy
import json
import os

import numpy as np
import pandas as pd

from common import PARAMS, results_path, write_csv
from core import release_times, single_node_rhs, solve_impulsive
from network_experiments import laplacian, make_topology, solve_network_data
from network_stability import master_lyapunov_spectrum

OUT = results_path(os.path.join("plan", "phase3"))
os.makedirs(OUT, exist_ok=True)


def effective_params(Dd, Dt, md, mt, k, params=PARAMS):
    p = copy.deepcopy(params)
    p[0][0] = params[0][0] - Dd * md * k
    p[2][0] = params[2][0] - Dt * mt * k
    return p


def reduction_check():
    out = []
    for topo in ("Cycle", "Complete"):
        g = make_topology(topo, 8)
        k = 2 if topo == "Cycle" else 7
        for D, m in [(0.03, 0.2), (0.3, 0.4)]:
            ics = [[50] * 8, [10] * 8, [15] * 8, [30] * 8, 100]
            net = solve_network_data(g, D, D, 35, 15, PARAMS, ics, "All", md=m, mt=m)
            pe = effective_params(D, D, m, m, k)

            def jump(y):
                y = y.copy(); y[2] += 35; return y
            ref = solve_impulsive(single_node_rhs(pe), [50, 10, 15, 30], (0.0, 100.0),
                                  release_times(15, 100), jump)
            ts = np.linspace(0, 100, 401)
            Y, R = net["Solution"](ts), ref(ts)
            err = max(np.max(np.abs(Y[v * 8:(v + 1) * 8] - R[v])) / np.max(np.abs(R[v])) for v in range(4))
            out.append({"Topology": topo, "D": D, "m": m, "max_rel_err": float(err)})
    return out


if __name__ == "__main__":
    chk = reduction_check()
    print("reduction check:", chk)
    json.dump({"cases": chk, "pass": all(c["max_rel_err"] < 1e-6 for c in chk)},
              open(os.path.join(OUT, "phase3c_reduction_check.json"), "w"), indent=1)
    rows = []
    Dd = Dt = 0.03
    for topo in ("Cycle", "Complete"):
        g = make_topology(topo, 8)
        k = 2 if topo == "Cycle" else 7
        mu = np.sort(np.linalg.eigvalsh(laplacian(g)))[1:]
        mu = np.unique(np.round(mu, 10))
        for mt in (0.0, 0.05, 0.1, 0.2, 0.4):
            spec = master_lyapunov_spectrum(effective_params(Dd, Dt, 0.0, mt, k), Dd, (1 - mt) * Dt,
                                            mu, 35, 15)
            lam = max(s["LyapunovExponent"] for s in spec)
            rows.append({"Topology": topo, "k": k, "mt": mt, "r3_eff": PARAMS[2][0] - Dt * mt * k,
                         "MaxTransverseLyapunov": lam, "PredictedSyncStable": bool(lam < 0),
                         "AllReliable": all(s["Reliable"] for s in spec)})
            print(rows[-1], flush=True)
    write_csv(os.path.join(OUT, "phase3c_stability.csv"), rows)
    ref = pd.read_csv(results_path("experiment4_stability_prediction.csv"))
    for r in rows:
        if r["mt"] == 0:
            e = ref[ref.Topology == r["Topology"]].MaxTransverseLyapunov.iloc[0]
            print(f"{r['Topology']}: m_t=0 gives {r['MaxTransverseLyapunov']:.6g}, experiment 4 gives {e:.6g}")
