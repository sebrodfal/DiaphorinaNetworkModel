"""Fase 4 (optional): transit mortality that depends on the edge.

Only meaningful if Fase 3 passes its gate (Results/plan/phase3/gate_phase4.json).

Model change: m_t A is replaced by the matrix M_t with entries m_ij A_ij, where
m_ij = 1 - exp(-d_ij / ell) and d_ij is the Euclidean distance between nodes
placed uniformly at random in the unit square (random geometric graph, 8 nodes,
connection radius 0.5). Disconnected graphs are discarded and the next seed is
used, until 10 connected graphs are collected. ell in {0.5, 1, 2}.

Design: D_t in 15 log-spaced values in [0.03, 3] (the Fase-3 crest lies in
[0.1, 1]), D_d = 0.03, CI-B, r3 = -0.3, m_d = 0, delta = 35, T = 15, t_max = 100.
10 graphs x 3 ell x 15 D_t = 450 combinations = 900 runs (None and Targeted).

For each (graph, ell) the observed D_t* (argmax of E_T, parabolic refinement in
log D_t) is compared with the Fase-3 rule D_t* = c |r3| / (m_bar k_s), with c the
median of D_t* m_t k_s / |r3| over the real Fase-3 maxima (CI-B, r3 = -0.3), and two
choices of m_bar:
  (1) the mean of m_ij over the edges leaving the focus (the plan's first choice);
  (2) the mean of m_ij over all edges weighted by the number of parasitoids that
      cross each edge in the targeted run at D_t*, sum_ij A_ij D_t int x3_i dt.
      (Weighting only the focus edges by flux would give back (1), because the
      outflow of the focus is split equally among its edges.)
A third prediction (not in the plan) uses the variant that Fase 3 favoured,
with the mean degree <k> of the network instead of k_s and the mean m_ij over
all edges: D_t* = c_k |r3| / (m_all <k>), c_k = median of D_t* m_t <k> / |r3|.
Output: Results/plan/phase4/phase4_grid.csv, phase4_summary.csv, fig4_1_observed_vs_predicted.png
"""

import json
import os
from multiprocessing import Pool

import numpy as np
import networkx as nx
import pandas as pd

import core
from common import outbreak_ics, results_path
from core import n_releases
from network_experiments import adjacency, solve_network_data
import plan_common as pc

OUT = results_path(os.path.join("plan", "phase4"))
ELLS = [0.5, 1.0, 2.0]
DT_GRID = list(np.logspace(np.log10(0.03), np.log10(3), 15))
DD, R3, RADIUS, N = 0.03, -0.3, 0.5, 8


def geometric_graphs(n_graphs=10):
    out, seed = [], 0
    while len(out) < n_graphs:
        seed += 1
        g = nx.random_geometric_graph(N, RADIUS, seed=seed)
        if nx.is_connected(g):
            pos = np.array([g.nodes[i]["pos"] for i in range(N)])
            out.append((seed, g, pos))
    return out


def mortality_matrix(pos, ell):
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=2)
    return 1 - np.exp(-d / ell)


def run_one(args):
    seed, ell, Dt = args
    core.RTOL, core.ATOL = pc.SWEEP_RTOL, pc.SWEEP_ATOL
    g, pos = next((g, p) for s, g, p in GRAPHS if s == seed)
    A = adjacency(g)
    M = mortality_matrix(pos, ell) * A
    ics = outbreak_ics(N, 100, "B")
    params = pc.params_with_r3(R3)
    res = {}
    for strat, nodes, delta in (("None", [], 0), ("Targeted", [0], pc.DELTA)):
        data = solve_network_data(g, DD, Dt, delta, pc.T_REL, params, ics, nodes, md=0.0, mt=M)
        v = data["Solution"].integrate(0.0, 100.0)
        res[strat] = (float(v[:N].sum()), v[2 * N:3 * N])
    k_s = int(A[0].sum())
    Q = pc.DELTA * n_releases(pc.T_REL, 100)
    m_focus = float(M[0][A[0] > 0].mean())
    x3int = res["Targeted"][1]
    flux = A * x3int[:, None]            # parasitoids leaving i towards j (up to D_t)
    m_flux = float((M * flux).sum() / flux.sum())
    return {"Seed": seed, "ell": ell, "Dt": Dt, "Dd": DD, "r3": R3, "k_s": k_s,
            "AUC_None": res["None"][0], "AUC_Targeted": res["Targeted"][0],
            "E_T": (res["None"][0] - res["Targeted"][0]) / Q,
            "m_focus_mean": m_focus, "m_flux_weighted": m_flux,
            "m_all_edges_mean": float(M[A > 0].mean())}


GRAPHS = geometric_graphs()

if __name__ == "__main__":
    from analyze_phase2 import optimum
    os.makedirs(OUT, exist_ok=True)
    gate = json.load(open(results_path(os.path.join("plan", "phase3", "gate_phase4.json"))))
    print("Phase-3 gate pass:", gate["pass"])
    t3 = pd.read_csv(results_path(os.path.join("plan", "phase3", "table3_1_optimum.csv")))
    # constants from the Fase-3 maxima with the same r3 and CI as here: the crest
    # constant depends on r3 (see docs), so pooling both r3 values would bias it
    real = t3[(t3.CI == "B") & (t3.md == 0) & (t3.mt > 0) & t3.real_max & (t3.r3 == R3)]
    c = float(real.rule_ratio.median())
    c_k = float(real.rule_ratio_kmean.median())
    print(f"Fase-3 constant c_k (mean degree) = {c_k:.3f}")
    print(f"Fase-3 constant c = median(D_t* m_t k_s/|r3|) = {c:.3f} over {len(real)} real maxima")
    jobs = [(s, ell, dt) for s, _, _ in GRAPHS for ell in ELLS for dt in DT_GRID]
    with Pool() as p:
        rows = p.map(run_one, jobs, chunksize=4)
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(OUT, "phase4_grid.csv"), index=False)
    summ = []
    for (seed, ell), g in d.groupby(["Seed", "ell"]):
        o = optimum(g, "E_T")
        gs = g.sort_values("Dt").reset_index(drop=True)
        i = o["i_star"]
        k = int(gs.k_s[0])
        kbar = float(np.mean([kk for _, kk in next(g_ for s_, g_, _ in GRAPHS if s_ == seed).degree()]))
        summ.append({"Seed": seed, "ell": ell, "k_s": k, **o,
                     "real_max": bool(o["interior"] and o["prominence"] > 0.05),
                     "m_focus_mean": gs.m_focus_mean[0], "m_flux_weighted": gs.m_flux_weighted[i],
                     "pred_focus_mean": c * abs(R3) / (gs.m_focus_mean[0] * k),
                     "pred_flux_weighted": c * abs(R3) / (gs.m_flux_weighted[i] * k),
                     "k_mean": kbar, "m_all_edges_mean": gs.m_all_edges_mean[0],
                     "pred_mean_degree": c_k * abs(R3) / (gs.m_all_edges_mean[0] * kbar)})
    s = pd.DataFrame(summ)
    s["c_used"], s["c_k_used"] = c, c_k
    s.to_csv(os.path.join(OUT, "phase4_summary.csv"), index=False)
    ok = s[s.real_max]
    for col in ("pred_focus_mean", "pred_flux_weighted", "pred_mean_degree"):
        err = np.log10(ok.Dt_star / ok[col])
        print(f"{col}: n={len(ok)}, median |log10 error| = {np.median(np.abs(err)):.3f}, "
              f"bias = {np.median(err):+.3f}, corr(log obs, log pred) = "
              f"{np.corrcoef(np.log10(ok.Dt_star), np.log10(ok[col]))[0, 1]:.3f}")

    import matplotlib.pyplot as plt
    from plan_plots import SERIES
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True, sharey=True)
    for ax, col, lab in zip(axes, ("pred_focus_mean", "pred_flux_weighted", "pred_mean_degree"),
                            ("(1) $c|r_3|/(\\bar m_{foco} k_s)$",
                             "(2) $c|r_3|/(\\bar m_{flujo} k_s)$",
                             "(3) $c_k|r_3|/(\\bar m_{red} \\langle k \\rangle)$ (no está en el plan)")):
        for k, ell in enumerate(ELLS):
            h = s[s.ell == ell]
            ax.scatter(h[col][h.real_max], h.Dt_star[h.real_max], color=SERIES[k], s=30, label=f"$\\ell$={ell}")
            ax.scatter(h[col][~h.real_max], h.Dt_star[~h.real_max], facecolor="none", edgecolor=SERIES[k], s=30)
        lim = [0.02, 5]
        ax.plot(lim, lim, "k-", lw=1)
        ax.axhline(DT_GRID[-1], color="0.6", ls=":", lw=1)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("$D_t^*$ predicho"); ax.set_title(lab)
    axes[0].set_ylabel("$D_t^*$ observado")
    axes[0].legend()
    fig.suptitle(f"Figura 4.1. Mortalidad por arista en grafos geométricos (CI-B, $D_d$={DD}, $c$={c:.2f}, $c_k$={c_k:.2f}). "
                 "Relleno: máximo interior con prominencia > 5 %; vacío: no. Línea punteada: borde de la malla.",
                 fontsize=9.5)
    fig.savefig(os.path.join(OUT, "fig4_1_observed_vs_predicted.png"))
    print(s[["Seed", "ell", "k_s", "Dt_star", "real_max", "prominence", "m_focus_mean", "m_flux_weighted",
             "pred_focus_mean", "pred_flux_weighted", "pred_mean_degree"]].round(3).to_string(index=False))
