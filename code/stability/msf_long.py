"""Transverse stability of the synchronous state (section 2(h)).

Computes the largest Lyapunov exponent of every Laplacian mode along the
uniform-release single-patch trajectory, with standard errors from batch
means, for the seven seed-1 topologies and, on the two regular graphs, for
parasitoid transit mortality. Writes CSV files to the output directory.

Usage (from code/):
    python stability/msf_long.py --tend 30000          # writes Results/stability/ (~15 min)
    python stability/msf_long.py --tend 1200 --out Results/stability_short   # original horizon
"""
import argparse
import os

import networkx as nx
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.stats import spearmanr

# Local model (table 1)
R = np.array([-0.20, -0.72, -0.30, 0.72])
A_CROWD = np.full(4, 5.55e-4)
C_CROWD = np.full(4, 5e-3)
B = np.zeros((4, 4))
B[0, 1], B[0, 2] = 0.040, -0.010
B[1, 0], B[1, 3] = -0.010, 0.054
B[2, 0] = 0.010
B[3, 1] = -0.072

X0 = (50.0, 10.0, 15.0, 30.0)
DELTA, T = 35.0, 15.0
D = 0.03
DISCARD = 0.3
N_BATCH = 20

# Synchrony index of Experiment 2 (uniform releases, seed-1 graphs), table 3
SYNC = {"path": 0.731, "cycle": 0.744, "star": 0.851, "complete": 0.851,
        "random": 0.752, "small world": 0.823, "scale free": 0.820}


def f(x, r):
    beta = B @ x
    return x * (r + beta - (A_CROWD + C_CROWD * beta) * x)


def jac(x, r):
    beta = B @ x
    J = x[:, None] * B * (1.0 - C_CROWD * x)[:, None]
    J[np.diag_indices(4)] = r + beta - 2.0 * (A_CROWD + C_CROWD * beta) * x
    return J


def exponents(mus, gamma, r, tend):
    """Per-segment log growth rates for each mode in `mus`."""
    k = len(mus)
    x = np.array(X0, float)
    Y = np.tile([1.0, 0, 0, 0], (k, 1))
    t0 = DISCARD * tend
    logs = [[] for _ in range(k)]

    def rhs(t, z):
        xs, Ys = z[:4], z[4:].reshape(k, 4)
        Jx = jac(xs, r)
        dY = np.array([(Jx - m * gamma) @ Ys[i] for i, m in enumerate(mus)])
        return np.concatenate([f(xs, r), dY.ravel()])

    t = 0.0
    for tk in np.arange(T, tend + 1e-9, T):
        sol = solve_ivp(rhs, (t, tk), np.concatenate([x, Y.ravel()]),
                        method="DOP853", rtol=1e-10, atol=1e-12)
        z = sol.y[:, -1]
        x, Y = z[:4].copy(), z[4:].reshape(k, 4)
        if t >= t0 - 1e-9:
            nrm = np.linalg.norm(Y, axis=1)
            for i in range(k):
                logs[i].append(np.log(nrm[i]) / T)
            Y = Y / nrm[:, None]
        else:
            Y = np.tile([1.0, 0, 0, 0], (k, 1))
        if tk < tend - 1e-9:          # releases only for t_k < t_end
            x[2] += DELTA
        t = tk
    return [np.asarray(s) for s in logs]


def mean_se(s):
    nb = N_BATCH if len(s) >= 2 * N_BATCH else max(2, len(s) // 2)
    b = s[: len(s) // nb * nb].reshape(nb, -1).mean(1)
    return s.mean(), b.std(ddof=1) / np.sqrt(nb)


def graphs():
    return {
        "path": nx.path_graph(8),
        "cycle": nx.cycle_graph(8),
        "star": nx.star_graph(7),
        "complete": nx.complete_graph(8),
        "random": nx.gnp_random_graph(8, 0.4, seed=1),
        "small world": nx.watts_strogatz_graph(8, 4, 0.2, seed=1),
        "scale free": nx.barabasi_albert_graph(8, 2, seed=1),
    }


def laplacian_eigs(g):
    L = nx.laplacian_matrix(g, nodelist=sorted(g)).toarray().astype(float)
    return np.sort(np.linalg.eigvalsh(L))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tend", type=float, default=30000)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "..", "Results", "stability"))
    ap.add_argument("--mc", type=int, default=5000)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    G = graphs()
    eig = {n: laplacian_eigs(g) for n, g in G.items()}
    mus = sorted({0.0} | {round(float(e), 10) for e in np.concatenate([v[1:] for v in eig.values()])})
    gamma = np.diag([D, 0, D, 0])
    S = exponents(mus, gamma, R, args.tend)
    by_mu = {m: mean_se(s) for m, s in zip(mus, S)}

    pd.DataFrame([{"mu": m, "Lambda": v[0], "se": v[1], "n_segments": len(s)}
                  for (m, v), s in zip(by_mu.items(), S)]).to_csv(
        os.path.join(args.out, "modes.csv"), index=False)

    rows = []
    for n, e in eig.items():
        cand = [(by_mu[round(float(mu), 10)], mu) for mu in e[1:]]
        (lam, se), mu_max = max(cand, key=lambda c: c[0][0])
        rows.append({"topology": n, "lambda2": e[1], "mu_of_max": mu_max,
                     "Lambda_max": lam, "se": se,
                     "amplification_100": np.exp(100 * lam),
                     "efolding_time": 1 / abs(lam), "synchrony_index": SYNC[n]})
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(args.out, "table3_stability.csv"), index=False)

    rho, p = spearmanr(tab.Lambda_max, tab.synchrony_index)
    rng = np.random.default_rng(0)
    sims = [spearmanr(tab.Lambda_max + rng.normal(size=len(tab)) * tab.se,
                      tab.synchrony_index)[0] for _ in range(args.mc)]
    rho_l2, p_l2 = spearmanr(tab.lambda2, tab.synchrony_index)
    m1, se1 = by_mu[0.0]
    pd.DataFrame([{
        "spearman_Lambda_sync": rho, "p": p,
        "mc_2.5": np.percentile(sims, 2.5), "mc_97.5": np.percentile(sims, 97.5),
        "spearman_lambda2_sync": rho_l2, "p_lambda2": p_l2,
        "mode1_Lambda": m1, "mode1_se": se1,
    }]).to_csv(os.path.join(args.out, "summary.csv"), index=False)

    # Transit mortality on the regular graphs (section 2(h), last paragraph)
    mort = []
    for n, k in (("cycle", 2), ("complete", 7)):
        e = [m for m in eig[n][1:]]
        e = sorted(set(np.round(e, 10)))
        for mt in (0.0, 0.05, 0.1, 0.2, 0.4):
            r = R.copy()
            r[2] = R[2] - D * mt * k
            g_mt = np.diag([D, 0, (1 - mt) * D, 0])
            Sm = exponents(e, g_mt, r, args.tend)
            lam, se = max((mean_se(s) for s in Sm), key=lambda v: v[0])
            mort.append({"topology": n, "m_t": mt, "Lambda_max": lam, "se": se})
    pd.DataFrame(mort).to_csv(os.path.join(args.out, "mortality_regular.csv"), index=False)
    print(tab.round(5).to_string(index=False))
    print(f"rho={rho:.3f} p={p:.4f}  mode1={m1*1e3:.3f}e-3 +- {se1*1e3:.3f}e-3")


if __name__ == "__main__":
    main()
