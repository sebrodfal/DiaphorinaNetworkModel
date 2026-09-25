"""Network experiment infrastructure for DiaphorinaModel (port of NetworkExperiments.wl).

Reuses the ODE structure of system_evolve_network (functions.py) but
decouples solving from plotting, so batches of experiments can run headlessly
and results can be exported as data.

Node i of a Mathematica graph (1-based) is node i-1 here (0-based).
"""

import numpy as np
import networkx as nx

from core import local_field, release_times, solve_impulsive, subdivide

# Variable handles, used where Mathematica passes the symbols x1..x4.
x1, x2, x3, x4 = 0, 1, 2, 3


# ---------- Solver (Model + Solver, no plotting) ----------

def laplacian(graph):
    """Normal @ KirchhoffMatrix[graph], rows/columns in vertex order."""
    return nx.laplacian_matrix(graph, nodelist=sorted(graph.nodes())).toarray().astype(float)


def network_rhs(L, Dd, Dt, params):
    n = L.shape[0]

    def f(t, y):
        y1, y2, y3, y4 = y[:n], y[n:2 * n], y[2 * n:3 * n], y[3 * n:]
        f1, f2, f3, f4 = local_field(params, y1, y2, y3, y4)
        return np.concatenate((f1 - Dd * (L @ y1), f2, f3 - Dt * (L @ y3), f4))

    return f


def solve_network_data(graph, Dd, Dt, delta, T, params, init, release_nodes="All"):
    """SolveNetworkData[graph, Dd, Dt, delta, T, params, {Di0,B0,x30,x40,tt}, releaseNodes].

    release_nodes: "All" or a list of 0-based node indices ([] = no releases)."""
    Di0, B0, x30, x40, tt = init
    L = laplacian(graph)
    n = graph.number_of_nodes()
    y0 = np.concatenate([np.asarray(v, dtype=float) for v in (Di0, B0, x30, x40)])
    targets = list(range(n)) if release_nodes == "All" else list(release_nodes)

    if delta == 0 or len(targets) == 0:
        jumps, jump = [], None
    else:
        idx = 2 * n + np.array(targets)

        def jump(y):
            y = y.copy()
            y[idx] += delta
            return y
        jumps = release_times(T, tt)

    sol = solve_impulsive(network_rhs(L, Dd, Dt, params), y0, (0.0, float(tt)), jumps, jump)
    return {"Solution": sol, "Graph": graph, "N": n, "Dd": Dd, "Dt": Dt,
            "Delta": delta, "T": T, "tt": tt, "ReleaseNodes": targets}


# ---------- Topology generators ----------

def make_topology(name, n, seed=1):
    """MakeTopology[name, n, seed].

    Deterministic topologies match Mathematica exactly (StarGraph's hub is
    vertex 1 -> node 0 here). Stochastic ones are drawn from the same
    distributions with numpy's RNG, so the particular graph instances differ
    from Mathematica's for the same seed.
    WattsStrogatzGraphDistribution[n,p,k] starts from a 2k-regular ring, which
    is networkx's watts_strogatz_graph(n, 2k, p)."""
    if name == "Path":
        return nx.path_graph(n)
    if name == "Cycle":
        return nx.cycle_graph(n)
    if name == "Star":
        return nx.star_graph(n - 1)
    if name == "Complete":
        return nx.complete_graph(n)
    if name == "Random":
        return nx.gnp_random_graph(n, 0.4, seed=seed)
    if name == "SmallWorld":
        k = max(1, min(2, n // 2 - 1))
        return nx.watts_strogatz_graph(n, 2 * k, 0.2, seed=seed)
    if name == "ScaleFree":
        return nx.barabasi_albert_graph(n, min(2, n - 1), seed=seed)
    raise ValueError(f"UnknownTopology: {name}")


# ---------- Metrics ----------

def sample_times(data, npts=400):
    return subdivide(0.0, data["tt"], npts)


def _row(data, var, i):
    return var * data["N"] + i


def node_if(data, var, i):
    """NodeIF: the interpolating function var[i](t)."""
    sol, k = data["Solution"], _row(data, var, i)
    return lambda t: sol(t)[k]


def auc(data, var, i):
    return data["Solution"].integrate(0.0, data["tt"], _row(data, var, i))


def total_auc(data, var):
    return sum(auc(data, var, i) for i in range(data["N"]))


def last_cycle_mean(data, var, i):
    T, tt = data["T"], data["tt"]
    if T <= 0 or T >= tt:
        return node_if(data, var, i)(tt)
    t0 = tt - np.mod(tt, T) - T
    if t0 < 0:
        t0 = 0
    return data["Solution"].integrate(t0, tt, _row(data, var, i)) / (tt - t0)


def window_mean(data, var, i, win):
    """Mean over an explicit [t0,t1] window, independent of the release period T.
    Matches the single-node paper's reporting convention (mean over [tt/2,tt])."""
    return data["Solution"].integrate(win[0], win[1], _row(data, var, i)) / (win[1] - win[0])


def peak_value(data, var, i):
    return float(np.max(node_if(data, var, i)(sample_times(data))))


def node_min_value(data, var, i):
    return float(np.min(node_if(data, var, i)(sample_times(data))))


def spatial_cv(data, var):
    vals = np.array([last_cycle_mean(data, var, i) for i in range(data["N"])])
    if np.mean(vals) == 0:
        return 0.0
    return float(np.std(vals, ddof=1) / np.mean(vals))


def sync_index(data, var):
    n = data["N"]
    if n < 2:
        return 1.0
    tv = sample_times(data)
    series = data["Solution"](tv)[var * n:(var + 1) * n]
    C = np.corrcoef(series)
    iu = np.triu_indices(n, 1)
    return float(np.mean(C[iu]))


def algebraic_connectivity(graph):
    ev = np.sort(np.linalg.eigvalsh(laplacian(graph)))
    return 0.0 if len(ev) < 2 else float(ev[1])


def experiment_metrics(data, window=None):
    """ExperimentMetrics. window defaults to (tt/2, tt), matching ReleaseSweep/
    ScenarioMetrics in the single-node code."""
    win = (data["tt"] / 2, data["tt"]) if window is None else window
    N = data["N"]
    g = data["Graph"]
    return {
        "PestPeak": max(peak_value(data, x1, i) for i in range(N)),
        "PestMin": min(node_min_value(data, x1, i) for i in range(N)),
        "PestAUC": float(total_auc(data, x1)),
        "PestLastCycleMean": float(np.mean([last_cycle_mean(data, x1, i) for i in range(N)])),
        "PestWindowMean": float(np.mean([window_mean(data, x1, i, win) for i in range(N)])),
        "PestSpatialCV": spatial_cv(data, x1),
        "PestSync": sync_index(data, x1),
        "ParasitoidLastCycleMean": float(np.mean([last_cycle_mean(data, x3, i) for i in range(N)])),
        "ParasitoidWindowMean": float(np.mean([window_mean(data, x3, i, win) for i in range(N)])),
        "ParasitoidSync": sync_index(data, x3),
        "AlgebraicConnectivity": algebraic_connectivity(g),
        "MeanDegree": float(np.mean([d for _, d in g.degree()])),
    }
