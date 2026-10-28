"""Shared definitions for the work plan "dispersión relativa, hipótesis Goldilocks
y mortalidad en tránsito" (phases 0-4).

Everything the phase scripts need lives here: the set of focal networks, the
dimensionless groups (Pi, Pi_t), the efficiency of a release and a resumable
parallel runner that writes one CSV row per parameter combination.

Notation (as in the plan)
-------------------------
m_d, m_t   transit mortality of the psyllid (x1) and the parasitoid (x3)
n_rel(T)   number of releases actually applied, ceil(t_max/T) - 1
Q          parasitoids released = delta * n_rel * |R|
E          (AUC_none - AUC_strategy) / Q, AUC = network-wide integral of x1
E_T, E_U   efficiency of the release targeted at the focus / uniform release
tau_t      1/|r3|;  tau_d = 1/|r1| = 5 or T0 = 7.8
Pi         (D_t tau_t) / (D_d tau_d)
Pi_t       D_t k_s / (|r3| + m_t D_t k_s), k_s = degree of the focus (node 0)
"""

import copy
import csv
import os
import time
from multiprocessing import Pool

import numpy as np
import networkx as nx

import core
from common import PARAMS, outbreak_ics
from core import n_releases
from network_experiments import make_topology, solve_network_data

N_NODES = 8
DELTA = 35
T_REL = 15
T_MAX = 100
TAU_D = {"r1": 1.0 / abs(PARAMS[0][0]),  # 5.0
         "T0": 7.8}                       # intrinsic period quoted in the plan
STRATEGIES = ("None", "Uniform", "Targeted")

# Integrator tolerances for the sweeps. The 336-run reproduction test of Fase 0
# uses the original tolerances (core.RTOL=1e-10, core.ATOL=1e-12). The sweeps use
# 1e-8/1e-10; tests/phase0 checks that this changes E_T by < 1e-5 (relative).
SWEEP_RTOL = 1e-8
SWEEP_ATOL = 1e-10


# ---------------------------------------------------------------------------
# Focal networks: 8 nodes, focus (node 0) with degree 1..7
# ---------------------------------------------------------------------------
# The four stochastic instances are graphs of Experiment 3 (same generator and
# seed). They are the first seeds (1..10) whose node 0 has degree 3, 4, 5 and 6;
# no Random instance reaches degree 6, so that one is a ScaleFree graph.
FOCAL_NETWORKS = [
    # (name, generator, seed)
    ("Path", "Path", 1),             # focus at one end, degree 1
    ("Cycle", "Cycle", 1),           # degree 2
    ("Random_s9", "Random", 9),      # degree 3
    ("Random_s3", "Random", 3),      # degree 4
    ("Random_s4", "Random", 4),      # degree 5
    ("ScaleFree_s2", "ScaleFree", 2),  # degree 6
    ("Star", "Star", 1),             # hub, degree 7
]
# The complete graph also has a focus of degree 7. The plan uses it for the
# pilots and for Fase 3c, but it is not one of the seven networks. With m = 0 it
# is dynamically identical to the star (proposition of the manuscript, checked
# in Fase 0); with m_t > 0 it is not, so it is added as an eighth network there.
COMPLETE = ("Complete", "Complete", 1)

_GRAPH_CACHE = {}


def network(name):
    if name not in _GRAPH_CACHE:
        spec = {n: (g, s) for n, g, s in FOCAL_NETWORKS + [COMPLETE]}[name]
        g = make_topology(spec[0], N_NODES, spec[1])
        assert nx.is_connected(g)
        _GRAPH_CACHE[name] = g
    return _GRAPH_CACHE[name]


def focus_degree(name):
    return int(network(name).degree(0))


def network_names(include_complete=False):
    names = [n for n, _, _ in FOCAL_NETWORKS]
    return names + ["Complete"] if include_complete else names


# ---------------------------------------------------------------------------
# Dimensionless groups
# ---------------------------------------------------------------------------

def pi_relative(Dd, Dt, r3, tau_d):
    return (Dt / abs(r3)) / (Dd * tau_d)


def pi_t(Dt, k_s, r3, mt=0.0):
    return Dt * k_s / (abs(r3) + mt * Dt * k_s)


def logspace_grid(n, lo=1e-3, hi=1.0):
    return list(np.logspace(np.log10(lo), np.log10(hi), n))


# ---------------------------------------------------------------------------
# One parameter combination = the three strategies under identical parameters
# ---------------------------------------------------------------------------

def params_with_r3(r3):
    p = copy.deepcopy(PARAMS)
    p[2][0] = r3
    return p


def _solve(spec, release_nodes, delta):
    g = network(spec["Network"])
    ics = outbreak_ics(N_NODES, spec["tmax"], spec["CI"])
    return solve_network_data(g, spec["Dd"], spec["Dt"], delta, spec.get("T", T_REL),
                              params_with_r3(spec["r3"]), ics, release_nodes,
                              md=spec.get("md", 0.0), mt=spec.get("mt", 0.0))


def _aucs(data):
    """Network pest AUC, focus pest AUC and network parasitoid AUC."""
    n, tt = data["N"], data["tt"]
    v = data["Solution"].integrate(0.0, tt)
    return float(v[:n].sum()), float(v[0]), float(v[2 * n:3 * n].sum())


def run_combo(spec):
    """spec keys: Network, Dd, Dt, r3, CI, tmax, md, mt (optional: T, Delta).
    Returns one flat row with the inputs, the three AUCs and all derived
    quantities (k_s, n_rel, Q_T, Q_U, E_T, E_U, E_T/E_U, Pi, Pi_t)."""
    core.RTOL, core.ATOL = spec.get("rtol", SWEEP_RTOL), spec.get("atol", SWEEP_ATOL)
    delta = spec.get("Delta", DELTA)
    T = spec.get("T", T_REL)
    t0 = time.time()
    res = {}
    for strat, nodes, d in (("None", [], 0), ("Uniform", "All", delta), ("Targeted", [0], delta)):
        res[strat] = _aucs(_solve(spec, nodes, d))
    k_s = focus_degree(spec["Network"])
    nrel = n_releases(T, spec["tmax"])
    QT, QU = delta * nrel * 1, delta * nrel * N_NODES
    ET = (res["None"][0] - res["Targeted"][0]) / QT
    EU = (res["None"][0] - res["Uniform"][0]) / QU
    mt = spec.get("mt", 0.0)
    row = {k: spec[k] for k in ("Network", "CI", "Dd", "Dt", "r3", "tmax")}
    row.update({"md": spec.get("md", 0.0), "mt": mt, "Delta": delta, "T": T,
                "k_s": k_s, "n_rel": nrel, "Q_T": QT, "Q_U": QU,
                "AUC_None": res["None"][0], "AUC_Uniform": res["Uniform"][0],
                "AUC_Targeted": res["Targeted"][0],
                "AUCfocus_None": res["None"][1], "AUCfocus_Targeted": res["Targeted"][1],
                "ParAUC_None": res["None"][2], "ParAUC_Targeted": res["Targeted"][2],
                "E_T": ET, "E_U": EU, "E_ratio": ET / EU if EU != 0 else np.nan,
                "Pi_r1": pi_relative(spec["Dd"], spec["Dt"], spec["r3"], TAU_D["r1"])
                if spec["Dd"] > 0 else np.inf,
                "Pi_T0": pi_relative(spec["Dd"], spec["Dt"], spec["r3"], TAU_D["T0"])
                if spec["Dd"] > 0 else np.inf,
                "Pi_t": pi_t(spec["Dt"], k_s, spec["r3"], mt),
                "Seconds": time.time() - t0})
    return row


# ---------------------------------------------------------------------------
# Resumable parallel runner
# ---------------------------------------------------------------------------

KEY_COLS = ("Network", "CI", "Dd", "Dt", "r3", "tmax", "md", "mt")


def spec_key(d):
    return tuple(str(d.get(k, 0.0)) if k not in ("Dd", "Dt", "r3", "md", "mt", "tmax")
                 else repr(float(d.get(k, 0.0))) for k in KEY_COLS)


def run_specs(specs, csv_path, processes=None, chunksize=4, verbose=True):
    """Runs every spec not yet present in csv_path and appends rows as they
    finish, so an interrupted sweep can be resumed by calling it again."""
    done = set()
    if os.path.exists(csv_path):
        with open(csv_path, newline="") as fh:
            for r in csv.DictReader(fh):
                done.add(spec_key(r))
    todo = [s for s in specs if spec_key(s) not in done]
    if verbose:
        print(f"{csv_path}: {len(specs)} combos, {len(done)} already done, {len(todo)} to run")
    if not todo:
        return
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    new_file = not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0
    t0 = time.time()
    with open(csv_path, "a", newline="") as fh, Pool(processes or os.cpu_count()) as pool:
        writer = None
        for i, row in enumerate(pool.imap_unordered(run_combo, todo, chunksize=chunksize), 1):
            if writer is None:
                writer = csv.DictWriter(fh, fieldnames=list(row.keys()))
                if new_file:
                    writer.writeheader()
            writer.writerow(row)
            if i % 50 == 0:
                fh.flush()
                if verbose:
                    el = time.time() - t0
                    print(f"  {i}/{len(todo)}  {el:.0f}s elapsed, ~{el / i * (len(todo) - i):.0f}s left",
                          flush=True)
