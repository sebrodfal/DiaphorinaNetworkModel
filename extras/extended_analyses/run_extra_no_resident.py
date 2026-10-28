"""Extra (after the 2026-10-17 review): releases into a landscape WITHOUT a resident
parasitoid population.

In the plan every node starts with x3 = 15. With transit mortality that resident
population collapses in the no-release run as D_t grows, so the Fase-3 crest of
E_T is an optimum relative to a deteriorating counterfactual. Here x3(0) = 0 at
every node (Tamarixia not yet established; pest as in CI-B: 200 at the focus, 0
elsewhere). The no-release run then has no parasitoid at all and does not depend
on D_t, so a maximum of E_T is a genuine minimum of the pest AUC.

Grid: D_t in 25 log-spaced values in [1e-3, 10]; m_t in {0, 0.05, 0.1, 0.2, 0.4};
D_d in {0.01, 0.03, 0.1}; 7 focal networks + complete graph; r3 in {-0.3, -0.45};
m_d = 0, delta = 35, T = 15, t_max = 100. Strategies None, Uniform, Targeted.
25 x 5 x 3 x 8 x 2 = 6000 combinations = 18 000 runs.

Output: Results/no_resident_grid.csv next to this script (resumable). Same columns
as the phase grids, with CI = "B0". The model code is imported from code/.
"""

import csv
import os
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))

import core
from core import n_releases
from network_experiments import solve_network_data
import plan_common as pc

DT_GRID = pc.logspace_grid(25, 1e-3, 10.0)
MT_GRID = [0.0, 0.05, 0.1, 0.2, 0.4]
DD_GRID = [0.01, 0.03, 0.1]
R3_GRID = [-0.3, -0.45]
N = pc.N_NODES
CSV = os.path.join(HERE, "Results", "no_resident_grid.csv")


def ics_no_resident(tt):
    return [[200] + [0] * (N - 1), [10] * N, [0] * N, [30] * N, tt]


def run_combo(spec):
    core.RTOL, core.ATOL = pc.SWEEP_RTOL, pc.SWEEP_ATOL
    t0 = time.time()
    g = pc.network(spec["Network"])
    params = pc.params_with_r3(spec["r3"])
    res = {}
    for strat, nodes, d in (("None", [], 0), ("Uniform", "All", pc.DELTA), ("Targeted", [0], pc.DELTA)):
        data = solve_network_data(g, spec["Dd"], spec["Dt"], d, pc.T_REL, params,
                                  ics_no_resident(spec["tmax"]), nodes, mt=spec["mt"])
        v = data["Solution"].integrate(0.0, spec["tmax"])
        res[strat] = (float(v[:N].sum()), float(v[0]), float(v[2 * N:3 * N].sum()))
    k_s = pc.focus_degree(spec["Network"])
    nrel = n_releases(pc.T_REL, spec["tmax"])
    QT, QU = pc.DELTA * nrel, pc.DELTA * nrel * N
    ET = (res["None"][0] - res["Targeted"][0]) / QT
    EU = (res["None"][0] - res["Uniform"][0]) / QU
    row = {k: spec[k] for k in ("Network", "CI", "Dd", "Dt", "r3", "tmax")}
    row.update({"md": 0.0, "mt": spec["mt"], "Delta": pc.DELTA, "T": pc.T_REL, "k_s": k_s,
                "n_rel": nrel, "Q_T": QT, "Q_U": QU,
                "AUC_None": res["None"][0], "AUC_Uniform": res["Uniform"][0],
                "AUC_Targeted": res["Targeted"][0],
                "AUCfocus_None": res["None"][1], "AUCfocus_Targeted": res["Targeted"][1],
                "ParAUC_None": res["None"][2], "ParAUC_Targeted": res["Targeted"][2],
                "E_T": ET, "E_U": EU, "E_ratio": ET / EU if EU != 0 else float("nan"),
                "Pi_t": pc.pi_t(spec["Dt"], k_s, spec["r3"], spec["mt"]),
                "Seconds": time.time() - t0})
    return row


def specs():
    return [{"Network": net, "CI": "B0", "Dd": dd, "Dt": dt, "r3": r3, "tmax": 100, "mt": mt}
            for net in pc.network_names(include_complete=True) for r3 in R3_GRID
            for mt in MT_GRID for dd in DD_GRID for dt in DT_GRID]


def main():
    done = set()
    if os.path.exists(CSV):
        with open(CSV, newline="") as fh:
            done = {pc.spec_key(r) for r in csv.DictReader(fh)}
    todo = [s for s in specs() if pc.spec_key(s) not in done]
    print(f"{CSV}: {len(done)} done, {len(todo)} to run", flush=True)
    if not todo:
        return
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    new_file = not os.path.exists(CSV) or os.path.getsize(CSV) == 0
    t0 = time.time()
    with open(CSV, "a", newline="") as fh, Pool() as pool:
        writer = None
        for i, row in enumerate(pool.imap_unordered(run_combo, todo, chunksize=4), 1):
            if writer is None:
                writer = csv.DictWriter(fh, fieldnames=list(row.keys()))
                if new_file:
                    writer.writeheader()
            writer.writerow(row)
            if i % 200 == 0:
                fh.flush()
                el = time.time() - t0
                print(f"  {i}/{len(todo)}  {el:.0f}s, ~{el / i * (len(todo) - i):.0f}s left", flush=True)


if __name__ == "__main__":
    main()
