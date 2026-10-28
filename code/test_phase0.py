"""Fase 0: tests that must pass before the sweeps of phases 1-3 are run.

  1. With m_d = m_t = 0 the modified integrator reproduces the PestAUC of the
     336 runs of experiment2_outbreak_containment.csv (relative error < 1e-8).
  2. Mass balance: with the local dynamics switched off, the total population of
     a species decays at rate D m sum_i k_i x_i.
  3. With D = 0, any value of m gives the same result as m = 0.
  4. With m = 0, the star and the complete graph give the same AUC (proposition
     of the manuscript). With m_t > 0 they stop coinciding; the difference is
     recorded.
Additional checks (not in the plan, needed by the sweeps):
  5. n_rel(T) = ceil(t_max/T) - 1 equals the number of releases actually applied.
  6. The sweep tolerance (rtol 1e-8) changes E_T and E_U by < 1e-5 (relative)
     with respect to the reference tolerance (rtol 1e-10).
  7. Edge-dependent mortality (Fase 4): a matrix m_ij equal to a constant m on
     every entry gives the same solution as the scalar m.

Run:  python test_phase0.py      (about 4 min on 2 cores)
Output: Results/plan/phase0/phase0_tests.json and the printed report.
"""

import json
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd

import core
from common import PARAMS, base_outbreak_ics, outbreak_ics, results_path
from core import n_releases, release_times
from network_experiments import (adjacency, experiment_metrics, laplacian, make_topology,
                                 network_rhs, solve_network_data)
import plan_common as pc

OUT = results_path(os.path.join("plan", "phase0"))
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------- test 1
def _exp2_row(args):
    topo, Dd, Dt, strategy = args
    core.RTOL, core.ATOL = 1e-10, 1e-12
    n, tt, T, delta = 8, 100, 15, 35
    st = {"None": (0, []), "Uniform": (delta, "All"), "TargetedAtSource": (delta, [0])}[strategy]
    data = solve_network_data(make_topology(topo, n), Dd, Dt, st[0], T, PARAMS,
                              base_outbreak_ics(n, tt), st[1], md=0.0, mt=0.0)
    m = experiment_metrics(data)
    return {"Topology": topo, "Dd": Dd, "Dt": Dt, "Strategy": strategy, **m}


def test1():
    ref = pd.read_csv(results_path("experiment2_outbreak_containment.csv"), keep_default_na=False)
    args = [(r.Topology, r.Dd, r.Dt, r.Strategy) for r in ref.itertuples()]
    with Pool() as p:
        new = pd.DataFrame(p.map(_exp2_row, args, chunksize=4))
    rel = np.abs(new.PestAUC.values - ref.PestAUC.values) / np.abs(ref.PestAUC.values)
    others = {}
    for c in ["PestPeak", "PestMin", "PestLastCycleMean", "PestWindowMean", "PestSpatialCV",
              "PestSync", "ParasitoidLastCycleMean", "ParasitoidWindowMean", "ParasitoidSync"]:
        a, b = new[c].values, ref[c].values.astype(float)
        others[c] = float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))
    return {"runs": int(len(ref)), "max_rel_err_PestAUC": float(rel.max()),
            "max_rel_err_other_metrics": others, "pass": bool(rel.max() < 1e-8)}


# ---------------------------------------------------------------- test 2
def test2():
    rng = np.random.default_rng(0)
    out, ok = [], True
    for name in ["Path", "Cycle", "Random_s4", "Star", "Complete"]:
        g = pc.network(name)
        L, A = laplacian(g), adjacency(g)
        k = A.sum(axis=1)
        n = len(k)
        for D, m in [(0.05, 0.1), (0.5, 0.3), (1.0, 1.0)]:
            # (a) algebraic identity on random states, both species
            f = network_rhs(L, D, D, PARAMS, md=m, mt=m, A=A, local=False)
            y = rng.uniform(0, 100, 4 * n)
            dy = f(0.0, y)
            alg = max(abs(dy[s * n:(s + 1) * n].sum() + D * m * k @ y[s * n:(s + 1) * n])
                      / (D * m * k @ y[s * n:(s + 1) * n]) for s in (0, 2))
            # (b) integrated: S(t) - S(0) = -D m int k.x dt, with releases off
            data = solve_network_data(g, D, D, 0, 15, PARAMS, outbreak_ics(n, 50, "B"), [],
                                      md=m, mt=m, local=False)
            sol = data["Solution"]
            ts = np.linspace(0, 50, 11)
            worst = 0.0
            for s in (0, 2):
                S = sol(ts)[s * n:(s + 1) * n].sum(axis=0)
                for t in ts[1:]:
                    integral = sol.integrate(0.0, t)[s * n:(s + 1) * n] @ k
                    lhs, rhs = S[ts == t][0] - S[0], -D * m * integral
                    worst = max(worst, abs(lhs - rhs) / max(abs(rhs), 1e-12))
            passed = alg < 1e-12 and worst < 1e-7
            ok &= passed
            out.append({"Network": name, "D": D, "m": m, "algebraic_rel_err": float(alg),
                        "integrated_rel_err": float(worst), "pass": bool(passed)})
    return {"cases": out, "pass": bool(ok)}


# ---------------------------------------------------------------- test 3
def test3():
    out, ok = [], True
    for name in ["Path", "Star", "Complete"]:
        g = pc.network(name)
        for ci in ("A", "B"):
            ics = outbreak_ics(8, 100, ci)
            base = solve_network_data(g, 0, 0, 35, 15, PARAMS, ics, [0]).copy()
            b = base["Solution"].integrate(0, 100)
            for m in (0.1, 0.5, 1.0):
                d = solve_network_data(g, 0, 0, 35, 15, PARAMS, ics, [0], md=m, mt=m)
                diff = float(np.max(np.abs(d["Solution"].integrate(0, 100) - b)))
                ok &= diff == 0.0
                out.append({"Network": name, "CI": ci, "m": m, "max_abs_diff_AUC": diff})
    return {"cases": out, "pass": bool(ok)}


# ---------------------------------------------------------------- test 4
def _pair(ci, D, strat, mt, md=0.0):
    nodes = {"None": [], "Uniform": "All", "Targeted": [0]}[strat]
    delta = 0 if strat == "None" else 35
    aucs = []
    for name in ("Star", "Complete"):
        d = solve_network_data(pc.network(name), D, D, delta, 15, PARAMS, outbreak_ics(8, 100, ci),
                               nodes, md=md, mt=mt)
        aucs.append(d["Solution"].integrate(0, 100)[:8].sum())
    return aucs


def test4():
    eq, ok = [], True
    for ci in ("A", "B"):
        for D in (0.01, 0.1, 1.0):
            for strat in ("None", "Uniform", "Targeted"):
                s, c = _pair(ci, D, strat, 0.0)
                rel = abs(s - c) / c
                ok &= rel < 1e-8
                eq.append({"CI": ci, "D": D, "Strategy": strat, "AUC_Star": s, "AUC_Complete": c,
                           "rel_diff": rel})
    diff = []
    for ci in ("A", "B"):
        for mt in (0.05, 0.1, 0.2, 0.4):
            for strat in ("None", "Targeted"):
                s, c = _pair(ci, 0.1, strat, mt)
                diff.append({"CI": ci, "D": 0.1, "mt": mt, "Strategy": strat, "AUC_Star": s,
                             "AUC_Complete": c, "rel_diff_pct": 100 * (s - c) / c})
    return {"m0_equal": eq, "mt_positive": diff, "pass": bool(ok)}


# ---------------------------------------------------------------- tests 5, 6
def test5():
    rows = []
    for T in (5, 7, 10, 15, 20, 30, 40):
        for tt in (100, 300):
            rows.append({"T": T, "tmax": tt, "n_rel": n_releases(T, tt),
                         "applied": len(release_times(T, tt)), "floor": int(np.floor(tt / T))})
    ok = all(r["n_rel"] == r["applied"] for r in rows)
    return {"cases": rows, "pass": bool(ok)}


def test6():
    out = []
    for name, ci, D in [("Complete", "B", 0.001), ("Complete", "B", 1.0), ("Path", "A", 0.03),
                        ("Random_s4", "B", 0.1)]:
        spec = {"Network": name, "CI": ci, "Dd": D, "Dt": D, "r3": -0.3, "tmax": 100}
        lo = pc.run_combo(dict(spec))
        hi = pc.run_combo(dict(spec, rtol=1e-10, atol=1e-12))
        out.append({"Network": name, "CI": ci, "D": D,
                    "rel_err_E_T": abs(lo["E_T"] - hi["E_T"]) / abs(hi["E_T"]),
                    "rel_err_E_U": abs(lo["E_U"] - hi["E_U"]) / abs(hi["E_U"]),
                    "seconds_sweep_tol": lo["Seconds"], "seconds_ref_tol": hi["Seconds"]})
    ok = all(o["rel_err_E_T"] < 1e-5 and o["rel_err_E_U"] < 1e-5 for o in out)
    return {"cases": out, "pass": bool(ok)}


def test7():
    out, ok = [], True
    for name in ["Path", "Random_s4", "Complete"]:
        g = pc.network(name)
        for m in (0.1, 0.4):
            ics = outbreak_ics(8, 100, "B")
            a = solve_network_data(g, 0.1, 0.3, 35, 15, PARAMS, ics, [0], md=m, mt=m)
            b = solve_network_data(g, 0.1, 0.3, 35, 15, PARAMS, ics, [0], md=np.full((8, 8), m),
                                   mt=np.full((8, 8), m))
            ia, ib = a["Solution"].integrate(0, 100), b["Solution"].integrate(0, 100)
            rel = float(np.max(np.abs(ia - ib) / np.abs(ia)))
            ok &= rel < 1e-10
            out.append({"Network": name, "m": m, "max_rel_diff": rel})
    return {"cases": out, "pass": bool(ok)}


if __name__ == "__main__":
    report = {}
    for name, fn in [("test1_reproduce_exp2", test1), ("test2_mass_balance", test2),
                     ("test3_D0_any_m", test3), ("test4_star_vs_complete", test4),
                     ("test5_n_rel", test5), ("test6_sweep_tolerance", test6),
                     ("test7_edge_matrix", test7)]:
        report[name] = fn()
        print(f"{name}: {'PASS' if report[name]['pass'] else 'FAIL'}", flush=True)
    t1 = report["test1_reproduce_exp2"]
    print(f"  test1: {t1['runs']} runs, max rel err PestAUC = {t1['max_rel_err_PestAUC']:.2e}")
    print("  test2: worst algebraic / integrated rel err =",
          f"{max(c['algebraic_rel_err'] for c in report['test2_mass_balance']['cases']):.1e} /",
          f"{max(c['integrated_rel_err'] for c in report['test2_mass_balance']['cases']):.1e}")
    print("  test4 (m_t > 0, D = 0.1):")
    for r in report["test4_star_vs_complete"]["mt_positive"]:
        print(f"    CI-{r['CI']} mt={r['mt']:<5} {r['Strategy']:<9} star - complete = {r['rel_diff_pct']:+.3f} %")
    print("  test6:", [(c["Network"], c["D"], f"{c['rel_err_E_T']:.1e}") for c in report["test6_sweep_tolerance"]["cases"]])
    with open(os.path.join(OUT, "phase0_tests.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=float)
    all_ok = all(v["pass"] for v in report.values())
    print("ALL PHASE-0 TESTS PASS" if all_ok else "SOME PHASE-0 TEST FAILED")
    raise SystemExit(0 if all_ok else 1)
