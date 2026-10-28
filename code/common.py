"""Base parameters (tabla_latex.txt) and CSV export in the same style as
Mathematica's Export[..., "CSV"]: strings quoted, booleans as true/false,
reals with full precision."""

import os

import numpy as np

P1 = [-0.20, 0.040, -0.010, 5.55e-4, 0.005]
P2 = [-0.72, -0.010, 0.054, 5.55e-4, 0.005]
P3 = [-0.30, 0.010, 5.55e-4, 0.005]
P4 = [0.72, -0.072, 5.55e-4, 0.005]
PARAMS = [P1, P2, P3, P4]

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "Results")


def results_path(name):
    os.makedirs(RESULTS, exist_ok=True)
    return os.path.join(RESULTS, name)


def _fmt(v):
    if v is None:
        return "Missing[]"
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    return v


def write_csv(path, rows):
    """rows: list of dicts sharing the keys of rows[0] (column order = key order)."""
    cols = list(rows[0].keys())

    def cell(v):
        f = _fmt(v)
        if isinstance(v, str):
            return '"' + v.replace('"', '""') + '"'
        return f

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(",".join('"' + c + '"' for c in cols) + "\n")
        for r in rows:
            fh.write(",".join(cell(r[c]) for c in cols) + "\n")


def outbreak_ics(n, tt, kind="A"):
    """Initial conditions (plan, Fase 0, cambio 2). Node 0 is the focus.
    CI-A ("A"): 200 psyllids at the focus and 20 at every other node (original).
    CI-B ("B"): point outbreak, 200 at the focus and 0 elsewhere; the pest can
                only reach the other nodes by dispersal, so Dd = 0 makes no sense.
    Both use x2 = 10, x3 = 15, x4 = 30 at every node."""
    other = {"A": 20, "B": 0}[kind]
    return [[200] + [other] * (n - 1), [10] * n, [15] * n, [30] * n, tt]


def base_outbreak_ics(n, tt):
    """Point-source outbreak: node 1 starts at 200, the rest at 20 (CI-A)."""
    return outbreak_ics(n, tt, "A")
