"""Fase 1: relative dispersion without transit cost (m_d = m_t = 0).

Grid: D_t, D_d in 9 log-spaced values in [1e-3, 1]; r3 in {-0.2, -0.3, -0.45};
7 focal networks (focus degree 1..7); CI-A and CI-B; strategies None, Uniform,
Targeted (the three are solved inside each combination).
81 x 3 x 7 x 2 = 3402 combinations = 10 206 runs.

Pilot (run first with --pilot): complete graph, CI-B, r3 = -0.3, 81 combos = 243 runs.

Output: Results/plan/phase1/phase1_pilot.csv, Results/plan/phase1/phase1_grid.csv
"""

import os
import sys

from common import results_path
import plan_common as pc

D_GRID = pc.logspace_grid(9)
R3_GRID = [-0.2, -0.3, -0.45]


def specs(pilot=False):
    if pilot:
        return [{"Network": "Complete", "CI": "B", "Dd": dd, "Dt": dt, "r3": -0.3, "tmax": 100}
                for dd in D_GRID for dt in D_GRID]
    return [{"Network": net, "CI": ci, "Dd": dd, "Dt": dt, "r3": r3, "tmax": 100}
            for net in pc.network_names() for ci in ("A", "B") for r3 in R3_GRID
            for dd in D_GRID for dt in D_GRID]


if __name__ == "__main__":
    pilot = "--pilot" in sys.argv
    name = "phase1_pilot.csv" if pilot else "phase1_grid.csv"
    pc.run_specs(specs(pilot), results_path(os.path.join("plan", "phase1", name)))
