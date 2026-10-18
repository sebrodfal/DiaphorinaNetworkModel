"""Fase 2: Goldilocks hypothesis without dispersal cost (m_d = m_t = 0).

Grid: D_t in 25 log-spaced values in [1e-3, 1]; D_d in {0.01, 0.03, 0.1};
7 focal networks; CI-A and CI-B; t_max in {100, 300}; strategies None,
Uniform and Targeted inside each combination (r3 = -0.3).
25 x 3 x 7 x 2 x 2 = 2100 combinations = 6300 runs.

Verification run (--check-dd0): D_d = 0 with CI-A, same 25 D_t values, both
horizons, 7 networks (350 combinations). It tests whether the drop of E_T from
5.30 to 4.91 in the D_d = 0 row of Experiment 2 is real or a coarse-grid effect.

Output: Results/plan/phase2/phase2_grid.csv, Results/plan/phase2/phase2_dd0_check.csv
"""

import os
import sys

from common import results_path
import plan_common as pc

DT_GRID = pc.logspace_grid(25)
DD_GRID = [0.01, 0.03, 0.1]
TMAX = [100, 300]


def specs():
    return [{"Network": net, "CI": ci, "Dd": dd, "Dt": dt, "r3": -0.3, "tmax": tm}
            for tm in TMAX for net in pc.network_names() for ci in ("A", "B")
            for dd in DD_GRID for dt in DT_GRID]


def specs_dd0():
    return [{"Network": net, "CI": "A", "Dd": 0.0, "Dt": dt, "r3": -0.3, "tmax": tm}
            for tm in TMAX for net in pc.network_names() for dt in [0.0] + DT_GRID]


if __name__ == "__main__":
    if "--check-dd0" in sys.argv:
        pc.run_specs(specs_dd0(), results_path(os.path.join("plan", "phase2", "phase2_dd0_check.csv")))
    else:
        pc.run_specs(specs(), results_path(os.path.join("plan", "phase2", "phase2_grid.csv")))
