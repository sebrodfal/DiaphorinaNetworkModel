"""Fase 3a: transit mortality of the parasitoid only (m_d = 0).

Grid: m_t in {0, 0.05, 0.1, 0.2, 0.4}; D_t in 13 log-spaced values in [1e-3, 1];
D_d in {0.01, 0.03, 0.1}; 7 focal networks plus the complete graph (needed by
Fase 3c, and no longer equivalent to the star once m_t > 0); CI-B (main) and
CI-A (control); r3 in {-0.3, -0.45}; strategies None, Uniform, Targeted.
5 x 13 x 3 x 8 x 2 x 2 = 6240 combinations = 18 720 runs
(the plan's 16 380 runs are the 7-network subset).

Pilot (--pilot): complete graph, CI-B, r3 = -0.3: 195 combinations = 585 runs.

Fase 3b (--phase3b): m_d in {0.1, 0.2} (m_d = 0 comes from 3a) on a reduced grid
around the crest: m_t in {0.05, 0.1, 0.2, 0.4}, 7 values of D_t in [0.03, 1],
D_d in {0.01, 0.03, 0.1}, 8 networks, CI-B, r3 = -0.3: 1344 combinations.

Output: Results/plan/phase3/phase3a_pilot.csv, phase3a_grid.csv, phase3b_grid.csv
"""

import os
import sys

from common import results_path
import plan_common as pc

MT_GRID = [0.0, 0.05, 0.1, 0.2, 0.4]
DT_GRID = pc.logspace_grid(13)
DD_GRID = [0.01, 0.03, 0.1]
R3_GRID = [-0.3, -0.45]


def specs_3a(pilot=False):
    if pilot:
        return [{"Network": "Complete", "CI": "B", "Dd": dd, "Dt": dt, "r3": -0.3, "tmax": 100,
                 "mt": mt} for mt in MT_GRID for dd in DD_GRID for dt in DT_GRID]
    return [{"Network": net, "CI": ci, "Dd": dd, "Dt": dt, "r3": r3, "tmax": 100, "mt": mt}
            for ci in ("B", "A") for net in pc.network_names(include_complete=True)
            for r3 in R3_GRID for mt in MT_GRID for dd in DD_GRID for dt in DT_GRID]


# Fase 3b grid, chosen after inspecting the 3a pilot: for m_t > 0 the crest of
# E_T lies in D_t in [0.03, 1], so the 7 values of the 3a grid in that range are
# used (stored in phase3b_dt_values.json). m_d = 0 is taken from 3a itself.
MD_GRID = [0.1, 0.2]
MT_3B = [0.05, 0.1, 0.2, 0.4]


def specs_3b(dt_values):
    return [{"Network": net, "CI": "B", "Dd": dd, "Dt": dt, "r3": -0.3, "tmax": 100,
             "mt": mt, "md": md}
            for net in pc.network_names(include_complete=True) for md in MD_GRID
            for mt in MT_3B for dd in DD_GRID for dt in dt_values]


if __name__ == "__main__":
    out = lambda n: results_path(os.path.join("plan", "phase3", n))
    if "--pilot" in sys.argv:
        pc.run_specs(specs_3a(True), out("phase3a_pilot.csv"))
    elif "--phase3b" in sys.argv:
        import json
        dt_values = [v for v in DT_GRID if v >= 0.03]
        json.dump(dt_values, open(out("phase3b_dt_values.json"), "w"))
        pc.run_specs(specs_3b(dt_values), out("phase3b_grid.csv"))
    else:
        pc.run_specs(specs_3a(), out("phase3a_grid.csv"))
