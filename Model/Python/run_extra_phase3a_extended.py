"""Extra (after the 2026-10-17 review): Fase 3a with the D_t grid extended above 1.

The Fase-3a grid stops at D_t = 1, where 255 of the 384 curves with m_t > 0 have
their maximum, so most slopes of log D_t* vs log m_t used only two points. This
run adds the six next values of the same grid (step 0.25 decades): D_t = 10^0.25
... 10^1.5 (1.78 to 31.6), for every other parameter of Fase 3a. Combined with
phase3a_grid.csv (read only, not modified) it gives 19 values in [1e-3, 31.6].

6 x 5 x 3 x 8 x 2 x 2 = 2880 combinations = 8640 runs (None, Uniform, Targeted).
Output: Results/plan/extra/phase3a_ext_grid.csv (resumable, like the phase runs)
"""

import os

import numpy as np

from common import results_path
import plan_common as pc
from run_phase3_transit_mortality import DD_GRID, MT_GRID, R3_GRID

DT_EXT = list(10 ** (0.25 * np.arange(1, 7)))


def specs():
    return [{"Network": net, "CI": ci, "Dd": dd, "Dt": dt, "r3": r3, "tmax": 100, "mt": mt}
            for ci in ("B", "A") for net in pc.network_names(include_complete=True)
            for r3 in R3_GRID for mt in MT_GRID for dd in DD_GRID for dt in DT_EXT]


if __name__ == "__main__":
    pc.run_specs(specs(), results_path(os.path.join("plan", "extra", "phase3a_ext_grid.csv")))
