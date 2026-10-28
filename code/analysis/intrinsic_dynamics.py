"""Figure S1 and the intrinsic-dynamics numbers quoted in the Results.

Needs the model code and figures/_dispersion_table.csv, which make_figures.py
writes. Run from code/:
    python analysis/intrinsic_dynamics.py [path/to/code]     (default: code/)

(a) single-node trajectory without releases (self-sustained oscillation) and with
    the reference schedule delta=35, T=15;
(b) log10 dispersion of the last 30 post-release states (Experiment 5 CSV);
plus the mu = 0 (synchronous-mode) Lyapunov exponent with the same settings as
Experiment 4.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODEL = os.path.abspath(sys.argv[1] if len(sys.argv) > 1
                        else os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, MODEL)
from common import PARAMS                                   # noqa: E402
from core import single_node_rhs, solve_impulsive, release_times  # noqa: E402
from network_stability import master_lyapunov_exponent      # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "legend.fontsize": 7, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})
CM = 1 / 2.54
f = single_node_rhs(PARAMS)

# Unforced system: late-time oscillation
free = solve_impulsive(f, [50, 10, 15, 30], (0.0, 3000.0))
ts = np.linspace(2000, 3000, 40001)
Y = free(ts)
z = Y[0] - Y[0].mean()
up = ts[1:][(z[:-1] < 0) & (z[1:] >= 0)]
print("unforced: x1 in [%.2f, %.2f], x3 in [%.2f, %.2f]" % (Y[0].min(), Y[0].max(), Y[2].min(), Y[2].max()))
print("unforced: mean period from upward crossings = %.3f (sd %.3f, n=%d)"
      % (np.diff(up).mean(), np.diff(up).std(), len(up) - 1))


def jump(y):
    y = y.copy()
    y[2] += 35
    return y


forced = solve_impulsive(f, [50, 10, 15, 30], (0.0, 300.0), release_times(15, 300), jump)
lam0 = master_lyapunov_exponent(PARAMS, 0.03, 0.03, 0.0, 35, 15)
print("synchronous mode (mu=0) exponent:", lam0)

fig, axs = plt.subplots(1, 2, figsize=(17.5 * CM, 6.0 * CM), gridspec_kw={"wspace": 0.35,
                                                                           "width_ratios": [1.4, 1]})
t = np.linspace(0, 150, 6001)
axs[0].plot(t, free(t)[0], color="0.4", lw=0.9, label="no release")
axs[0].plot(t, forced(t)[0], color="#1f77b4", lw=0.9, label=r"$\delta=35$, $T=15$")
axs[0].set_xlabel("time (a.u.)")
axs[0].set_ylabel(r"$D$ (psyllid density, a.u.)")
axs[0].legend(frameon=False, loc="upper right")
axs[0].text(-0.12, 1.04, "(a)", transform=axs[0].transAxes, fontweight="bold", fontsize=9)

D = pd.read_csv(os.path.join(OUT, "_dispersion_table.csv"), index_col=0)
im = axs[1].imshow(np.log10(D.values[::-1]), cmap="magma", aspect="auto")
axs[1].set_xticks(range(D.shape[1]), D.columns)
axs[1].set_yticks(range(D.shape[0]), D.index[::-1])
axs[1].set_xlabel(r"release period $T$")
axs[1].set_ylabel(r"release size $\delta$")
cb = fig.colorbar(im, ax=axs[1], fraction=0.046, pad=0.04)
cb.set_label(r"$\log_{10}$ dispersion")
axs[1].text(-0.2, 1.04, "(b)", transform=axs[1].transAxes, fontweight="bold", fontsize=9)
fig.savefig(os.path.join(OUT, "figS1_intrinsic_dynamics.pdf"), bbox_inches="tight")
