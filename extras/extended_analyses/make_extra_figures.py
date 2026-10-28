"""Figures 9 and S7 of the extended analyses (phase 5 and the extended phase-3 grid).

They are not part of the submitted manuscript. Reads the CSV files in Results/ next
to this script and writes figures/fig9_*.pdf and figures/figS7_*.pdf, with the style
of code/make_manuscript_figures.py. Run after analyze_extra.py.
Usage: python make_extra_figures.py [fig9 figS7]   (no arguments = both)
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from make_manuscript_figures import LABEL, SHORT, plain_log_ticks, tag
from plan_plots import NET_COLOR, NET_ORDER, SERIES

OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)


def read(name):
    return pd.read_csv(os.path.join(HERE, "Results", name))


# --------------------------------------------------------------------- figure 9 (phase 5)
def fig9():
    d = read("no_resident_grid.csv")
    t = read("table_B_no_resident_optimum.csv")
    h = d[(d.Dd == 0.03) & (d.r3 == -0.3)]
    mts = sorted(h.mt.unique())
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.4), constrained_layout=True)
    for ax, n, l in zip(axes[:3], ("Path", "Random_s3", "Complete"), "abc"):
        for k, mt in enumerate(mts):
            g = h[(h.Network == n) & (h.mt == mt)].sort_values("Dt")
            ax.plot(g.Dt, g.AUC_Targeted / g.AUC_None, color=SERIES[k], label=f"$m_t$ = {mt}")
            o = t[(t.Network == n) & (t.Dd == 0.03) & (t.r3 == -0.3) & (t.mt == mt)].iloc[0]
            if o.interior:
                ax.plot(o.Dt_star, o.AUC_T_min / o.AUC_None, "o", ms=3.5, color=SERIES[k],
                        mfc="white" if o.real_opt else SERIES[k])
        ax.set_xscale("log"); ax.set_ylim(0.4, 1.02)
        ax.set_xlabel("$D_t$"); ax.set_title(LABEL[n])
        ax.set_ylabel("pest AUC, targeted / no release" if l == "a" else "")
        tag(ax, l, x=-0.22)
    axes[0].legend(loc="upper left")
    ax = axes[3]
    s = t[(t.Dd == 0.03) & (t.r3 == -0.3) & (t.mt > 0) & t.interior]
    for n in NET_ORDER:
        g = s[s.Network == n].sort_values("mt")
        ax.plot(g.mt, g.Dt_star, color=NET_COLOR[n], marker="o", ms=2.5, lw=1, label=SHORT[n])
    ax.set_xscale("log"); ax.set_yscale("log")
    plain_log_ticks(ax, "x", [0.05, 0.1, 0.2, 0.4]); plain_log_ticks(ax, "y", [0.003, 0.01, 0.03, 0.1])
    ax.legend(fontsize=5, ncol=2, loc="upper right", handlelength=1.2)
    ax.set_xlabel("$m_t$"); ax.set_ylabel("control-optimal $D_t^*$"); ax.set_title("all networks")
    tag(ax, "d", x=-0.27)
    fig.savefig(os.path.join(OUT, "fig9_control_optimum_no_resident.pdf"), bbox_inches="tight")
    plt.close(fig)


def figS7():
    t = read("table_A_ext_optimum.csv")
    t["real_max"] = t.real_max.astype(str) == "True"
    pos = t[(t.mt > 0) & (t.md == 0)]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5), constrained_layout=True)
    for ax, r3, l in zip(axes[:2], (-0.3, -0.45), "ab"):
        h = pos[(pos.CI == "B") & (pos.Dd == 0.03) & (pos.r3 == r3)]
        for n in NET_ORDER:
            g = h[h.Network == n].sort_values("mt")
            ax.plot(g.mt, g.Dt_star, color=NET_COLOR[n], marker="o", ms=2.5, lw=1, label=SHORT[n])
        ax.set_xscale("log"); ax.set_yscale("log")
        plain_log_ticks(ax, "x", [0.05, 0.1, 0.2, 0.4]); plain_log_ticks(ax, "y", [0.1, 0.3, 1, 3])
        ax.set_xlabel("$m_t$"); ax.set_ylabel("$D_t^*$ (efficiency crest)")
        ax.set_title(f"IC-B, $D_d$ = 0.03, $r_3$ = {r3}")
        tag(ax, l, x=-0.25)
    axes[0].legend(fontsize=5.5, ncol=2, loc="lower left")
    ax = axes[2]
    nets = [n for n in NET_ORDER if n in pos.Network.unique()]
    for i, n in enumerate(nets):
        for r3, mk in ((-0.3, "o"), (-0.45, "^")):
            g = pos[pos.real_max & (pos.Network == n) & (pos.r3 == r3)]
            x = i + (0.18 if r3 == -0.45 else -0.18) + np.linspace(-0.08, 0.08, len(g))
            ax.scatter(x, g.rule_ratio_kmean, s=6, marker=mk, color=NET_COLOR[n], linewidths=0)
    ax.axhline(1, color="k", lw=0.8)
    ax.set_yscale("log"); ax.set_ylim(0.4, 3)
    plain_log_ticks(ax, "y", [0.5, 1, 2])
    ax.set_xticks(range(len(nets)), [SHORT[n] for n in nets], rotation=60, ha="right")
    ax.set_ylabel(r"$D_t^*m_t\langle k\rangle/|r_3|$")
    ax.legend(handles=[Line2D([], [], ls="", marker="o", color="0.3", ms=3, label="$r_3$ = -0.3"),
                       Line2D([], [], ls="", marker="^", color="0.3", ms=3, label="$r_3$ = -0.45")],
              loc="upper right")
    tag(ax, "c", x=-0.25)
    fig.savefig(os.path.join(OUT, "figS7_crest_rule_extended_grid.pdf"))
    plt.close(fig)


FIGURES = {f.__name__: f for f in (fig9, figS7)}

if __name__ == "__main__":
    for name in sys.argv[1:] or FIGURES:
        FIGURES[name]()
        print("done", name)
