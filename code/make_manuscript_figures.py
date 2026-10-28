"""Publication figures (English labels, PDF) for the manuscript sections on
relative dispersal, the cost-free optimum and transit mortality.

Reads only the CSV outputs of phases 1-4 (Results/plan/phase*/) and writes
Results/plan/manuscript/fig4_*.pdf ... fig8_*.pdf and figS3_*.pdf ... figS6_*.pdf.
Run after analyze_phase1.py, analyze_phase2.py, analyze_phase3.py and
run_phase4_edge_mortality.py.
Usage: python make_manuscript_figures.py [fig4 figS3 ...]   (no arguments = all)
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

from common import results_path
import plan_common as pc
from plan_plots import DIV, SEQ, SERIES, NET_COLOR, NET_ORDER, R3_COLOR, DEGREE_MARKER

R = lambda *p: results_path(os.path.join("plan", *p))
OUT = R("manuscript")
os.makedirs(OUT, exist_ok=True)

LABEL = {"Path": "path, $k_s$=1", "Cycle": "cycle, $k_s$=2", "Random_s9": "random s9, $k_s$=3",
         "Random_s3": "random s3, $k_s$=4", "Random_s4": "random s4, $k_s$=5",
         "ScaleFree_s2": "scale-free s2, $k_s$=6", "Star": "star, $k_s$=7", "Complete": "complete, $k_s$=7"}
SHORT = {"Path": "path", "Cycle": "cycle", "Random_s9": "rand. s9", "Random_s3": "rand. s3",
         "Random_s4": "rand. s4", "ScaleFree_s2": "SF s2", "Star": "star", "Complete": "complete"}
MEAN_DEGREE = {n: float(np.mean([k for _, k in pc.network(n).degree()]))
               for n in pc.network_names(include_complete=True)}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "legend.fontsize": 6.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "lines.linewidth": 1.3, "savefig.dpi": 300, "pdf.fonttype": 42})


def tag(ax, s, x=-0.16, y=1.04):
    ax.text(x, y, f"({s})", transform=ax.transAxes, fontweight="bold", fontsize=8.5, va="bottom")


def read(*p):
    # plan CSVs hold no "None" strings, so default NA handling keeps numeric columns numeric
    return pd.read_csv(R(*p))


def plain_log_ticks(ax, axis, ticks):
    from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    a.set_minor_formatter(NullFormatter())


def grid_map(ax, g, xcol, ycol, val, cmap, norm=None, vmin=None, vmax=None):
    xs, ys = np.sort(g[xcol].unique()), np.sort(g[ycol].unique())
    Z = g.pivot_table(index=ycol, columns=xcol, values=val).loc[ys, xs].to_numpy()
    pc_ = ax.pcolormesh(xs, ys, Z, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, shading="nearest",
                        rasterized=True)
    ax.grid(False)
    return pc_, xs, ys, Z


# --------------------------------------------------------------------- figure 4
def fig4():
    d = read("phase1", "phase1_grid.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.6), constrained_layout=True)
    s = d[d.r3 == -0.3]
    for j, ci in enumerate(("A", "B")):
        ax = axes[0, j]
        g = s[s.CI == ci].groupby(["Dd", "Dt"])[["E_T", "E_ratio"]].median().reset_index()
        p, xs, ys, _ = grid_map(ax, g, "Dd", "Dt", "E_T", SEQ)
        _, _, _, Zr = grid_map(ax, g, "Dd", "Dt", "E_ratio", SEQ)
        ax.collections[-1].remove()
        ax.contour(xs, ys, Zr, levels=[1.0], colors="k", linewidths=1.2)
        ax.plot([1e-3, 1], [1e-3, 1], ls=":", color="0.35", lw=0.9)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("psyllid dispersal $D_d$"); ax.set_ylabel("parasitoid dispersal $D_t$")
        ax.set_title(f"IC-{ci}: targeted efficiency $E_T$")
        fig.colorbar(p, ax=ax, shrink=0.9)
        tag(ax, "ab"[j])
        ax = axes[1, j]
        sc = d[d.CI == ci]
        for r3, gg in sc.groupby("r3"):
            ax.scatter(gg.Pi_r1, gg.E_T, s=3, color=R3_COLOR[r3], alpha=0.45, linewidths=0,
                       rasterized=True, label=f"$r_3$ = {r3}")
        ax.set_xscale("log"); ax.set_xlabel(r"$\Pi = D_t\tau_t/(D_d\tau_d)$"); ax.set_ylabel("$E_T$")
        ax.set_title(f"IC-{ci}: collapse test, all networks")
        tag(ax, "cd"[j])
    axes[1, 0].legend(markerscale=3, loc="upper left")
    fig.savefig(os.path.join(OUT, "fig4_relative_dispersal.pdf"))
    plt.close(fig)


# --------------------------------------------------------------------- figure 5
def fig5():
    d = read("phase2", "phase2_grid.csv")
    c = read("phase2", "phase2_dd0_check.csv")
    ex2 = pd.read_csv(results_path("experiment2_outbreak_containment.csv"), keep_default_na=False)
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.4), constrained_layout=True)
    panels = [("A", 100, 0.03), ("B", 100, 0.03), ("B", 300, 0.03)]
    for ax, (ci, tm, dd), t in zip(axes, panels, "abc"):
        for n in pc.network_names():
            g = d[(d.Network == n) & (d.CI == ci) & (d.tmax == tm) & (d.Dd == dd)].sort_values("Dt")
            ax.plot(g.Dt, g.E_T, color=NET_COLOR[n], label=LABEL[n])
        ax.axvline(dd, color="0.5", ls=":", lw=0.9)
        ax.set_xscale("log"); ax.set_xlabel("$D_t$"); ax.set_title(f"IC-{ci}, $D_d$={dd}, $t_{{max}}$={tm}")
        tag(ax, t, x=-0.22)
    axes[0].set_ylabel("$E_T$")
    ax = axes[3]
    for n in pc.network_names():
        g = c[(c.Network == n) & (c.tmax == 100) & (c.Dt > 0)].sort_values("Dt")
        ax.plot(g.Dt, g.E_T, color=NET_COLOR[n])
    e0 = c[(c.tmax == 100) & (c.Dt == 0)].E_T.iloc[0]
    ax.axhline(e0, color="0.5", ls="--", lw=0.9)
    pts = []
    s = ex2[ex2.Dd == 0]
    for D in (0.01, 0.03, 0.1):
        v = [(lambda r: (r["None"] - r["TargetedAtSource"]) / 210)(
            s[(s.Topology == tp) & (s.Dt == D)].set_index("Strategy").PestAUC) for tp in s.Topology.unique()]
        pts.append(np.mean(v))
    ax.plot([0.01, 0.03, 0.1], pts, "ks", ms=4, label="Exp. 2")
    ax.set_xscale("log"); ax.set_xlabel("$D_t$"); ax.set_title("IC-A, $D_d$=0, $t_{max}$=100")
    ax.legend(loc="lower right")
    tag(ax, "d", x=-0.22)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=7, bbox_to_anchor=(0.5, -0.1), frameon=False)
    fig.savefig(os.path.join(OUT, "fig5_cost_free.pdf"), bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------- figure 6
def fig6():
    d = read("phase3", "phase3a_grid.csv")
    t = read("phase3", "table3_1_optimum.csv")
    t["real_max"] = t.real_max.astype(str) == "True"
    fig = plt.figure(figsize=(7.2, 5.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    s = d[(d.CI == "B") & (d.r3 == -0.3) & (d.Dd == 0.03) & (d.md == 0)]
    vmax = s.E_T.max()
    mts = sorted(s.mt.unique())
    for j, n in enumerate(("Complete", "Random_s4", "Star")):
        ax = fig.add_subplot(gs[0, j])
        g = s[s.Network == n]
        Dt = np.sort(g.Dt.unique())
        Z = g.pivot_table(index="mt", columns="Dt", values="E_T").loc[mts, Dt].to_numpy()
        p = ax.pcolormesh(Dt, np.arange(len(mts)), Z, cmap=SEQ, vmin=0, vmax=vmax, shading="nearest",
                          rasterized=True)
        tt = t[(t.Network == n) & (t.CI == "B") & (t.r3 == -0.3) & (t.Dd == 0.03) & (t.md == 0)].set_index("mt")
        for i, mt in enumerate(mts):
            o = tt.loc[mt]
            ax.plot(o.Dt_star, i, marker="o" if o.real_max else "x", color="#e34948",
                    mfc="white" if o.real_max else "#e34948", ms=4.5, mew=1.2)
        ax.set_xscale("log"); ax.grid(False)
        ax.set_yticks(range(len(mts)), [str(m) for m in mts])
        ax.set_xlabel("$D_t$"); ax.set_ylabel("$m_t$" if j == 0 else "")
        ax.set_title(LABEL[n])
        tag(ax, "abc"[j], x=-0.2)
    fig.colorbar(p, ax=fig.axes[:3], shrink=0.85, label="$E_T$")
    m = d[(d.Network == "Complete") & (d.CI == "B") & (d.r3 == -0.3) & (d.Dd == 0.01) & (d.md == 0)]
    axs = [fig.add_subplot(gs[1, j]) for j in range(3)]
    for k, mt in enumerate(mts):
        g = m[m.mt == mt].sort_values("Dt")
        axs[0].plot(g.Dt, g.E_T, color=SERIES[k], label=f"$m_t$ = {mt}")
        axs[1].plot(g.Dt, g.AUC_None / 1e3, color=SERIES[k], ls="--")
        axs[1].plot(g.Dt, g.AUC_Targeted / 1e3, color=SERIES[k])
        axs[2].plot(g.Dt, g.ParAUC_None / 1e3, color=SERIES[k])
    axs[0].set_ylabel("$E_T$"); axs[0].legend()
    axs[1].set_ylabel("pest AUC ($10^3$)")
    axs[1].legend(handles=[Line2D([], [], color="0.3", ls="--", label="no release"),
                           Line2D([], [], color="0.3", label="targeted")], loc="upper left")
    axs[2].set_ylabel("parasitoid AUC, no release ($10^3$)")
    for ax, l in zip(axs, "def"):
        ax.set_xscale("log"); ax.set_xlabel("$D_t$"); tag(ax, l, x=-0.2)
    axs[1].set_title("complete graph, IC-B, $D_d$ = 0.01")
    fig.savefig(os.path.join(OUT, "fig6_transit_mortality.pdf"))
    plt.close(fig)


# --------------------------------------------------------------------- figure 7
def fig7():
    t = read("phase3", "table3_1_optimum.csv")
    r = t[t.real_max & (t.mt > 0) & (t.md == 0)]
    s4 = read("phase4", "phase4_summary.csv")
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.5), constrained_layout=True)
    nets = [n for n in NET_ORDER if n in r.Network.unique()]
    for ax, col, lab, l in ((axes[0], "rule_ratio", r"$D_t^*m_t k_s/|r_3|$", "a"),
                            (axes[1], "rule_ratio_kmean", r"$D_t^*m_t\langle k\rangle/|r_3|$", "b")):
        for i, n in enumerate(nets):
            for r3, mk in ((-0.3, "o"), (-0.45, "^")):
                h = r[(r.Network == n) & (r.r3 == r3)]
                x = i + (0.18 if r3 == -0.45 else -0.18) + np.linspace(-0.08, 0.08, len(h))
                ax.scatter(x, h[col], s=9, marker=mk, color=NET_COLOR[n], linewidths=0)
        ax.axhline(1, color="k", lw=0.8)
        ax.set_yscale("log"); ax.set_ylim(0.3, 8)
        plain_log_ticks(ax, "y", [0.3, 0.5, 1, 2, 5])
        ax.set_xticks(range(len(nets)), [SHORT[n] for n in nets], rotation=60, ha="right")
        ax.set_ylabel(lab); tag(ax, l, x=-0.25)
    axes[0].legend(handles=[Line2D([], [], ls="", marker="o", color="0.3", ms=4, label="$r_3$ = -0.3"),
                            Line2D([], [], ls="", marker="^", color="0.3", ms=4, label="$r_3$ = -0.45")],
                   loc="upper left")
    lim = [0.1, 2]
    for ax, col, lab, l in ((axes[2], "pred_focus_mean", r"$c\,|r_3|/(\bar m_s k_s)$", "c"),
                            (axes[3], "pred_mean_degree", r"$c_k|r_3|/(\bar m\langle k\rangle)$", "d")):
        for k, ell in enumerate(sorted(s4.ell.unique())):
            h = s4[s4.ell == ell]
            ax.scatter(h[col], h.Dt_star, s=12, color=SERIES[k], linewidths=0, label=f"$\\ell$ = {ell}")
        ax.plot(lim, lim, "k-", lw=0.8)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
        plain_log_ticks(ax, "x", [0.1, 0.2, 0.5, 1, 2]); plain_log_ticks(ax, "y", [0.1, 0.2, 0.5, 1, 2])
        ax.set_xlabel("predicted $D_t^*$, " + lab); ax.set_ylabel("observed $D_t^*$")
        tag(ax, l, x=-0.25)
    axes[2].legend(loc="upper left")
    fig.savefig(os.path.join(OUT, "fig7_crest_rule.pdf"))
    plt.close(fig)


# --------------------------------------------------------------------- figure 8
def fig8():
    sc = read("phase3", "table3_4_star_vs_complete.csv")
    cr = read("phase3", "table3_3_degree_correlation.csv")
    t6 = read("phase3", "table3_6_phase3b_optimum.csv")
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.4), constrained_layout=True)
    mts = sorted(sc.mt.unique())
    ax = axes[0]
    for k, mt in enumerate(mts):
        g = sc[(sc.CI == "B") & (sc.Dd == 0.03) & (sc.mt == mt)].sort_values("Dt")
        ax.plot(g.Dt, g.rel_diff_pct_Targeted, color=SERIES[k], label=f"$m_t$ = {mt}")
    ax.set_xscale("log"); ax.set_xlabel("$D_t$"); ax.set_ylabel("AUC star $-$ complete (%)")
    ax.set_title("IC-B, targeted"); ax.legend(loc="lower left"); tag(ax, "a", x=-0.27)
    for ax, ci, l in ((axes[1], "A", "b"), (axes[2], "B", "c")):
        for k, mt in enumerate(mts):
            g = cr[(cr.CI == ci) & (cr.Dd == 0.03) & (cr.mt == mt)].sort_values("Dt")
            ax.plot(g.Dt, g.corr_k_AUC_T, color=SERIES[k])
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xscale("log"); ax.set_ylim(-1.05, 1.05); ax.set_xlabel("$D_t$")
        ax.set_ylabel("corr($k_s$, AUC targeted)"); ax.set_title(f"IC-{ci}, $D_d$ = 0.03")
        tag(ax, l, x=-0.27)
    ax = axes[3]
    for k, md in enumerate(sorted(t6.md.unique())):
        for dd, ls in ((0.03, "-"), (0.1, "--")):
            g = t6[(t6.Network == "Complete") & (t6.Dd == dd) & (t6.md == md) & (t6.mt > 0)].sort_values("mt")
            ax.plot(g.mt, g.Dt_star, color=SERIES[k], ls=ls, marker="o", ms=2.5,
                    label=f"$m_d$ = {md}" if dd == 0.03 else None)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("$m_t$"); ax.set_ylabel("$D_t^*$")
    plain_log_ticks(ax, "x", [0.05, 0.1, 0.2, 0.4]); plain_log_ticks(ax, "y", [0.03, 0.1, 0.3, 1])
    ax.set_title("complete graph, IC-B")
    h = ax.get_legend_handles_labels()[0] + [Line2D([], [], color="0.3", label="$D_d$ = 0.03"),
                                               Line2D([], [], color="0.3", ls="--", label="$D_d$ = 0.1")]
    ax.legend(handles=h, loc="upper right")
    tag(ax, "d", x=-0.27)
    fig.savefig(os.path.join(OUT, "fig8_structural_consequences.pdf"))
    plt.close(fig)


# --------------------------------------------------------------------- ESM figures
def figS3():
    d = read("phase1", "phase1_grid.csv")
    s = d[d.r3 == -0.3]
    nets = pc.network_names()
    fig, axes = plt.subplots(2, 7, figsize=(10, 3.4), constrained_layout=True, sharex=True, sharey=True)
    for i, ci in enumerate(("A", "B")):
        for j, n in enumerate(nets):
            p, *_ = grid_map(axes[i, j], s[(s.CI == ci) & (s.Network == n)], "Dd", "Dt", "E_T", SEQ,
                             vmin=s.E_T.min(), vmax=s.E_T.max())
            axes[i, j].set_xscale("log"); axes[i, j].set_yscale("log")
            axes[i, j].set_title(f"{LABEL[n]}, IC-{ci}", fontsize=6.5)
            if i == 1:
                axes[i, j].set_xlabel("$D_d$")
        axes[i, 0].set_ylabel("$D_t$")
    fig.colorbar(p, ax=axes, shrink=0.8, label="$E_T$")
    fig.savefig(os.path.join(OUT, "figS3_relative_dispersal_by_network.pdf"))
    plt.close(fig)


def figS4():
    d = read("phase2", "phase2_grid.csv")
    fig, axes = plt.subplots(4, 3, figsize=(7.2, 8.4), constrained_layout=True, sharex=True)
    for i, (ci, tm) in enumerate([("A", 100), ("A", 300), ("B", 100), ("B", 300)]):
        for j, dd in enumerate(sorted(d.Dd.unique())):
            ax = axes[i, j]
            for n in pc.network_names():
                g = d[(d.Network == n) & (d.CI == ci) & (d.tmax == tm) & (d.Dd == dd)].sort_values("Dt")
                ax.plot(g.Dt, g.E_T, color=NET_COLOR[n], label=LABEL[n])
            ax.axvline(dd, color="0.5", ls=":", lw=0.9)
            ax.set_xscale("log"); ax.set_title(f"IC-{ci}, $t_{{max}}$={tm}, $D_d$={dd}")
            if j == 0:
                ax.set_ylabel("$E_T$")
            if i == 3:
                ax.set_xlabel("$D_t$")
    axes[0, 0].legend()
    fig.savefig(os.path.join(OUT, "figS4_cost_free_all_curves.pdf"))
    plt.close(fig)


def figS5():
    d = read("phase3", "phase3a_grid.csv")
    t = read("phase3", "table3_1_optimum.csv")
    t["real_max"] = t.real_max.astype(str) == "True"
    nets = pc.network_names(include_complete=True)
    fig, axes = plt.subplots(4, 4, figsize=(7.2, 7.6), constrained_layout=True, sharex=True, sharey=True)
    for b, ci in enumerate(("B", "A")):
        s = d[(d.CI == ci) & (d.r3 == -0.3) & (d.Dd == 0.03) & (d.md == 0)]
        mts = sorted(s.mt.unique())
        for j, n in enumerate(nets):
            ax = axes[2 * b + j // 4, j % 4]
            g = s[s.Network == n]
            Dt = np.sort(g.Dt.unique())
            Z = g.pivot_table(index="mt", columns="Dt", values="E_T").loc[mts, Dt].to_numpy()
            p = ax.pcolormesh(Dt, np.arange(len(mts)), Z, cmap=SEQ, vmin=0, vmax=d.E_T.max(),
                              shading="nearest", rasterized=True)
            tt = t[(t.Network == n) & (t.CI == ci) & (t.r3 == -0.3) & (t.Dd == 0.03) & (t.md == 0)].set_index("mt")
            for i, mt in enumerate(mts):
                o = tt.loc[mt]
                ax.plot(o.Dt_star, i, marker="o" if o.real_max else "x", color="#e34948",
                        mfc="white" if o.real_max else "#e34948", ms=4, mew=1.1)
            ax.set_xscale("log"); ax.grid(False)
            ax.set_yticks(range(len(mts)), [str(m) for m in mts])
            ax.set_title(f"{LABEL[n]}, IC-{ci}", fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel("$m_t$")
    for ax in axes[-1]:
        ax.set_xlabel("$D_t$")
    fig.colorbar(p, ax=axes, shrink=0.6, label="$E_T$")
    fig.savefig(os.path.join(OUT, "figS5_transit_mortality_maps.pdf"))
    plt.close(fig)


def figS6():
    t6 = read("phase3", "table3_6_phase3b_optimum.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), constrained_layout=True, sharey=True)
    for ax, dd in zip(axes, sorted(t6.Dd.unique())):
        for n in NET_ORDER:
            for k, md in enumerate(sorted(t6.md.unique())):
                g = t6[(t6.Network == n) & (t6.Dd == dd) & (t6.md == md) & (t6.mt > 0)].sort_values("mt")
                if len(g):
                    ax.plot(g.mt, g.Dt_star, color=NET_COLOR[n], ls=["-", "--", ":"][k], lw=1,
                            label=LABEL[n] if k == 0 else None)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("$m_t$"); ax.set_title(f"$D_d$ = {dd}")
        plain_log_ticks(ax, "x", [0.05, 0.1, 0.2, 0.4]); plain_log_ticks(ax, "y", [0.03, 0.1, 0.3, 1])
    axes[0].set_ylabel("$D_t^*$ (restricted to $D_t \\geq 0.03$)")
    axes[0].legend(fontsize=5.5)
    fig.savefig(os.path.join(OUT, "figS6_pest_transit_mortality.pdf"))
    plt.close(fig)


FIGURES = {f.__name__: f for f in (fig4, fig5, fig6, fig7, fig8, figS3, figS4, figS5, figS6)}

if __name__ == "__main__":
    import sys
    for name in sys.argv[1:] or FIGURES:
        FIGURES[name]()
        print("done", name)
