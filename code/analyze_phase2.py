"""Fase 2 analysis: is there an interior D_t that maximises E_T without
transit cost?

Input : Results/plan/phase2/phase2_grid.csv, phase2_dd0_check.csv
Output: fig2_1_ET_vs_Dt.png, fig2_2_dd0_check.png, table2_1_optimum.csv

For every curve E_T(D_t) (network, D_d, CI, t_max):
  D_t*       argmax on the grid, refined by a parabola in log10(D_t) through the
             maximum and its two neighbours (only when the maximum is interior)
  interior   the grid maximum is not at D_t = 1e-3 nor at D_t = 1
  prominence (E_T(D_t*) - E_T(D_t=1)) / E_T(D_t*)          [definition of the plan]
  prom_left  (E_T(D_t*) - E_T(D_t=1e-3)) / E_T(D_t*)       [reported in addition]
  Pi_t*      D_t* k_s / |r3|
A maximum counts as real when it is interior, its prominence exceeds 5 % and
both conditions hold for t_max = 100 and t_max = 300.
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import results_path
from plan_plots import NET_COLOR, NET_LABEL, NET_ORDER

OUT = results_path(os.path.join("plan", "phase2"))


def optimum(g, y="E_T"):
    g = g.sort_values("Dt")
    x, e = np.log10(g.Dt.to_numpy()), g[y].to_numpy()
    i = int(np.argmax(e))
    interior = 0 < i < len(e) - 1
    xs, es = x[i], e[i]
    if interior:
        a, b, c = np.polyfit(x[i - 1:i + 2], e[i - 1:i + 2], 2)
        if a < 0:
            xs = -b / (2 * a)
            es = c - b * b / (4 * a)
    return {"Dt_star": 10 ** xs, "E_star": es, "interior": interior, "i_star": i,
            "prominence": (es - e[-1]) / es, "prom_left": (es - e[0]) / es,
            "E_at_Dt1": e[-1], "E_at_Dtmin": e[0]}


def table(d, y="E_T", extra_keys=()):
    rows = []
    keys = ["Network", "CI", "Dd", "tmax", *extra_keys]
    for key, g in d.groupby(keys):
        o = optimum(g, y)
        k_s, r3 = g.k_s.iloc[0], g.r3.iloc[0]
        mt = g.mt.iloc[0] if "mt" in g else 0.0
        row = dict(zip(keys, key))
        row.update({"k_s": k_s, **o,
                    "Pi_t_star": o["Dt_star"] * k_s / (abs(r3) + mt * o["Dt_star"] * k_s)})
        rows.append(row)
    return pd.DataFrame(rows)


def real_maxima(t, keys=("Network", "CI", "Dd")):
    t = t.copy()
    t["ok"] = t.interior & (t.prominence > 0.05)
    both = t.groupby(list(keys)).ok.all().rename("real_in_both_horizons")
    return t.merge(both.reset_index(), on=list(keys))


def fig_curves(d, t):
    fig, axes = plt.subplots(4, 3, figsize=(13, 13), sharex=True, constrained_layout=True)
    for i, (ci, tm) in enumerate([("A", 100), ("A", 300), ("B", 100), ("B", 300)]):
        for j, dd in enumerate(sorted(d.Dd.unique())):
            ax = axes[i, j]
            for n in NET_ORDER:
                g = d[(d.Network == n) & (d.CI == ci) & (d.tmax == tm) & (d.Dd == dd)].sort_values("Dt")
                if g.empty:
                    continue
                ax.plot(g.Dt, g.E_T, color=NET_COLOR[n], label=NET_LABEL[n])
                o = t[(t.Network == n) & (t.CI == ci) & (t.tmax == tm) & (t.Dd == dd)].iloc[0]
                ax.plot(o.Dt_star, o.E_star, marker="o" if o.interior else "x", ms=6,
                        color=NET_COLOR[n], mfc="white" if o.interior else NET_COLOR[n])
            ax.axvline(dd, color="0.5", ls=":", lw=1)
            ax.set_xscale("log")
            ax.set_title(f"CI-{ci}, $t_{{max}}$={tm}, $D_d$={dd}")
            if j == 0:
                ax.set_ylabel("$E_T$")
            if i == 3:
                ax.set_xlabel("$D_t$")
    axes[0, 0].legend(fontsize=7.5, loc="best")
    fig.suptitle("Figura 2.1. $E_T$ contra $D_t$ sin costo de dispersión ($r_3=-0.3$). "
                 "Círculo: máximo interior; cruz: máximo en el borde de la malla; "
                 "línea punteada: $D_t=D_d$.", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig2_1_ET_vs_Dt.png"))
    plt.close(fig)


def fig_dd0(c):
    ex2 = pd.read_csv(results_path("experiment2_outbreak_containment.csv"), keep_default_na=False)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True, sharey=True)
    for ax, tm in zip(axes, (100, 300)):
        for n in NET_ORDER:
            g = c[(c.Network == n) & (c.tmax == tm)].sort_values("Dt")
            if g.empty:
                continue
            g0 = g[g.Dt == 0]
            g = g[g.Dt > 0]
            ax.plot(g.Dt, g.E_T, color=NET_COLOR[n], label=NET_LABEL[n])
            if not g0.empty:
                ax.axhline(g0.E_T.iloc[0], color="0.6", lw=0.8, ls="--")
        if tm == 100:
            s = ex2[ex2.Dd == 0]
            pts = []
            for D in (0.01, 0.03, 0.1):
                vals = []
                for topo in s.Topology.unique():
                    r = s[(s.Topology == topo) & (s.Dt == D)].set_index("Strategy").PestAUC
                    vals.append((r["None"] - r["TargetedAtSource"]) / 210)
                pts.append(np.mean(vals))
            ax.plot([0.01, 0.03, 0.1], pts, "ks", ms=7, label="Exp. 2, media de 7 topologías")
        ax.set_xscale("log"); ax.set_xlabel("$D_t$"); ax.set_title(f"$D_d=0$, CI-A, $t_{{max}}$={tm}")
    axes[0].set_ylabel("$E_T$")
    axes[0].legend(fontsize=7.5)
    fig.suptitle("Figura 2.2. Verificación con $D_d=0$. Línea gris discontinua: $E_T$ con $D_t=0$ "
                 "(5.30 en $t_{max}=100$).", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig2_2_dd0_check.png"))
    plt.close(fig)


if __name__ == "__main__":
    d = pd.read_csv(os.path.join(OUT, "phase2_grid.csv"), keep_default_na=False)
    print(len(d), "combinations,", 3 * len(d), "runs")
    t = real_maxima(table(d))
    t.to_csv(os.path.join(OUT, "table2_1_optimum.csv"), index=False)
    pd.set_option("display.width", 220); pd.set_option("display.max_rows", 200)
    cols = ["Network", "CI", "Dd", "tmax", "k_s", "Dt_star", "Pi_t_star", "interior",
            "prominence", "prom_left", "real_in_both_horizons"]
    print(t[cols].round(4).to_string(index=False))
    tu = real_maxima(table(d, "E_U"))
    tu.to_csv(os.path.join(OUT, "table2_1b_optimum_uniform.csv"), index=False)
    print("\nUniform release: real interior maxima:", int(tu.real_in_both_horizons.sum()), "of", len(tu))
    fig_curves(d, t)
    p = os.path.join(OUT, "phase2_dd0_check.csv")
    if os.path.exists(p):
        c = pd.read_csv(p, keep_default_na=False)
        fig_dd0(c)
        tc = table(c[c.Dt > 0])
        tc.to_csv(os.path.join(OUT, "table2_2_dd0_check.csv"), index=False)
        print("\nD_d = 0 check:\n", tc[["Network", "tmax", "Dt_star", "interior", "prominence",
                                        "prom_left", "E_at_Dtmin", "E_at_Dt1"]].round(3).to_string(index=False))
