"""Fase 3 analysis (3a, 3b and the AUC part of 3c).

Input : Results/plan/phase3/phase3a_grid.csv (and phase3b_grid.csv if present)
Output: fig3_1_maps_Dt_mt.png   E_T over (D_t, m_t) per network, with the crest
        fig3_2_rule.png         D_t* m_t k_s / |r3| in every condition
        fig3_3_degree.png       correlation between focus degree and AUC vs m_t
        fig3_4_star_vs_complete.png   (3c) AUC difference star - complete vs m_t
        fig3_5_phase3b.png      (3b) shift of D_t* with m_d
        table3_1_optimum.csv    D_t*, Pi_t*, prominence, rule ratio
        gate_phase4.json        criterion to go on to phase 4
"""

import json
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from common import results_path
from analyze_phase2 import optimum
from plan_plots import NET_COLOR, NET_LABEL, NET_ORDER, SEQ, SERIES

OUT = results_path(os.path.join("plan", "phase3"))
import plan_common as pc
MEAN_DEGREE = {n: float(np.mean([k for _, k in pc.network(n).degree()]))
               for n in pc.network_names(include_complete=True)}
KEYS = ["Network", "CI", "Dd", "r3", "mt", "md"]


def table(d):
    rows = []
    for key, g in d.groupby(KEYS):
        o = optimum(g, "E_T")
        row = dict(zip(KEYS, key))
        k_s = int(g.k_s.iloc[0])
        row.update(k_s=k_s, **o)
        Dt, mt, r3 = o["Dt_star"], row["mt"], row["r3"]
        row["Pi_t_star"] = Dt * k_s / (abs(r3) + mt * Dt * k_s)
        row["Pi_t_ceiling"] = 1 / mt if mt > 0 else np.inf
        row["rule_ratio"] = Dt * mt * k_s / abs(r3) if mt > 0 else np.nan
        kbar = MEAN_DEGREE[row["Network"]]
        row["k_mean"] = kbar
        row["rule_ratio_kmean"] = Dt * mt * kbar / abs(r3) if mt > 0 else np.nan
        row["real_max"] = bool(o["interior"] and o["prominence"] > 0.05)
        # mechanism diagnostics (not in the plan)
        gs = g.sort_values("Dt").reset_index(drop=True)
        i = o["i_star"]
        row["ParAUC_None_at_star"] = gs.ParAUC_None[i]
        row["ParAUC_Targeted_at_star"] = gs.ParAUC_Targeted[i]
        row["AUC_None_at_star"] = gs.AUC_None[i]
        j = int(np.argmin(gs.AUC_Targeted))
        row["Dt_min_AUC_Targeted"] = gs.Dt[j]
        row["AUC_Targeted_min_interior"] = bool(0 < j < len(gs) - 1)
        row["AUC_Targeted_min_depth"] = (gs.AUC_Targeted.iloc[[0, -1]].min() - gs.AUC_Targeted[j]) \
            / gs.AUC_Targeted.iloc[[0, -1]].min()
        rows.append(row)
    return pd.DataFrame(rows)


def fig_maps(d, t, ci="B", r3=-0.3, dd=0.03):
    nets = [n for n in NET_ORDER if n in d.Network.unique()]
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), constrained_layout=True, sharex=True, sharey=True)
    s = d[(d.CI == ci) & (d.r3 == r3) & (d.Dd == dd) & (d.md == 0)]
    vmin, vmax = s.E_T.min(), s.E_T.max()
    mts = sorted(s.mt.unique())
    for ax, n in zip(axes.ravel(), nets):
        g = s[s.Network == n]
        Dt = np.sort(g.Dt.unique())
        Z = g.pivot(index="mt", columns="Dt", values="E_T").loc[mts, Dt].to_numpy()
        pc = ax.pcolormesh(Dt, np.arange(len(mts)), Z, cmap=SEQ, vmin=vmin, vmax=vmax, shading="nearest")
        tt = t[(t.Network == n) & (t.CI == ci) & (t.r3 == r3) & (t.Dd == dd) & (t.md == 0)].set_index("mt")
        for j, mt in enumerate(mts):
            o = tt.loc[mt]
            ax.plot(o.Dt_star, j, marker="o" if o.real_max else "x", color="#e34948",
                    mfc="white" if o.real_max else "#e34948", ms=7, mew=1.6)
        k = int(g.k_s.iloc[0])
        mm = np.array([m for m in mts if m > 0])
        ax.plot(abs(r3) / (mm * k), [mts.index(m) for m in mm], "k--", lw=1.2)
        ax.set_xscale("log"); ax.grid(False)
        ax.set_xlim(Dt.min() / 1.2, Dt.max() * 1.2)
        ax.set_yticks(range(len(mts)), [str(m) for m in mts])
        ax.set_title(NET_LABEL[n], fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("$m_t$")
    for ax in axes[-1]:
        ax.set_xlabel("$D_t$")
    fig.colorbar(pc, ax=axes, shrink=0.8, label="$E_T$")
    fig.suptitle(f"Figura 3.1. $E_T$ en el plano ($D_t$, $m_t$). CI-{ci}, $r_3$={r3}, $D_d$={dd}, $m_d$=0. "
                 "Círculo rojo: máximo interior con prominencia > 5 %; cruz: máximo en el borde o poco prominente; "
                 "línea discontinua: regla $D_t m_t k_s = |r_3|$.", fontsize=9.5)
    fig.savefig(os.path.join(OUT, "fig3_1_maps_Dt_mt.png"))
    plt.close(fig)


def fig_rule(t):
    s = t[(t.mt > 0) & (t.md == 0)]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True, sharey=True)
    for row, (col, lab) in enumerate((("rule_ratio", r"$D_t^*\,m_t\,k_s\,/\,|r_3|$  (grado del foco)"),
                                      ("rule_ratio_kmean", r"$D_t^*\,m_t\,\langle k\rangle\,/\,|r_3|$  (grado medio)"))):
        for ax, ci in zip(axes[row], ("B", "A")):
            g = s[s.CI == ci]
            mts = sorted(g.mt.unique())
            nets = [n for n in NET_ORDER if n in g.Network.unique()]
            for i, n in enumerate(nets):
                h = g[g.Network == n]
                x = np.array([mts.index(m) for m in h.mt]) + (i - len(nets) / 2) * 0.07 \
                    + h.r3.map({-0.3: -0.015, -0.45: 0.015}).to_numpy()
                real = h.real_max.to_numpy()
                ax.scatter(x[real], h[col][real], color=NET_COLOR[n], s=20, marker="o",
                           label=NET_LABEL[n])
                ax.scatter(x[~real], h[col][~real], facecolor="none", edgecolor=NET_COLOR[n], s=20,
                           marker="o", alpha=0.6)
            r = g[g.real_max][col]
            ax.axhline(1, color="k", lw=1)
            ax.set_yscale("log")
            ax.set_xticks(range(len(mts)), [str(m) for m in mts])
            ax.set_xlabel("$m_t$")
            ax.set_ylabel(lab)
            ax.set_title(f"CI-{ci}" + (" (caso principal)" if ci == "B" else " (control)")
                         + f": máximos reales, mediana {np.exp(np.log(r).median()):.2f}, "
                         f"rango {r.min():.2f}–{r.max():.2f}", fontsize=9)
    axes[0, 0].legend(fontsize=7.5, ncol=2)
    fig.suptitle("Figura 3.2. Prueba de la regla $D_t^* m_t k \\approx |r_3|$ en todas las condiciones "
                 "(3 valores de $D_d$, 2 de $r_3$). Relleno: máximo interior con prominencia > 5 %; "
                 "vacío: máximo en el borde ($D_t^*=1$) o poco prominente, donde el cociente es solo una cota inferior.",
                 fontsize=9.5)
    fig.savefig(os.path.join(OUT, "fig3_2_rule.png"))
    plt.close(fig)


def degree_correlations(d, ci, r3=-0.3):
    s = d[(d.CI == ci) & (d.r3 == r3) & (d.md == 0) & (d.Network != "Complete")]
    rows = []
    for (dd, mt, dt), g in s.groupby(["Dd", "mt", "Dt"]):
        rows.append({"Dd": dd, "mt": mt, "Dt": dt,
                     "corr_k_AUC_T": np.corrcoef(g.k_s, g.AUC_Targeted)[0, 1],
                     "corr_k_E_T": np.corrcoef(g.k_s, g.E_T)[0, 1]})
    return pd.DataFrame(rows)


def fig_degree(d):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True, sharex=True, sharey=True)
    out = []
    for i, ci in enumerate(("B", "A")):
        c = degree_correlations(d, ci)
        c["CI"] = ci
        out.append(c)
        for j, dd in enumerate(sorted(c.Dd.unique())):
            ax = axes[i, j]
            for k, mt in enumerate(sorted(c.mt.unique())):
                h = c[(c.Dd == dd) & (c.mt == mt)].sort_values("Dt")
                ax.plot(h.Dt, h.corr_k_AUC_T, color=SERIES[k], label=f"$m_t$={mt}")
            ax.axhline(0, color="k", lw=0.8)
            ax.set_xscale("log"); ax.set_ylim(-1.05, 1.05)
            ax.set_title(f"CI-{ci}, $D_d$={dd}")
            if i == 1:
                ax.set_xlabel("$D_t$")
            if j == 0:
                ax.set_ylabel("corr($k_s$, AUC liberación dirigida)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Figura 3.3. Correlación de Pearson entre el grado del foco y el AUC de plaga con liberación "
                 "dirigida, sobre las 7 redes ($r_3=-0.3$, $m_d=0$).", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig3_3_degree.png"))
    plt.close(fig)
    return pd.concat(out)


def fig_star_complete(d):
    s = d[d.Network.isin(["Star", "Complete"]) & (d.md == 0) & (d.r3 == -0.3)]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True, sharey=True)
    rows = []
    for ax, ci in zip(axes, ("B", "A")):
        for j, dd in enumerate(sorted(s.Dd.unique())):
            for k, mt in enumerate(sorted(s.mt.unique())):
                st = s[(s.Network == "Star") & (s.CI == ci) & (s.Dd == dd) & (s.mt == mt)].set_index("Dt").sort_index()
                co = s[(s.Network == "Complete") & (s.CI == ci) & (s.Dd == dd) & (s.mt == mt)].set_index("Dt").sort_index()
                rel = 100 * (st.AUC_Targeted - co.AUC_Targeted) / co.AUC_Targeted
                for dt, v in rel.items():
                    rows.append({"CI": ci, "Dd": dd, "mt": mt, "Dt": dt, "rel_diff_pct_Targeted": v,
                                 "rel_diff_pct_None": 100 * (st.AUC_None[dt] - co.AUC_None[dt]) / co.AUC_None[dt]})
                if dd == 0.03:
                    ax.plot(rel.index, rel.values, color=SERIES[k], label=f"$m_t$={mt}")
        ax.set_xscale("log"); ax.set_xlabel("$D_t$")
        ax.set_title(f"CI-{ci}, $D_d$=0.03, liberación dirigida")
    axes[0].set_ylabel("AUC estrella − AUC completo (% del completo)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Figura 3.4 (3c). Pérdida de la equivalencia estrella–completo con mortalidad en tránsito.",
                 fontsize=10)
    fig.savefig(os.path.join(OUT, "fig3_4_star_vs_complete.png"))
    plt.close(fig)
    return pd.DataFrame(rows)


def fig_mechanism(d, net="Complete", ci="B", r3=-0.3, dd=0.01):
    s = d[(d.Network == net) & (d.CI == ci) & (d.r3 == r3) & (d.Dd == dd) & (d.md == 0)]
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2), constrained_layout=True, sharex=True)
    panels = [("E_T", "$E_T$"), ("AUC_None", "AUC plaga sin liberación"),
              ("AUC_Targeted", "AUC plaga con liberación dirigida"),
              ("ParAUC_None", "AUC parasitoide sin liberación")]
    for k, mt in enumerate(sorted(s.mt.unique())):
        g = s[s.mt == mt].sort_values("Dt")
        for ax, (col, lab) in zip(axes, panels):
            ax.plot(g.Dt, g[col], color=SERIES[k], marker="o", ms=3, label=f"$m_t$={mt}")
    for ax, (col, lab) in zip(axes, panels):
        ax.set_xscale("log"); ax.set_xlabel("$D_t$"); ax.set_title(lab)
    axes[0].legend(fontsize=8)
    fig.suptitle(f"Figura 3.5. Mecanismo del máximo de $E_T$ ({NET_LABEL[net]}, CI-{ci}, $r_3$={r3}, "
                 f"$D_d$={dd}). El pico de $E_T$ coincide con el colapso del parasitoide residente en la "
                 "corrida sin liberación, no con un mínimo del AUC con liberación.", fontsize=9.5)
    fig.savefig(os.path.join(OUT, "fig3_5_mechanism.png"))
    plt.close(fig)


def fig_phase3b(t3a, b):
    tb = table(b)
    base = t3a[(t3a.CI == "B") & (t3a.r3 == -0.3) & (t3a.md == 0) & (t3a.mt > 0)]
    allt = pd.concat([base, tb[tb.mt > 0]])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True, sharey=True)
    for ax, dd in zip(axes, sorted(allt.Dd.unique())):
        for n in NET_ORDER:
            for k, md in enumerate(sorted(allt.md.unique())):
                h = allt[(allt.Network == n) & (allt.Dd == dd) & (allt.md == md)].sort_values("mt")
                if h.empty:
                    continue
                ax.plot(h.mt, h.Dt_star, color=NET_COLOR[n], ls=["-", "--", ":"][k], marker="o", ms=3,
                        label=f"{NET_LABEL[n]}" if k == 0 else None)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("$m_t$"); ax.set_title(f"$D_d$={dd}")
    axes[0].set_ylabel("$D_t^*$")
    axes[0].legend(fontsize=7, ncol=1)
    fig.suptitle("Figura 3.6 (3b). $D_t^*$ contra $m_t$; línea continua $m_d$=0, discontinua $m_d$=0.1, "
                 "punteada $m_d$=0.2 (CI-B, $r_3=-0.3$). Restringido a $D_t \\geq 0.03$.", fontsize=9.5)
    fig.savefig(os.path.join(OUT, "fig3_6_phase3b.png"))
    plt.close(fig)
    return allt


def gate(t):
    s = t[(t.CI == "B") & (t.md == 0) & (t.mt > 0)]
    per = s.groupby(["Network", "Dd", "r3"]).real_max.sum().rename("n_mt_with_real_max").reset_index()
    ok_combo = per[per.n_mt_with_real_max >= 2]
    # relation D_t*(m_t): slope of log D_t* vs log m_t over the real maxima, per network/Dd/r3
    slopes = []
    for key, g in s[s.real_max].groupby(["Network", "Dd", "r3"]):
        if len(g) >= 2:
            b = np.polyfit(np.log10(g.mt), np.log10(g.Dt_star), 1)[0]
            slopes.append(dict(zip(["Network", "Dd", "r3"], key), slope=b, n=len(g),
                               rule_ratio_median=float(g.rule_ratio.median())))
    slopes = pd.DataFrame(slopes)
    nets_ok = sorted(ok_combo.Network.unique())
    res = {"conditions_with_>=2_mt_real_max": int(len(ok_combo)),
           "of_conditions": int(len(per)),
           "networks_with_>=2_mt_real_max": nets_ok,
           "slopes_logDt*_vs_logmt": slopes.to_dict("records") if len(slopes) else [],
           }
    res["slope_median"] = float(slopes.slope.median()) if len(slopes) else None
    res["slope_iqr"] = [float(slopes.slope.quantile(.25)), float(slopes.slope.quantile(.75))] if len(slopes) else None
    res["pass"] = bool(len(nets_ok) >= 3 and len(slopes) >= 3
                       and (slopes.slope < 0).mean() >= 0.8)
    return res


if __name__ == "__main__":
    d = pd.read_csv(os.path.join(OUT, "phase3a_grid.csv"), keep_default_na=False)
    print(len(d), "combinations,", 3 * len(d), "runs")
    t = table(d)
    t.to_csv(os.path.join(OUT, "table3_1_optimum.csv"), index=False)
    fig_maps(d, t)
    fig_maps(d, t, ci="A")
    os.replace(os.path.join(OUT, "fig3_1_maps_Dt_mt.png"), os.path.join(OUT, "fig3_1b_maps_Dt_mt_CIA.png"))
    fig_maps(d, t)
    fig_rule(t)
    c = fig_degree(d)
    c.to_csv(os.path.join(OUT, "table3_3_degree_correlation.csv"), index=False)
    sc = fig_star_complete(d)
    sc.to_csv(os.path.join(OUT, "table3_4_star_vs_complete.csv"), index=False)
    fig_mechanism(d)
    pb = os.path.join(OUT, "phase3b_grid.csv")
    if os.path.exists(pb):
        b = pd.read_csv(pb, keep_default_na=False)
        b0 = d[(d.CI == "B") & (d.r3 == -0.3) & (d.Dt >= b.Dt.min() - 1e-12)]
        b0 = b0[b0.mt.isin(b.mt.unique())]
        allt = fig_phase3b(table(b0), b)
        allt.to_csv(os.path.join(OUT, "table3_6_phase3b_optimum.csv"), index=False)
        piv = allt.pivot_table(index=["Network", "Dd", "mt"], columns="md", values="Dt_star")
        print("\n3b: D_t* by m_d\n", piv.round(3).to_string())
        red = b.merge(b0[["Network", "Dd", "Dt", "mt", "AUC_None"]], on=["Network", "Dd", "Dt", "mt"],
                      suffixes=("", "_md0"))
        red["pest_reduction_by_md_pct"] = 100 * (1 - red.AUC_None / red.AUC_None_md0)
        red.groupby(["md"]).pest_reduction_by_md_pct.describe().to_csv(
            os.path.join(OUT, "table3_7_pest_reduction_by_md.csv"))
        print(red.groupby(["md"]).pest_reduction_by_md_pct.describe().round(2))
    g = gate(t)
    json.dump(g, open(os.path.join(OUT, "gate_phase4.json"), "w"), indent=1, default=float)
    pd.set_option("display.width", 220); pd.set_option("display.max_rows", 400)
    s = t[(t.CI == "B") & (t.r3 == -0.3) & (t.Dd == 0.03)]
    print(s[["Network", "mt", "k_s", "Dt_star", "interior", "prominence", "prom_left", "Pi_t_star",
             "rule_ratio", "real_max"]].round(4).to_string(index=False))
    print("gate:", {k: v for k, v in g.items() if k != "slopes_logDt*_vs_logmt"})
