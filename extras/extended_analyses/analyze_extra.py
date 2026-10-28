"""Analysis of the extra runs added after the 2026-10-17 review.

Inputs (read only): code/Results/plan/phase3/phase3a_grid.csv, code/Results/plan/phase1/phase1_grid.csv
                    Results/phase3a_ext_grid.csv, Results/no_resident_grid.csv (next to this script)
Outputs in Results/ next to this script:
  A. Fase 3a on the extended grid D_t in [1e-3, 31.6]
     table_A_ext_optimum.csv, table_A_ext_slopes.csv, table_A_ext_rule.csv, figA_ext_rule.png
  B. Releases without a resident parasitoid (control optimum of D_t)
     table_B_no_resident_optimum.csv, table_B_no_resident_slopes.csv, figB_no_resident.png
  C. Fase-1 collapse test with the fraction of pest AUC avoided as response
     table_C_collapse_relative.csv
  summary_extra.json
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import results_path
from plan_plots import NET_COLOR, NET_LABEL, NET_ORDER, SERIES
from analyze_phase2 import optimum
import analyze_phase3 as a3
import analyze_phase1 as a1

OUT = os.path.join(HERE, "Results")
P = lambda *p: results_path(os.path.join("plan", *p))


def slopes(t, keys):
    rows = []
    for key, g in t[t.real_max & (t.mt > 0)].groupby(keys):
        if len(g) >= 2:
            rows.append({**dict(zip(keys, key)), "n": len(g),
                         "slope": np.polyfit(np.log10(g.mt), np.log10(g.Dt_star), 1)[0]})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ A
def part_a(summary):
    base = pd.read_csv(P("phase3", "phase3a_grid.csv"))
    ext = pd.read_csv(os.path.join(OUT, "phase3a_ext_grid.csv"))
    d = pd.concat([base, ext], ignore_index=True)
    t = a3.table(d)
    t.to_csv(os.path.join(OUT, "table_A_ext_optimum.csv"), index=False)
    pos = t[(t.mt > 0) & (t.md == 0)]
    s = slopes(t, ["Network", "CI", "Dd", "r3"])
    s.to_csv(os.path.join(OUT, "table_A_ext_slopes.csv"), index=False)
    rule = []
    for (ci, r3), g in pos[pos.real_max].groupby(["CI", "r3"]):
        for col in ("rule_ratio", "rule_ratio_kmean"):
            rule.append({"CI": ci, "r3": r3, "version": col, "n": len(g), "median": g[col].median(),
                         "min": g[col].min(), "max": g[col].max(),
                         "sd_log10": np.log10(g[col]).std(),
                         "loss_at_crest_median": (g[col] * abs(r3)).median()})
    rule = pd.DataFrame(rule)
    rule.to_csv(os.path.join(OUT, "table_A_ext_rule.csv"), index=False)
    m0 = t[(t.mt == 0) & (t.md == 0)]
    summary["A_extended_grid"] = {
        "Dt_grid": [float(x) for x in np.sort(d.Dt.unique())],
        "curves_mt_pos": int(len(pos)),
        "real_maxima_mt_pos": int(pos.real_max.sum()),
        "at_right_edge_mt_pos": int((pos.i_star == d.Dt.nunique() - 1).sum()),
        "real_maxima_by_mt": {str(k): int(v) for k, v in pos[pos.real_max].groupby("mt").size().items()},
        "mt0_real_maxima": int(m0.real_max.sum()), "mt0_curves": int(len(m0)),
        "slopes_n_points": {str(k): int(v) for k, v in s.n.value_counts().items()},
        "slope_median": float(s.slope.median()),
        "slope_iqr": [float(s.slope.quantile(.25)), float(s.slope.quantile(.75))],
        "slope_median_CIB_r3_-0.3": float(s[(s.CI == "B") & (s.r3 == -0.3)].slope.median()),
    }
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for ax, r3 in zip(axes[:2], (-0.3, -0.45)):
        h = pos[(pos.CI == "B") & (pos.Dd == 0.03) & (pos.r3 == r3)]
        for n in NET_ORDER:
            g = h[h.Network == n].sort_values("mt")
            ax.plot(g.mt, g.Dt_star, color=NET_COLOR[n], marker="o", ms=4, label=NET_LABEL[n])
        mm = np.array([0.05, 0.4])
        ax.plot(mm, 0.1 / mm, "k--", lw=1, label="pendiente −1")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("$m_t$"); ax.set_ylabel("$D_t^*$")
        ax.set_title(f"CI-B, $D_d$=0.03, $r_3$={r3}")
    axes[0].legend(fontsize=7)
    ax = axes[2]
    for i, n in enumerate(NET_ORDER):
        g = pos[pos.real_max & (pos.Network == n)]
        ax.scatter(np.full(len(g), i) + np.random.default_rng(i).uniform(-0.2, 0.2, len(g)),
                   g.rule_ratio_kmean, s=10, color=NET_COLOR[n], alpha=0.7)
    ax.axhline(1, color="k", lw=0.8)
    ax.set_yscale("log"); ax.set_xticks(range(len(NET_ORDER)), NET_ORDER, rotation=45, ha="right")
    ax.set_ylabel(r"$D_t^*\,m_t\,\langle k\rangle/|r_3|$")
    ax.set_title("Todos los máximos reales (CI-A y CI-B, ambos $r_3$)")
    fig.suptitle("Fase 3a con malla extendida ($D_t$ hasta 31.6): posición de la cresta de $E_T$", fontsize=10)
    fig.savefig(os.path.join(OUT, "figA_ext_rule.png"))
    plt.close(fig)
    return t


# ------------------------------------------------------------------ B
def part_b(summary):
    d = pd.read_csv(os.path.join(OUT, "no_resident_grid.csv"))
    d["neg_AUC_T"] = -d.AUC_Targeted
    const = d.groupby(["Network", "Dd", "r3", "mt"]).AUC_None.agg(lambda x: (x.max() - x.min()) / x.mean())
    rows = []
    kbar = a3.MEAN_DEGREE
    for key, g in d.groupby(["Network", "Dd", "r3", "mt"]):
        g = g.sort_values("Dt")
        o = optimum(g, "neg_AUC_T")
        a = g.AUC_Targeted.to_numpy()
        amin = -o["E_star"]
        n, dd, r3, mt = key
        rows.append({"Network": n, "Dd": dd, "r3": r3, "mt": mt, "k_s": int(g.k_s.iloc[0]),
                     "k_mean": kbar[n], "Dt_star": o["Dt_star"], "interior": o["interior"],
                     "AUC_T_min": amin, "AUC_None": g.AUC_None.iloc[0],
                     "gain_vs_low_pct": 100 * (a[0] - amin) / a[0],
                     "gain_vs_high_pct": 100 * (a[-1] - amin) / a[-1],
                     "pest_avoided_at_opt_pct": 100 * (g.AUC_None.iloc[0] - amin) / g.AUC_None.iloc[0]})
    t = pd.DataFrame(rows)
    t["real_opt"] = t.interior & (t.gain_vs_low_pct > 5) & (t.gain_vs_high_pct > 5)
    t["loss_rate_at_opt"] = t.Dt_star * t.mt * t.k_mean
    t.to_csv(os.path.join(OUT, "table_B_no_resident_optimum.csv"), index=False)
    t2 = t.rename(columns={"real_opt": "real_max"})
    s = slopes(t2, ["Network", "Dd", "r3"])
    s.to_csv(os.path.join(OUT, "table_B_no_resident_slopes.csv"), index=False)
    summary["B_no_resident"] = {
        "AUC_None_max_rel_variation_over_Dt": float(const.max()),
        "curves": int(len(t)), "interior_minimum": int(t.interior.sum()),
        "real_optimum_(>5%_both_sides)": int(t.real_opt.sum()),
        "real_optimum_by_mt": {str(k): int(v) for k, v in t[t.real_opt].groupby("mt").size().items()},
        "real_optimum_by_network": {k: int(v) for k, v in t[t.real_opt].groupby("Network").size().items()},
        "median_gain_vs_low_pct_by_mt": {str(k): float(v) for k, v in
                                         t[t.interior].groupby("mt").gain_vs_low_pct.median().items()},
        "median_gain_vs_high_pct_by_mt": {str(k): float(v) for k, v in
                                          t[t.interior].groupby("mt").gain_vs_high_pct.median().items()},
        "slope_median": float(s.slope.median()) if len(s) else None,
        "slope_iqr": [float(s.slope.quantile(.25)), float(s.slope.quantile(.75))] if len(s) else None,
    }
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.5), constrained_layout=True, sharex=True, sharey=True)
    h = d[(d.Dd == 0.03) & (d.r3 == -0.3)]
    mts = sorted(h.mt.unique())
    for ax, n in zip(axes.ravel(), NET_ORDER):
        for k, mt in enumerate(mts):
            g = h[(h.Network == n) & (h.mt == mt)].sort_values("Dt")
            ax.plot(g.Dt, g.AUC_Targeted / g.AUC_None, color=SERIES[k], label=f"$m_t$={mt}")
            o = t[(t.Network == n) & (t.Dd == 0.03) & (t.r3 == -0.3) & (t.mt == mt)].iloc[0]
            if o.interior:
                ax.plot(o.Dt_star, o.AUC_T_min / o.AUC_None, "o", color=SERIES[k],
                        mfc="white" if o.real_opt else SERIES[k], ms=6)
        ax.set_xscale("log"); ax.set_title(NET_LABEL[n], fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("AUC plaga con liberación dirigida / sin liberación")
    for ax in axes[-1]:
        ax.set_xlabel("$D_t$")
    axes[0, 0].legend(fontsize=7.5)
    fig.suptitle("Liberación dirigida sin parasitoide residente ($x_3(0)=0$), plaga CI-B, $D_d$=0.03, "
                 "$r_3$=−0.3. Círculo: mínimo interior (blanco: >5 % mejor que ambos extremos).", fontsize=10)
    fig.savefig(os.path.join(OUT, "figB_no_resident.png"))
    plt.close(fig)
    return t


# ------------------------------------------------------------------ C
def part_c(summary):
    d = a1.load()
    d["Rel_T"] = (d.AUC_None - d.AUC_Targeted) / d.AUC_None
    tab = a1.table(d, "Rel_T")
    tab.to_csv(os.path.join(OUT, "table_C_collapse_relative.csv"), index=False)
    summary["C_collapse_relative"] = tab[["CI", "Candidate", "R2_cv_min", "Gain_vs_ratio"]].to_dict("records")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    summary = {}
    part_a(summary)
    part_b(summary)
    part_c(summary)
    json.dump(summary, open(os.path.join(OUT, "summary_extra.json"), "w"), indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))
