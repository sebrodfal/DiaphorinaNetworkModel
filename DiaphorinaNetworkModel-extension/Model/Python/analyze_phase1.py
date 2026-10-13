"""Fase 1 analysis: maps of E_T and E_T/E_U, collapse curves and the
cross-validated R^2 of each candidate variable.

Input : Results/plan/phase1/phase1_grid.csv
Output: Results/plan/phase1/fig1_1_maps.png, fig1_1b_maps_by_network.png,
        fig1_2_collapse.png, table1_1_cv_r2.csv (+ table1_1_cv_r2_Eratio.csv)

Collapse test. For each candidate a smooth curve (cubic polynomial in log10 of
the candidate; for the pair (Pi, Pi_t) a bivariate cubic) is fitted and its R^2
is computed on held-out data, leaving out one value of r3 at a time (3 folds)
or one whole network at a time (7 folds). Held-out values of the candidate
outside the training range are clamped to it, so the fitted curve is extended
as a constant instead of extrapolating a polynomial.
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.lines import Line2D

from common import results_path
from plan_plots import DEGREE_MARKER, DIV, NET_LABEL, NET_ORDER, R3_COLOR, SEQ

OUT = results_path(os.path.join("plan", "phase1"))

CANDIDATES = {
    "D_t": ["Dt"],
    "D_t/D_d": ["ratio"],
    "Pi (tau_d=1/|r1|)": ["Pi_r1"],
    "Pi (tau_d=T0)": ["Pi_T0"],
    "(Pi, Pi_t)": ["Pi_r1", "Pi_t"],
    "ref: (D_d, D_t)": ["Dd", "Dt"],          # not in the plan: upper reference
}


def poly_features(X, deg=3):
    X = np.atleast_2d(X)
    cols = [np.ones(len(X))]
    if X.shape[1] == 1:
        cols += [X[:, 0] ** p for p in range(1, deg + 1)]
    else:
        a, b = X[:, 0], X[:, 1]
        cols += [a ** i * b ** j for i in range(deg + 1) for j in range(deg + 1) if 0 < i + j <= deg]
    return np.column_stack(cols)


def fit_predict(Xtr, ytr, Xte, deg=3):
    lo, hi = Xtr.min(axis=0), Xtr.max(axis=0)
    Xte = np.clip(Xte, lo, hi)
    coef, *_ = np.linalg.lstsq(poly_features(Xtr, deg), ytr, rcond=None)
    return poly_features(Xte, deg) @ coef


def r2(y, yhat):
    return 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)


def cv_r2(df, cols, y, group):
    X = np.log10(df[cols].to_numpy(float))
    yy = df[y].to_numpy(float)
    pred = np.empty_like(yy)
    for g in df[group].unique():
        te = (df[group] == g).to_numpy()
        pred[te] = fit_predict(X[~te], yy[~te], X[te])
    return r2(yy, pred), pred


def load():
    d = pd.read_csv(os.path.join(OUT, "phase1_grid.csv"), keep_default_na=False)
    d["ratio"] = d.Dt / d.Dd
    return d


def table(d, y="E_T"):
    rows = []
    for ci in ("A", "B"):
        s = d[d.CI == ci].reset_index(drop=True)
        X_all = {}
        for name, cols in CANDIDATES.items():
            X = np.log10(s[cols].to_numpy(float))
            insample = r2(s[y].to_numpy(), fit_predict(X, s[y].to_numpy(), X))
            r_r3, _ = cv_r2(s, cols, y, "r3")
            r_net, _ = cv_r2(s, cols, y, "Network")
            rows.append({"CI": ci, "Candidate": name, "R2_insample": insample,
                         "R2_cv_leave_r3": r_r3, "R2_cv_leave_network": r_net,
                         "R2_cv_min": min(r_r3, r_net)})
    t = pd.DataFrame(rows)
    for ci in ("A", "B"):
        base = t[(t.CI == ci) & (t.Candidate == "D_t/D_d")].R2_cv_min.iloc[0]
        m = t.CI == ci
        t.loc[m, "Gain_vs_ratio"] = t.loc[m, "R2_cv_min"] - base
        t.loc[m, "Collapses"] = (t.loc[m, "R2_cv_min"] >= 0.9) & (t.loc[m, "Gain_vs_ratio"] >= 0.1) \
            & ~t.loc[m, "Candidate"].str.startswith("ref")
    return t


def fig_maps(d):
    s = d[d.r3 == -0.3]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.4), constrained_layout=True)
    for i, ci in enumerate(("A", "B")):
        g = s[s.CI == ci].groupby(["Dd", "Dt"])[["E_T", "E_ratio"]].median().reset_index()
        Dd, Dt = np.sort(g.Dd.unique()), np.sort(g.Dt.unique())
        for j, (col, cmap, lab) in enumerate((("E_T", SEQ, "$E_T$ (AUC evitado / parasitoide)"),
                                              ("E_ratio", DIV, "$E_T/E_U$"))):
            Z = g.pivot(index="Dt", columns="Dd", values=col).loc[Dt, Dd].to_numpy()
            ax = axes[i, j]
            norm = TwoSlopeNorm(1.0, min(Z.min(), 0.99), max(Z.max(), 1.01)) if col == "E_ratio" else None
            pc = ax.pcolormesh(Dd, Dt, Z, cmap=cmap, norm=norm, shading="nearest")
            if col == "E_ratio" and Z.min() < 1 < Z.max():
                ax.contour(Dd, Dt, Z, levels=[1.0], colors="k", linewidths=1.2)
            ax.plot([1e-3, 1], [1e-3, 1], ls=":", color="0.3", lw=1)
            ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(False)
            ax.set_xlabel("$D_d$ (psílido)"); ax.set_ylabel("$D_t$ (parasitoide)")
            ax.set_title(f"CI-{ci}: {lab}")
            fig.colorbar(pc, ax=ax, shrink=0.9)
    fig.suptitle("Figura 1.1. Mediana sobre las 7 redes, $r_3=-0.3$, $m_d=m_t=0$. "
                 "Línea punteada: $D_t=D_d$; línea negra: $E_T=E_U$.", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig1_1_maps.png"))
    plt.close(fig)

    s = d[d.r3 == -0.3]
    nets = [n for n in NET_ORDER if n in s.Network.unique()]
    fig, axes = plt.subplots(2, len(nets), figsize=(2.3 * len(nets), 5.2), constrained_layout=True,
                             sharex=True, sharey=True)
    vmin, vmax = s.E_T.min(), s.E_T.max()
    for i, ci in enumerate(("A", "B")):
        for j, n in enumerate(nets):
            g = s[(s.CI == ci) & (s.Network == n)]
            Dd, Dt = np.sort(g.Dd.unique()), np.sort(g.Dt.unique())
            Z = g.pivot(index="Dt", columns="Dd", values="E_T").loc[Dt, Dd].to_numpy()
            ax = axes[i, j]
            pc = ax.pcolormesh(Dd, Dt, Z, cmap=SEQ, vmin=vmin, vmax=vmax, shading="nearest")
            ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(False)
            ax.set_title(f"{NET_LABEL[n]}\nCI-{ci}", fontsize=8)
            if i == 1:
                ax.set_xlabel("$D_d$")
            if j == 0:
                ax.set_ylabel("$D_t$")
    fig.colorbar(pc, ax=axes, shrink=0.8, label="$E_T$")
    fig.suptitle("Figura 1.1b. $E_T$ por red ($r_3=-0.3$)", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig1_1b_maps_by_network.png"))
    plt.close(fig)


def fig_collapse(d):
    one_d = [("D_t", "Dt", "$D_t$"), ("D_t/D_d", "ratio", "$D_t/D_d$"),
             ("Pi (tau_d=1/|r1|)", "Pi_r1", r"$\Pi$ ($\tau_d=1/|r_1|$)"),
             ("Pi (tau_d=T0)", "Pi_T0", r"$\Pi$ ($\tau_d=T_0$)")]
    fig, axes = plt.subplots(2, 5, figsize=(17, 7), constrained_layout=True)
    for i, ci in enumerate(("A", "B")):
        s = d[d.CI == ci].reset_index(drop=True)
        for j, (name, col, lab) in enumerate(one_d):
            ax = axes[i, j]
            for (r3, k), g in s.groupby(["r3", "k_s"]):
                ax.scatter(g[col], g.E_T, s=7, color=R3_COLOR[r3], marker=DEGREE_MARKER[k],
                           alpha=0.45, linewidths=0)
            X = np.log10(s[[col]].to_numpy(float))
            xs = np.linspace(X.min(), X.max(), 200)[:, None]
            ax.plot(10 ** xs[:, 0], fit_predict(X, s.E_T.to_numpy(), xs), color="k", lw=1.6)
            rr, _ = cv_r2(s, [col], "E_T", "r3")
            rn, _ = cv_r2(s, [col], "E_T", "Network")
            ax.set_xscale("log"); ax.set_xlabel(lab); ax.set_ylabel("$E_T$")
            ax.set_title(f"CI-{ci}   $R^2_{{cv}}$ r3: {rr:.2f} · red: {rn:.2f}", fontsize=9)
        ax = axes[i, 4]
        rr, pred = cv_r2(s, ["Pi_r1", "Pi_t"], "E_T", "Network")
        for (r3, k), g in s.groupby(["r3", "k_s"]):
            ax.scatter(pred[g.index], g.E_T, s=7, color=R3_COLOR[r3], marker=DEGREE_MARKER[k],
                       alpha=0.45, linewidths=0)
        lim = [s.E_T.min(), s.E_T.max()]
        ax.plot(lim, lim, color="k", lw=1)
        ax.set_xlabel(r"$E_T$ predicha por $(\Pi,\Pi_t)$ (validación: red fuera)")
        ax.set_ylabel("$E_T$ observada")
        rr3, _ = cv_r2(s, ["Pi_r1", "Pi_t"], "E_T", "r3")
        ax.set_title(f"CI-{ci}   $R^2_{{cv}}$ r3: {rr3:.2f} · red: {rr:.2f}", fontsize=9)
    handles = [Line2D([], [], ls="", marker="o", color=c, label=f"$r_3$={r}") for r, c in R3_COLOR.items()]
    handles += [Line2D([], [], ls="", marker=m, color="0.35", label=f"$k_s$={k}") for k, m in DEGREE_MARKER.items()]
    fig.legend(handles=handles, loc="lower center", ncol=10, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Figura 1.2. Prueba de colapso de $E_T$ ($m_d=m_t=0$). Curva negra: ajuste cúbico en escala log "
                 "con todos los datos.", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig1_2_collapse.png"))
    plt.close(fig)


if __name__ == "__main__":
    d = load()
    print(len(d), "combinations,", 3 * len(d), "runs")
    t = table(d, "E_T")
    t.to_csv(os.path.join(OUT, "table1_1_cv_r2.csv"), index=False)
    pd.set_option("display.width", 200)
    print(t.round(3).to_string(index=False))
    t2 = table(d, "E_ratio")
    t2.to_csv(os.path.join(OUT, "table1_1_cv_r2_Eratio.csv"), index=False)
    print("\nResponse E_T/E_U:\n", t2.round(3).to_string(index=False))
    fig_maps(d)
    fig_collapse(d)
    print("figures saved in", OUT)
