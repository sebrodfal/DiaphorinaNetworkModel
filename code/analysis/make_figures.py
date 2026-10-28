"""Figures 2, 3 and S2, and the numbers of Experiments 2-7 quoted in the Results,
from the CSV files written by code/run_experiment*.py. Writes code/figures/.

Usage (from code/):
    python analysis/make_figures.py [path/to/Results]     (default: code/Results)

Only numpy, pandas, scipy, networkx and matplotlib are needed. The graph generators
below replicate network_experiments.make_topology so that the degree of the outbreak
node can be recovered for every graph instance used in Experiment 3.
"""

import os
import sys

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "Results")
OUT = os.path.join(HERE, "..", "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                     "font.family": "sans-serif", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.unicode_minus": False})
CM = 1 / 2.54
COL = {"None": "0.45", "Uniform": "#1f77b4", "TargetedAtSource": "#d95f02"}
LAB = {"None": "No release", "Uniform": "Uniform", "TargetedAtSource": "Targeted at source"}


def read(name):
    # keep_default_na=False: the strategy label "None" must not become NaN
    return pd.read_csv(os.path.join(RES, name), keep_default_na=False)


def make_topology(name, n, seed=1):
    if name == "Path":
        return nx.path_graph(n)
    if name == "Cycle":
        return nx.cycle_graph(n)
    if name == "Star":
        return nx.star_graph(n - 1)
    if name == "Complete":
        return nx.complete_graph(n)
    if name == "Random":
        return nx.gnp_random_graph(n, 0.4, seed=seed)
    if name == "SmallWorld":
        k = max(1, min(2, n // 2 - 1))
        return nx.watts_strogatz_graph(n, 2 * k, 0.2, seed=seed)
    if name == "ScaleFree":
        return nx.barabasi_albert_graph(n, min(2, n - 1), seed=seed)
    raise ValueError(name)


def panel(ax, letter):
    ax.text(-0.18, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=9)


def n_releases(T, tt=100):
    """Release events actually applied by core.release_times: kT < tt, k >= 1."""
    return int(np.ceil(tt / T) - 1)


# =============================================================== Experiment 2
e2 = read("experiment2_outbreak_containment.csv")
key = ["Topology", "Dd", "Dt"]
base = e2[e2.Strategy == "None"].set_index(key).PestAUC
print("== Experiment 2 ==")
print("runs:", len(e2))
print("PestPeak unique values:", e2.PestPeak.unique())
print("No-release AUC range: %.0f - %.0f" % (base.min(), base.max()))
eff = {}
for s, released in (("Uniform", 35 * 6 * 8), ("TargetedAtSource", 35 * 6 * 1)):
    x = e2[e2.Strategy == s].set_index(key).PestAUC
    red = 100 * (1 - x / base)
    per = ((base - x) / released).rename("per").reset_index()
    eff[s] = per.groupby(["Dd", "Dt"]).per.mean().unstack()
    print(f"{s}: reduction {red.min():.2f}-{red.max():.2f}% (mean {red.mean():.2f});"
          f" AUC avoided per unit released {per.per.min():.2f}-{per.per.max():.2f}")
    print(eff[s].round(2))
for s in ("None", "Uniform", "TargetedAtSource"):
    g = e2[e2.Strategy == s].groupby(["Dd", "Dt"]).PestAUC
    print(f"{s}: max spread across topologies within a (Dd,Dt) cell ="
          f" {(100 * (g.max() - g.min()) / g.mean()).max():.2f}%")
    print(e2[e2.Strategy == s].groupby(["Dd", "Dt"]).PestAUC.mean().unstack().round(0))
st = e2[(e2.Dd == 0.03) & (e2.Dt == 0.03)]
print(st[["Topology", "Strategy", "PestAUC", "PestSpatialCV", "PestSync"]].to_string())

# =============================================================== Experiment 3
e3 = read("experiment3_replicates_raw.csv")
rows = []
for topo in ["Path", "Cycle", "Star", "Complete", "Random", "SmallWorld", "ScaleFree"]:
    seeds = range(1, 11) if topo in ("Random", "SmallWorld", "ScaleFree") else [1]
    for seed in seeds:
        g = make_topology(topo, 8, seed)
        rows.append({"Topology": topo, "Seed": seed, "SourceDegree": g.degree(0),
                     "Connected": nx.is_connected(g),
                     "Lambda2": float(np.sort(np.linalg.eigvalsh(
                         nx.laplacian_matrix(g, nodelist=sorted(g)).toarray()))[1])})
G = pd.DataFrame(rows)
m = e3[(e3.Dd == 0.03) & (e3.Dt == 0.03)].merge(G, on=["Topology", "Seed"])
print("\n== Experiment 3 ==")
print("runs:", len(e3), " graphs:", len(G), " all connected:", G.Connected.all())
print("max |lambda2 recomputed - CSV|:", np.abs(m.Lambda2 - m.AlgebraicConnectivity).max())
print("min lambda2:", G.Lambda2.min())
grp = e3.groupby(["Topology", "Dd", "Dt", "Strategy"]).PestAUC
print("max CV of PestAUC across seeds (%):", 100 * (grp.std() / grp.mean()).max())
for s in ("None", "Uniform", "TargetedAtSource"):
    y = m[m.Strategy == s]
    print(s, "Spearman rho: AUC~source degree %.3f | AUC~lambda2 %.3f | AUC~mean degree %.3f"
          " | sync~source degree %.3f" % (spearmanr(y.SourceDegree, y.PestAUC)[0],
                                          spearmanr(y.Lambda2, y.PestAUC)[0],
                                          spearmanr(y.MeanDegree, y.PestAUC)[0],
                                          spearmanr(y.SourceDegree, y.PestSync)[0]))
    gg = y.groupby("SourceDegree").PestAUC.agg(["count", "mean", "min", "max"])
    print(gg.round(1))
    print("  change deg 1 -> 7: %.2f%%" % (100 * (gg["mean"].iloc[-1] / gg["mean"].iloc[0] - 1)))

# =============================================================== Figure 2
fig = plt.figure(figsize=(17.5 * CM, 6.2 * CM))
ax_a = fig.add_axes([0.07, 0.17, 0.20, 0.70])
ax_b = fig.add_axes([0.35, 0.17, 0.20, 0.70])
ax_cb = fig.add_axes([0.565, 0.17, 0.012, 0.70])
ax_c = fig.add_axes([0.73, 0.17, 0.26, 0.70])
vmin = min(eff[s].values.min() for s in eff)
vmax = max(eff[s].values.max() for s in eff)
for ax, s, letter in ((ax_a, "Uniform", "(a)"), (ax_b, "TargetedAtSource", "(b)")):
    M = eff[s]
    im = ax.imshow(M.values[::-1], cmap="viridis", vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(4), [f"{v:g}" for v in M.columns])
    ax.set_yticks(range(4), [f"{v:g}" for v in M.index[::-1]])
    ax.set_xlabel(r"parasitoid dispersal $D_t$")
    ax.set_ylabel(r"psyllid dispersal $D_d$")
    ax.set_title(LAB[s])
    for a in range(4):
        for b in range(4):
            v = M.values[::-1][a, b]
            ax.text(b, a, f"{v:.2f}", ha="center", va="center", fontsize=6,
                    color="white" if v < (vmin + vmax) / 2 else "black")
    panel(ax, letter)
cb = fig.colorbar(im, cax=ax_cb)
cb.set_label("AUC avoided per unit released")
ax = ax_c
for s in ("None", "Uniform", "TargetedAtSource"):
    y = m[m.Strategy == s]
    ref = y[y.SourceDegree == 1].PestAUC.mean()
    ax.plot(y.SourceDegree + {"None": -0.15, "Uniform": 0, "TargetedAtSource": 0.15}[s],
            100 * (y.PestAUC / ref - 1), "o", ms=3, color=COL[s], alpha=0.8, label=LAB[s])
    gm = y.groupby("SourceDegree").PestAUC.mean()
    ax.plot(gm.index, 100 * (gm / ref - 1), "-", color=COL[s], lw=1)
ax.set_xlabel("degree of the outbreak node")
ax.set_ylabel("change in network pest AUC (%)")
ax.set_xticks(range(1, 8))
ax.legend(frameon=False, loc="lower left")
panel(ax, "(c)")
fig.savefig(os.path.join(OUT, "fig2_strategy_topology.pdf"))
plt.close(fig)

# =============================================================== Experiment 6
e6 = read("experiment6_dose_frequency.csv")
e6["NRel"] = e6["T"].map(n_releases)
e6["NRelScript"] = np.floor(100 / e6["T"]).astype(int)
e6["Nodes"] = np.where(e6.Strategy == "Uniform", 8, 1)
e6["Released"] = e6.Delta * e6.NRel * e6.Nodes
e6["Efficiency"] = e6.AUCAvoided / e6.Released
e6["Reduction"] = 100 * e6.AUCAvoided / e6.BaselineAUC
print("\n== Experiment 6 ==")
print("runs:", len(e6), "+", e6.Topology.nunique(), "baselines")
print(e6[["T", "NRel", "NRelScript"]].drop_duplicates().to_string(index=False))
for (t, s), g in e6.groupby(["Topology", "Strategy"]):
    b0 = g.loc[g.CostEfficiency.idxmax()]
    b1 = g.loc[g.Efficiency.idxmax()]
    hi = g.loc[g.Reduction.idxmax()]
    print(f"{t}/{s}: script best (d={b0.Delta},T={b0['T']}) CE={b0.CostEfficiency:.2f};"
          f" corrected best (d={b1.Delta},T={b1['T']}) CE={b1.Efficiency:.2f} red={b1.Reduction:.1f}%;"
          f" max reduction {hi.Reduction:.1f}% at (d={hi.Delta},T={hi['T']}) CE={hi.Efficiency:.2f}")
print("Efficiency range by strategy:")
print(e6.groupby("Strategy").Efficiency.agg(["min", "max"]).round(2))
c = e6[e6.Topology == "Complete"]
for s in ("Uniform", "TargetedAtSource"):
    print(s, "reduction (%)")
    print(c[c.Strategy == s].pivot(index="Delta", columns="T", values="Reduction").round(1))
    print(s, "efficiency")
    print(c[c.Strategy == s].pivot(index="Delta", columns="T", values="Efficiency").round(2))
cu = c[c.Strategy == "Uniform"].copy()
cu["Rate"] = cu.Delta / cu["T"]
print("Uniform, Complete, schedules with delta/T = 2:")
print(cu[np.isclose(cu.Rate, 2)][["Delta", "T", "Reduction", "Released"]].to_string(index=False))
print("topology spread of reduction within a schedule (max, percentage points):",
      e6.groupby(["Strategy", "Delta", "T"]).Reduction.agg(lambda v: v.max() - v.min()).max())
print("min PestMin:", e6.PestMin.min())

# =============================================================== Figure 3
fig = plt.figure(figsize=(17.5 * CM, 6.2 * CM))
axs = [fig.add_axes([0.06, 0.17, 0.19, 0.70]), fig.add_axes([0.40, 0.17, 0.19, 0.70]),
       fig.add_axes([0.76, 0.17, 0.23, 0.70])]
caxs = [fig.add_axes([0.26, 0.17, 0.011, 0.70]), fig.add_axes([0.60, 0.17, 0.011, 0.70])]
for ax, cax, val, title, cmap in (
        (axs[0], caxs[0], "Reduction", "reduction in pest AUC (%)", "magma_r"),
        (axs[1], caxs[1], "Efficiency", "AUC avoided per unit released", "viridis")):
    P = cu.pivot(index="Delta", columns="T", values=val)
    im = ax.imshow(P.values[::-1], cmap=cmap, aspect="auto")
    ax.set_xticks(range(P.shape[1]), [str(v) for v in P.columns])
    ax.set_yticks(range(P.shape[0]), [str(v) for v in P.index[::-1]])
    ax.set_xlabel(r"release period $T$")
    ax.set_ylabel(r"release size $\delta$")
    for a in range(P.shape[0]):
        for b in range(P.shape[1]):
            v = P.values[::-1][a, b]
            dark = (v > 45) if val == "Reduction" else (v < 4.3)
            ax.text(b, a, f"{v:.1f}" if val == "Reduction" else f"{v:.2f}",
                    ha="center", va="center", fontsize=6, color="white" if dark else "black")
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(title)
panel(axs[0], "(a)")
panel(axs[1], "(b)")
ax = axs[2]
mk = {5: "o", 10: "s", 15: "^", 20: "D", 30: "v"}
for s in ("Uniform", "TargetedAtSource"):
    y = c[c.Strategy == s]
    for T, yy in y.groupby("T"):
        ax.plot(yy.Released, yy.Reduction, mk[T], ms=3.5, color=COL[s],
                mfc=COL[s] if s == "Uniform" else "none", ls="none",
                label=f"T = {T}" if s == "Uniform" else None)
ax.set_xscale("log")
ax.set_xlabel("total parasitoids released")
ax.set_ylabel("reduction in pest AUC (%)")
ax.legend(frameon=False, fontsize=6, loc="upper left")
ax.text(0.97, 0.42, "filled: Uniform\nopen: Targeted", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=6)
panel(ax, "(c)")
fig.savefig(os.path.join(OUT, "fig3_dose_frequency.pdf"))
plt.close(fig)

# =============================================================== Experiment 7 (figure S2)
e7 = read("experiment7_sensitivity.csv")
rank = read("experiment7_sensitivity_ranking.csv")
print("\n== Experiment 7 ==")
print(rank.round(2).to_string(index=False))
order = rank.Parameter.tolist()[::-1]
lo = [e7[(e7.Parameter == p) & np.isclose(e7.Factor, 0.9)].PctChangeAUC.iloc[0] for p in order]
hi = [e7[(e7.Parameter == p) & np.isclose(e7.Factor, 1.1)].PctChangeAUC.iloc[0] for p in order]
print("same-sign responses (non-monotone):",
      [p for p, a, b in zip(order, lo, hi) if np.sign(a) == np.sign(b)])
fig, ax = plt.subplots(figsize=(8.5 * CM, 8.0 * CM))
yv = np.arange(len(order))
ax.barh(yv + 0.2, lo, 0.4, color="#9ecae1", label=r"parameter $\times\,0.9$")
ax.barh(yv - 0.2, hi, 0.4, color="#08519c", label=r"parameter $\times\,1.1$")
TEXLAB = {p: "$" + p[0] + "_{" + p[1:] + "}$" for p in order}
ax.set_yticks(yv, [TEXLAB[p] for p in order])
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("change in network pest AUC (%)")
ax.legend(frameon=False, loc="lower right")
fig.savefig(os.path.join(OUT, "figS2_sensitivity.pdf"), bbox_inches="tight")
plt.close(fig)

# =============================================================== Experiment 5
e5 = read("experiment5_bifurcation_map.csv")
print("\n== Experiment 5 ==")
print("regimes:", e5.Regime.value_counts().to_dict())
print("dispersion range: %.2f - %.2f" % (e5.Dispersion.min(), e5.Dispersion.max()))
print(e5.nsmallest(5, "Dispersion").to_string(index=False))
print("mean dispersion by T:", e5.groupby("T").Dispersion.mean().round(1).to_dict())
e5.pivot(index="Delta", columns="T", values="Dispersion").to_csv(
    os.path.join(OUT, "_dispersion_table.csv"))

# =============================================================== Experiment 4
e4 = read("experiment4_stability_prediction.csv").set_index("Topology")
emp = e2[(e2.Dd == 0.03) & (e2.Dt == 0.03) & (e2.Strategy == "Uniform")].set_index("Topology")
e4["EmpSync_seed1"] = emp.PestSync
e4["EmpCV_seed1"] = emp.PestSpatialCV
e4["EfoldTime"] = 1 / e4.MaxTransverseLyapunov.abs()
e4["Factor_tt100"] = np.exp(100 * e4.MaxTransverseLyapunov)
print("\n== Experiment 4 ==")
print(e4.round(4).to_string())
print("Spearman lambda_max vs sync: %.3f (p=%.3f)" % tuple(spearmanr(e4.MaxTransverseLyapunov, e4.EmpSync_seed1)))
print("Spearman lambda_max vs CV:   %.3f (p=%.3f)" % tuple(spearmanr(e4.MaxTransverseLyapunov, e4.EmpCV_seed1)))
print("Spearman lambda2 vs sync:    %.3f" % spearmanr(e4.AlgebraicConnectivity, e4.EmpSync_seed1)[0])
