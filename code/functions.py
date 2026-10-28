"""Single-node model, network version for plotting, and the single-node paper's
analysis functions (port of Functions.m).

params = [[r1,b12,b13,a1,c1],[r2,b21,b24,a2,c2],[r3,b31,a3,c3],[r4,b42,a4,c4]]
init   = [Di0, B0, x30, x40, tt]

Solver functions return a PiecewiseSolution `sol`; sol(t) gives [x1,x2,x3,x4].
Plotting functions return matplotlib Figure objects.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from core import local_field, release_times, solve_impulsive, subdivide, unpack
from network_experiments import solve_network_data

LABELS = ["Diaphorina citri", "Plant shoots", "Tamarixia radiata", "Plant vigor"]

# Approximation of Mathematica's "TemperatureMap" color scheme.
TEMPERATURE_MAP = LinearSegmentedColormap.from_list(
    "TemperatureMap",
    ["#1f3fbf", "#3f8fdf", "#9fd3ef", "#f2f2f2", "#f7c68a", "#e8733a", "#c02a1d"])


def _style(ax, xlabel, ylabel, top_pad=0.3):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, linestyle="--", alpha=0.6)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + top_pad * (hi - lo))


def _periodic_jump(delta, T, tt, var=2):
    def jump(y):
        y = y.copy()
        y[var] += delta
        return y
    return release_times(T, tt), jump


# ---------------------------------------------------------------------------
# Differential Equation Model
# ---------------------------------------------------------------------------

def chart_report(init, param1, param2, param3, param4):
    """ChartReport: table of initial conditions and parameters."""
    cols = {
        "init": dict(zip(["x1[0]", "x2[0]", "x3[0]", "x4[0]", "t"], init)),
        "x1": dict(zip(["r1", "b12", "b13", "a1", "c1"], param1)),
        "x2": dict(zip(["r2", "b21", "b24", "a2", "c2"], param2)),
        "x3": dict(zip(["r3", "b31", "a3", "c3"], param3)),
        "x4": dict(zip(["r4", "b42", "a4", "c4"], param4)),
    }
    return pd.DataFrame({k: pd.Series(list(v.keys())) + " = " + pd.Series([str(x) for x in v.values()])
                         for k, v in cols.items()}).fillna("")


def solve_system(delta, T, params, init):
    """SolveSystem: single node with periodic Tamarixia releases."""
    Di0, B0, x30, x40, tt = init
    rhs = lambda t, y: np.array(local_field(params, *y))
    jumps, jump = _periodic_jump(delta, T, tt)
    return solve_impulsive(rhs, [Di0, B0, x30, x40], (0.0, float(tt)), jumps, jump)


def system_evolve(delta, T, params, init, label_position=(0.97, 0.95), n_plot=2000):
    sol = solve_system(delta, T, params, init)
    tt = init[-1]
    ts = np.linspace(0, tt, n_plot + 1)
    Y = sol(ts)
    fig, ax = plt.subplots(figsize=(10, 6))
    for k in range(4):
        ax.plot(ts, Y[k], label=LABELS[k])
    _style(ax, "time (a.u.)", "population")
    ax.legend(loc="upper right", bbox_to_anchor=label_position)
    return fig


# ---------------------------------------------------------------------------
# Network version (test)
# ---------------------------------------------------------------------------

def system_evolve_network(graph, Dd, Dt, delta, T, params, init, label_position=(0.97, 0.95),
                          n_plot=2000):
    """SystemEvolveNetwork: one figure per node (releases at every node)."""
    data = solve_network_data(graph, Dd, Dt, delta, T, params, init, "All")
    n, tt = data["N"], init[-1]
    ts = np.linspace(0, tt, n_plot + 1)
    Y = data["Solution"](ts)
    figs = []
    for i in range(n):
        fig, ax = plt.subplots(figsize=(10, 6))
        for k in range(4):
            ax.plot(ts, Y[k * n + i], label=LABELS[k])
        _style(ax, "time (a.u.)", "population")
        ax.legend(loc="upper right", bbox_to_anchor=label_position)
        ax.set_title(f"Node #{i + 1}")
        figs.append(fig)
    return figs


# ---------------------------------------------------------------------------
# Release Schedule Sensitivity: (T, delta) sweep
# ---------------------------------------------------------------------------

RELEASE_SWEEP_METRIC_LABELS = {
    "MeanX1": "Mean Diaphorina citri density",
    "PestDays": "Pest-days (\u222b x1 dt)",
    "FractionBelowThreshold": "Fraction of time below threshold",
}


def release_sweep(T_range, delta_range, params, init,
                  metrics=("MeanX1", "PestDays", "FractionBelowThreshold"),
                  window=None, threshold=10, samples=300):
    """ReleaseSweep: solves once per (T,delta) and computes all metrics.
    Returns {metric: [(T, delta, value), ...]}; value is None on failure.
    Loop order as in Mathematica: delta outer, T inner."""
    tt = init[-1]
    win = (tt / 2, tt) if window is None else window
    per_point = []
    for delta in delta_range:
        for T in T_range:
            try:
                sol = solve_system(delta, T, params, init)
                vals = sol(subdivide(win[0], win[1], samples))[0]
                vals_full = sol(subdivide(0, tt, samples))[0]
                here = {"MeanX1": float(np.mean(vals)),
                        "PestDays": float(np.mean(vals_full) * tt),
                        "FractionBelowThreshold": float(np.count_nonzero(vals < threshold) / len(vals))}
            except Exception:
                here = {k: None for k in ("MeanX1", "PestDays", "FractionBelowThreshold")}
            per_point.append((T, delta, here))
    return {m: [(p[0], p[1], p[2][m]) for p in per_point] for m in metrics}


def _numeric(data):
    return [row for row in data if row[2] is not None and np.isfinite(row[2])]


def release_sweep_plot(result, log_color=(), highlight=None, show_samples=True, boundary=None):
    """ReleaseSweepPlot: one linearly interpolated density plot per metric."""
    out = {}
    for metric, data in result.items():
        clean = np.array(_numeric(data), dtype=float)
        z = np.log10(clean[:, 2]) if metric in log_color else clean[:, 2]
        tri = mtri.Triangulation(clean[:, 0], clean[:, 1])
        fig, ax = plt.subplots(figsize=(7.5, 6))
        pc = ax.tripcolor(tri, z, shading="gouraud", cmap=TEMPERATURE_MAP)
        cb = fig.colorbar(pc, ax=ax)
        if metric in log_color:
            cb.set_label("log$_{10}$(value)")
        if show_samples:
            ax.plot(clean[:, 0], clean[:, 1], ".", color="0.25", ms=3)
        if boundary is not None:
            b = np.array(_numeric(boundary), dtype=float)
            if len(b) >= 3:
                ax.tricontour(b[:, 0], b[:, 1], b[:, 2], levels=[0], colors="k", linewidths=2)
        if highlight is not None:
            ax.plot(*highlight, marker="*", color="k", ms=18)
        ax.set_xlabel("Release period T")
        ax.set_ylabel("Release magnitude \u03b4")
        out[metric] = fig
    return out


# ---------------------------------------------------------------------------
# Scenario Comparison Metrics (no control / single introduction / periodic release)
# ---------------------------------------------------------------------------

def scenario_trajectories(params, init, scenarios):
    """scenarios: list of (name, delta, T) or (name, delta, T, x30override).
    Returns {name: x1(t)}."""
    out = {}
    for scn in scenarios:
        name, delta, T = scn[0], scn[1], scn[2]
        this_init = list(init)
        if len(scn) >= 4:
            this_init[2] = scn[3]
        sol = solve_system(delta, T, params, this_init)
        out[name] = (lambda s: (lambda t: s(t)[0]))(sol)
    return out


def scenario_metrics(params, init, scenarios, recovery_fraction=0.9, samples=1000, window=None):
    """ScenarioMetrics. MeanDensity/PctReductionMeanDensity use `window`
    (default [tt/2,tt]); PestDays and RecoveryTime use the full horizon.
    Missing values (Mathematica's Missing[...]) are returned as None."""
    tt = init[-1]
    win = (tt / 2, tt) if window is None else window
    names = [s[0] for s in scenarios]
    sols = scenario_trajectories(params, init, scenarios)
    ref_name = names[0]
    ref = sols[ref_name]
    ts_full = subdivide(0, tt, samples)
    ts_win = subdivide(win[0], win[1], samples)
    ref_full, ref_win = ref(ts_full), ref(ts_win)
    rows = []
    for k, name in enumerate(names):
        x1f = sols[name]
        vals_full, vals_win = x1f(ts_full), x1f(ts_win)
        mean_win = float(np.mean(vals_win))
        pest_days = float(np.mean(vals_full) * tt)
        delta_here = scenarios[k][1]
        t_release = scenarios[k][2] if delta_here > 0 else 0
        if name == ref_name:
            rec_t = 0
        else:
            post = ts_full > t_release
            post_ts = ts_full[post]
            ratio = vals_full[post] / ref_full[post]
            if len(ratio) == 0:
                rec_t = None
            else:
                dip = int(np.argmin(ratio))
                hits = np.nonzero(ratio[dip:] >= recovery_fraction)[0]
                rec_t = None if len(hits) == 0 else float(post_ts[dip + hits[0]])
        rows.append({"Scenario": name, "MeanDensity": mean_win, "PestDays": pest_days,
                     "RecoveryTime": rec_t,
                     "PctReductionMeanDensity": 0.0 if name == ref_name
                     else float(100 * (1 - mean_win / np.mean(ref_win)))})
    return pd.DataFrame(rows)


def scenario_comparison_plot(params, init, scenarios, n_plot=2000):
    sols = scenario_trajectories(params, init, scenarios)
    tt = init[-1]
    ts = np.linspace(0, tt, n_plot + 1)
    fig, ax = plt.subplots(figsize=(10, 6))
    for name, f in sols.items():
        ax.plot(ts, f(ts), label=name)
    _style(ax, "time (a.u.)", "Diaphorina citri")
    ax.legend(loc="upper right", bbox_to_anchor=(0.98, 0.97))
    return fig


# ---------------------------------------------------------------------------
# Self-Limitation Term Boundedness / Positivity Check
# ---------------------------------------------------------------------------

def _self_limitation_coefs(params, Y):
    r1, b12, b13, a1, c1, r2, b21, b24, a2, c2, r3, b31, a3, c3, r4, b42, a4, c4 = unpack(params)
    y1, y2, y3, y4 = Y
    return [a1 + c1 * (b12 * y2 + b13 * y3), a2 + c2 * (b21 * y1 + b24 * y4),
            a3 + c3 * b31 * y1, a4 + c4 * b42 * y2]


def self_limitation_check(delta, T, params, init, samples=2000):
    tt = init[-1]
    sol = solve_system(delta, T, params, init)
    ts = subdivide(0, tt, samples)
    coefs = _self_limitation_coefs(params, sol(ts))
    out = {}
    for k, c in enumerate(coefs):
        out[f"x{k + 1}"] = {"MinCoefficient": float(c.min()),
                            "ArgMinTime": float(ts[int(np.argmin(c))]),
                            "StaysPositive": bool(c.min() > 0),
                            "FractionNegative": float(np.count_nonzero(c < 0) / len(c))}
    return out


def self_limitation_plot(delta, T, params, init, n_plot=2000):
    tt = init[-1]
    sol = solve_system(delta, T, params, init)
    ts = np.linspace(0, tt, n_plot + 1)
    coefs = _self_limitation_coefs(params, sol(ts))
    fig, ax = plt.subplots(figsize=(10, 6))
    for k, c in enumerate(coefs):
        ax.plot(ts, c, label="coef. " + LABELS[k])
    ax.axhline(0, color="red", linestyle="--")
    _style(ax, "time (a.u.)",
           r"self-limitation coefficient  $a_i + c_i\,\sigma\beta\,x_j$", top_pad=0.35)
    ax.legend(loc="upper right", bbox_to_anchor=(0.98, 0.97))
    return fig


# ---------------------------------------------------------------------------
# Saturating (Holling type II) Predation Response
# ---------------------------------------------------------------------------

def solve_system_holling(delta, T, params, ks, init):
    r1, b12, b13, a1, c1, r2, b21, b24, a2, c2, r3, b31, a3, c3, r4, b42, a4, c4 = unpack(params)
    k13, k31 = ks
    Di0, B0, x30, x40, tt = init

    def rhs(t, y):
        y1, y2, y3, y4 = y
        h13 = b13 * y3 / (1 + k13 * y3)
        h31 = b31 * y1 / (1 + k31 * y1)
        return np.array([
            y1 * (r1 + b12 * y2 + h13 - (a1 + c1 * (b12 * y2 + h13)) * y1),
            y2 * (r2 + b21 * y1 + b24 * y4 - (a2 + c2 * (b21 * y1 + b24 * y4)) * y2),
            y3 * (r3 + h31 - (a3 + c3 * h31) * y3),
            y4 * (r4 + b42 * y2 - (a4 + c4 * b42 * y2) * y4)])
    jumps, jump = _periodic_jump(delta, T, tt)
    return solve_impulsive(rhs, [Di0, B0, x30, x40], (0.0, float(tt)), jumps, jump)


def holling_comparison_plot(delta, T, params, ks, init, n_plot=2000):
    tt = init[-1]
    base = solve_system(delta, T, params, init)
    holl = solve_system_holling(delta, T, params, ks, init)
    ts = np.linspace(0, tt, n_plot + 1)
    B, H = base(ts), holl(ts)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(ts, B[0], "b-", label="Diaphorina citri mass-action (baseline)")
    ax.plot(ts, H[0], "b--", label=f"Diaphorina citri Holling II ($k_{{13}}$={ks[0]}, $k_{{31}}$={ks[1]})")
    ax.plot(ts, B[2], "g-", label="Tamarixia radiata mass-action (baseline)")
    ax.plot(ts, H[2], "g--", label="Tamarixia radiata Holling II")
    _style(ax, "time (a.u.)", "population", top_pad=0.35)
    ax.legend(loc="upper right", bbox_to_anchor=(0.98, 0.97))
    return fig


# ---------------------------------------------------------------------------
# Pest-Free Periodic Solution & Invasion Exponent
# ---------------------------------------------------------------------------

def solve_pest_free_subsystem(delta, T, params, init3):
    """Returns sol with sol(t) = [x2, x3, x4]."""
    r1, b12, b13, a1, c1, r2, b21, b24, a2, c2, r3, b31, a3, c3, r4, b42, a4, c4 = unpack(params)
    B0, x30, x40, tt = init3

    def rhs(t, y):
        y2, y3, y4 = y
        return np.array([y2 * (r2 + b24 * y4 - (a2 + c2 * b24 * y4) * y2),
                         y3 * (r3 - a3 * y3),
                         y4 * (r4 + b42 * y2 - (a4 + c4 * b42 * y2) * y4)])
    jumps, jump = _periodic_jump(delta, T, tt, var=1)
    return solve_impulsive(rhs, [B0, x30, x40], (0.0, float(tt)), jumps, jump)


def invasion_exponent(params, init3, delta, T, total_time=None, burn_in_fraction=0.3, samples=3000):
    tt = max(400, 60 * T) if total_time is None else total_time
    burn = burn_in_fraction * tt
    sol = solve_pest_free_subsystem(delta, T, params, list(init3) + [tt])
    r1, b12, b13 = params[0][0], params[0][1], params[0][2]
    ts = subdivide(burn, tt, samples)
    Y = sol(ts)
    integrand = r1 + b12 * Y[0] + b13 * Y[1]
    lam = float(np.mean(integrand))
    h = len(integrand) // 2
    conv = float(abs(np.mean(integrand[:h]) - np.mean(integrand[h:])))
    return {"InvasionExponent": lam, "PestFreeStable": lam < 0, "TotalTime": tt,
            "BurnIn": burn, "ConvergenceGap": conv, "Reliable": conv < 0.02}


def invasion_exponent_sweep(T_range, delta_range, params, init3, **opts):
    out = []
    for delta in delta_range:
        for T in T_range:
            try:
                out.append((T, delta, invasion_exponent(params, init3, delta, T, **opts)["InvasionExponent"]))
            except Exception:
                out.append((T, delta, None))
    return out


def invasion_exponent_sweep_plot(data, highlight=None, show_samples=True):
    """Symmetric Blue-White-Red scale centered at zero, with the zero contour."""
    clean = np.array(_numeric(data), dtype=float)
    max_abs = float(np.max(np.abs(clean[:, 2])))
    cmap = LinearSegmentedColormap.from_list("BWR", ["blue", "white", "red"])
    tri = mtri.Triangulation(clean[:, 0], clean[:, 1])
    fig, ax = plt.subplots(figsize=(7.5, 6))
    pc = ax.tripcolor(tri, clean[:, 2], shading="gouraud", cmap=cmap,
                      norm=TwoSlopeNorm(0.0, -max_abs, max_abs))
    fig.colorbar(pc, ax=ax)
    ax.tricontour(tri, clean[:, 2], levels=[0], colors="k", linewidths=2)
    if show_samples:
        ax.plot(clean[:, 0], clean[:, 1], ".", color="0.25", ms=3)
    if highlight is not None:
        ax.plot(*highlight, marker="*", color="k", ms=18)
    ax.set_xlabel("Release period T")
    ax.set_ylabel("Release magnitude \u03b4")
    return fig


# ---------------------------------------------------------------------------
# Robustness to Imperfect Release Synchrony (jitter)
# ---------------------------------------------------------------------------

def solve_system_jitter(delta, T, jitter, params, init, seed=1):
    """Release k happens at k*T + U(-jitter, jitter). numpy's RNG replaces
    Mathematica's, so the realization for a given seed differs from Mathematica's."""
    Di0, B0, x30, x40, tt = init
    rng = np.random.default_rng(seed)
    n_rel = int(np.floor(tt / T))
    times = T * np.arange(1, n_rel + 1) + rng.uniform(-jitter, jitter, n_rel)
    times = [t for t in times if 0 < t < tt]
    rhs = lambda t, y: np.array(local_field(params, *y))

    def jump(y):
        y = y.copy()
        y[2] += delta
        return y
    return solve_impulsive(rhs, [Di0, B0, x30, x40], (0.0, float(tt)), times, jump)


def jitter_realizations(delta, T, jitter, params, init, n_seeds, n_samp):
    ts = subdivide(0, init[-1], n_samp)
    vals = np.array([solve_system_jitter(delta, T, jitter, params, init, seed)(ts)[0]
                     for seed in range(1, n_seeds + 1)])
    return {"Times": ts, "Values": vals}


def jitter_robustness_stats(delta, T, jitter, params, init, seeds=50, samples=400, window=None):
    tt = init[-1]
    win = (tt / 2, tt) if window is None else window
    r = jitter_realizations(delta, T, jitter, params, init, seeds, samples)
    ts, V = r["Times"], r["Values"]
    in_win = (ts >= win[0]) & (ts <= win[1])
    means = V[:, in_win].mean(axis=1)
    return {"Seeds": seeds, "Window": win, "MeanOfMeans": float(means.mean()),
            "SDOfMeans": float(means.std(ddof=1))}


def jitter_robustness_plot(delta, T, jitter, params, init, seeds=50, samples=400):
    exact = solve_system(delta, T, params, init)
    r = jitter_realizations(delta, T, jitter, params, init, seeds, samples)
    ts, V = r["Times"], r["Values"]
    m, sd = V.mean(axis=0), V.std(axis=0, ddof=1)
    fig, ax = plt.subplots(figsize=(10, 10 / 1.618))
    ax.fill_between(ts, m - sd, m + sd, color="gray", alpha=0.4, linewidth=0)
    ax.plot(ts, m, color="gray", lw=2,
            label=f"Jitter: mean \u00b1 1 s.d. (N={seeds}, J={jitter})")
    ax.plot(ts, exact(ts)[0], color="black", lw=2, label=f"Exact period (T={T})")
    _style(ax, "time (a.u.)", r"Diaphorina citri ($x_1$)", top_pad=0.28)
    ax.legend(loc="upper right", bbox_to_anchor=(0.98, 0.97))
    return fig
