"""Linear (Floquet / Master Stability) analysis of the synchronized state of the
DiaphorinaModel network extension (port of NetworkStability.wl).

With a Uniform release strategy and identical local parameters at every node,
the synchronized state x_i(t) = s(t) is an exact solution of the network system
(Laplacian rows sum to zero). Projecting the variational equation onto the
Laplacian eigenbasis decouples it into n independent 4x4 problems, one per
eigenvalue mu_k:

    Y_k'(t) = (Df(s(t)) - mu_k * D) Y_k(t),   D = diag(Dd, 0, Dt, 0)

mu_1 = 0 is the synchronized mode itself; mu_2..mu_n are the transverse modes.
The jump map (x3 -> x3 + delta) has identity Jacobian, so it does not enter.

master_monodromy / master_stability_value / master_stability_curve require a
genuine period-T orbit (find_periodic_state). master_lyapunov_exponent /
master_lyapunov_spectrum do not: they integrate one perturbation vector along a
long reference trajectory, renormalizing every segment (Benettin), and report
the finite-time Lyapunov exponent of the mode.

Convention at the release instants: the reference trajectory is right-
continuous, so a segment starting at a release time uses the post-release
state. Mathematica's InterpolatingFunction is ambiguous exactly at the jump,
which is why its exponents differ from these (see chat notes).
"""

import numpy as np

from core import (local_field, local_jacobian_function, release_times, solve_impulsive,
                  single_node_rhs, METHOD, RTOL, ATOL)
from scipy.integrate import solve_ivp
from network_experiments import laplacian, network_rhs


# ---------- Reference periodic orbit (single node, converged) ----------

def _single_node_solution(params, delta, T, t_total, init):
    jumps = [] if delta == 0 else release_times(T, t_total)

    def jump(y):
        y = y.copy()
        y[2] += delta
        return y
    return solve_impulsive(single_node_rhs(params), init, (0.0, float(t_total)), jumps, jump)


def find_periodic_state(params, delta, T, n_periods=25):
    t_burn = n_periods * T
    sol = _single_node_solution(params, delta, T, t_burn, [50, 10, 15, 30])
    s0 = sol(t_burn - T + 1e-6)
    s_prev = sol(t_burn - 2 * T + 1e-6)
    return {"s0": s0, "ConvergenceGap": float(np.linalg.norm(s0 - s_prev)),
            "Solution": sol, "tBurn": t_burn}


# ---------- Master stability: Y_k(T) for a single scalar mu ----------

def master_monodromy(params, Dd, Dt, mu, s0, T):
    Jf = local_jacobian_function(params)
    Dmat = np.diag([Dd, 0.0, Dt, 0.0])

    def rhs(t, z):
        s = z[:4]
        Y = z[4:].reshape(4, 4)
        ds = np.array(local_field(params, *s))
        dY = (Jf(*s) - mu * Dmat) @ Y
        return np.concatenate((ds, dY.ravel()))

    z0 = np.concatenate((np.asarray(s0, float), np.eye(4).ravel()))
    r = solve_ivp(rhs, (0.0, T), z0, method=METHOD, rtol=RTOL, atol=ATOL)
    return r.y[4:, -1].reshape(4, 4)


def master_stability_value(params, Dd, Dt, mu, s0, T):
    return float(np.max(np.abs(np.linalg.eigvals(master_monodromy(params, Dd, Dt, mu, s0, T)))))


def master_stability_curve(params, Dd, Dt, s0, T, mu_values):
    return [{"mu": mu, "SpectralRadius": master_stability_value(params, Dd, Dt, mu, s0, T)}
            for mu in mu_values]


# ---------- Master stability without assuming periodicity ----------

def reference_trajectory(params, delta, T, t_total, init=(50, 10, 15, 30)):
    """Single-node trajectory over [0, t_total]: the path the whole synchronized
    network would follow. Returns a right-continuous callable s(t) -> (4,)."""
    return _single_node_solution(params, delta, T, t_total, list(init))


def perturbation_segment(Jf, Dmat, mu, s_fn, ta, tb, v0):
    """Integrates the transverse variational equation for one vector over [ta,tb]."""
    def rhs(t, v):
        return (Jf(*s_fn(t)) - mu * Dmat) @ v
    r = solve_ivp(rhs, (ta, tb), np.asarray(v0, float), method=METHOD, rtol=RTOL, atol=ATOL)
    return r.y[:, -1]


def lyapunov_exponent_from_reference(Jf, Dmat, mu, s_fn, burn, seg_len, n_seg, v0,
                                     convergence_threshold=0.02):
    t_c = burn
    v = np.asarray(v0, float)
    v = v / np.linalg.norm(v)
    log_sums = [0.0]
    for _ in range(n_seg):
        v_end = perturbation_segment(Jf, Dmat, mu, s_fn, t_c, t_c + seg_len, v)
        nrm = np.linalg.norm(v_end)
        v = v_end / nrm
        log_sums.append(log_sums[-1] + np.log(nrm))
        t_c += seg_len
    lam = log_sums[-1] / (n_seg * seg_len)
    half = n_seg // 2
    lam1 = log_sums[half] / (half * seg_len) if half >= 1 else None
    lam2 = (log_sums[-1] - log_sums[half]) / ((n_seg - half) * seg_len) if n_seg - half >= 1 else None
    gap = abs(lam1 - lam2) if (lam1 is not None and lam2 is not None) else None
    return {"LyapunovExponent": float(lam), "ConvergenceGap": gap,
            "Reliable": bool(gap is not None and gap < convergence_threshold)}


def _lyap_setup(params, Dd, Dt, T, total_time, burn_in_fraction, segment_length):
    t_total = max(600, 80 * (T if T > 0 else 10)) if total_time is None else total_time
    burn = burn_in_fraction * t_total
    seg_len = (T if T > 0 else 10) if segment_length is None else segment_length
    n_seg = int(np.floor((t_total - burn) / seg_len))
    return t_total, burn, seg_len, n_seg, np.diag([Dd, 0.0, Dt, 0.0]), local_jacobian_function(params)


def master_lyapunov_exponent(params, Dd, Dt, mu, delta, T, total_time=None, burn_in_fraction=0.3,
                             segment_length=None, initial_perturbation=(1.0, 0.0, 0.0, 0.0),
                             convergence_threshold=0.02, init=(50, 10, 15, 30)):
    t_total, burn, seg_len, n_seg, Dmat, Jf = _lyap_setup(
        params, Dd, Dt, T, total_time, burn_in_fraction, segment_length)
    ref = reference_trajectory(params, delta, T, t_total, init)
    res = lyapunov_exponent_from_reference(Jf, Dmat, mu, ref, burn, seg_len, n_seg,
                                           initial_perturbation, convergence_threshold)
    res.update({"ModeStable": res["LyapunovExponent"] < 0, "TotalTime": t_total,
                "BurnIn": burn, "NumSegments": n_seg, "SegmentLength": seg_len})
    return res


def master_lyapunov_spectrum(params, Dd, Dt, mu_values, delta, T, total_time=None,
                             burn_in_fraction=0.3, segment_length=None,
                             initial_perturbation=(1.0, 0.0, 0.0, 0.0),
                             convergence_threshold=0.02, init=(50, 10, 15, 30)):
    """Solves the reference trajectory once and reuses it for every mu."""
    t_total, burn, seg_len, n_seg, Dmat, Jf = _lyap_setup(
        params, Dd, Dt, T, total_time, burn_in_fraction, segment_length)
    ref = reference_trajectory(params, delta, T, t_total, init)
    out = []
    for mu in mu_values:
        res = lyapunov_exponent_from_reference(Jf, Dmat, mu, ref, burn, seg_len, n_seg,
                                               initial_perturbation, convergence_threshold)
        res.update({"mu": float(mu), "ModeStable": res["LyapunovExponent"] < 0})
        out.append(res)
    return out


# ---------- Brute-force full-network monodromy (finite differences) ----------
# Ground truth used only to validate master_monodromy; O(4n) integrations.

def continuous_flow_network(graph, Dd, Dt, params, X0flat, T):
    L = laplacian(graph)
    r = solve_ivp(network_rhs(L, Dd, Dt, params), (0.0, T), np.asarray(X0flat, float),
                  method=METHOD, rtol=RTOL, atol=ATOL)
    return r.y[:, -1]


def network_monodromy(graph, Dd, Dt, params, X0flat, T, eps=1e-5):
    X0 = np.asarray(X0flat, float)
    dim = X0.size
    F0 = continuous_flow_network(graph, Dd, Dt, params, X0, T)
    cols = [(continuous_flow_network(graph, Dd, Dt, params, X0 + eps * np.eye(dim)[d], T) - F0) / eps
            for d in range(dim)]
    return np.array(cols).T


def network_floquet_multipliers(graph, Dd, Dt, params, X0flat, T, eps=1e-5):
    return np.linalg.eigvals(network_monodromy(graph, Dd, Dt, params, X0flat, T, eps))
