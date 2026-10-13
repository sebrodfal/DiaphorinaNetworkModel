"""Shared numerical core for the Python port of the DiaphorinaNetworkModel.

Replaces NDSolve[..., WhenEvent[Mod[t,T]==0 && t>0, x3 -> x3 + delta],
Method -> {"TimeIntegration" -> "ExplicitRungeKutta"}].

Impulses are handled by integrating piecewise between release times and
applying the jump between pieces (exactly equivalent to WhenEvent, without
relying on event detection). Every piecewise solution is RIGHT-CONTINUOUS:
evaluated exactly at a release time it returns the post-release state.
"""

import numpy as np
from scipy.integrate import solve_ivp

# Integrator settings. Mathematica's defaults are roughly 8 significant digits;
# these are tighter so the Python results are converged to ~10 digits.
RTOL = 1e-10
ATOL = 1e-12
METHOD = "DOP853"

# Gauss-Legendre nodes used to integrate the dense output exactly over each
# internal step (DOP853's interpolant is a degree-7 polynomial per step).
_GL_X, _GL_W = np.polynomial.legendre.leggauss(8)


# ---------------------------------------------------------------------------
# Local vector field (one node) and its Jacobian
# ---------------------------------------------------------------------------

def unpack(params):
    (r1, b12, b13, a1, c1), (r2, b21, b24, a2, c2), (r3, b31, a3, c3), (r4, b42, a4, c4) = params
    return r1, b12, b13, a1, c1, r2, b21, b24, a2, c2, r3, b31, a3, c3, r4, b42, a4, c4


def local_field(params, y1, y2, y3, y4):
    """LocalField[params][{y1,y2,y3,y4}]. Works elementwise on arrays."""
    r1, b12, b13, a1, c1, r2, b21, b24, a2, c2, r3, b31, a3, c3, r4, b42, a4, c4 = unpack(params)
    return (
        y1 * (r1 + b12 * y2 + b13 * y3 - (a1 + c1 * (b12 * y2 + b13 * y3)) * y1),
        y2 * (r2 + b21 * y1 + b24 * y4 - (a2 + c2 * (b21 * y1 + b24 * y4)) * y2),
        y3 * (r3 + b31 * y1 - (a3 + c3 * b31 * y1) * y3),
        y4 * (r4 + b42 * y2 - (a4 + c4 * b42 * y2) * y4),
    )


def local_jacobian_function(params):
    """LocalJacobianFunction[params]: returns J(y1,y2,y3,y4) -> 4x4 array."""
    r1, b12, b13, a1, c1, r2, b21, b24, a2, c2, r3, b31, a3, c3, r4, b42, a4, c4 = unpack(params)

    def J(y1, y2, y3, y4):
        return np.array([
            [r1 + b12 * y2 + b13 * y3 - 2 * (a1 + c1 * (b12 * y2 + b13 * y3)) * y1,
             b12 * y1 - c1 * b12 * y1 ** 2,
             b13 * y1 - c1 * b13 * y1 ** 2,
             0.0],
            [b21 * y2 - c2 * b21 * y2 ** 2,
             r2 + b21 * y1 + b24 * y4 - 2 * (a2 + c2 * (b21 * y1 + b24 * y4)) * y2,
             0.0,
             b24 * y2 - c2 * b24 * y2 ** 2],
            [b31 * y3 - c3 * b31 * y3 ** 2,
             0.0,
             r3 + b31 * y1 - 2 * (a3 + c3 * b31 * y1) * y3,
             0.0],
            [0.0,
             b42 * y4 - c4 * b42 * y4 ** 2,
             0.0,
             r4 + b42 * y2 - 2 * (a4 + c4 * b42 * y2) * y4],
        ])

    return J


def single_node_rhs(params):
    def f(t, y):
        return np.array(local_field(params, y[0], y[1], y[2], y[3]))
    return f


# ---------------------------------------------------------------------------
# Piecewise (impulsive) solution
# ---------------------------------------------------------------------------

class PiecewiseSolution:
    """Continuous-time solution made of smooth pieces separated by jumps.

    sol(t) -> array of shape (dim,) for scalar t, (dim, len(t)) for array t.
    Right-continuous at the jumps.
    """

    def __init__(self, pieces, t0, t1):
        self.pieces = pieces            # list of (ta, tb, OdeSolution)
        self.starts = np.array([p[0] for p in pieces])
        self.t0, self.t1 = t0, t1
        self.dim = pieces[0][2](pieces[0][0]).shape[0]

    def __call__(self, t):
        scalar = np.ndim(t) == 0
        tv = np.atleast_1d(np.asarray(t, dtype=float))
        idx = np.clip(np.searchsorted(self.starts, tv, side="right") - 1, 0, len(self.pieces) - 1)
        out = np.empty((self.dim, tv.size))
        for k in np.unique(idx):
            m = idx == k
            ta, tb, so = self.pieces[k]
            out[:, m] = so(np.clip(tv[m], ta, tb))
        return out[:, 0] if scalar else out

    def component(self, i):
        return lambda t: self(t)[i]

    def integrate(self, a, b, comps=None):
        """Exact integral of the interpolant over [a,b] (NIntegrate of the
        InterpolatingFunction), for all components or the selected ones."""
        total = 0.0
        for ta, tb, so in self.pieces:
            lo, hi = max(a, ta), min(b, tb)
            if hi <= lo:
                continue
            ts = np.asarray(so.ts)
            if ts[0] > ts[-1]:
                ts = ts[::-1]
            grid = np.concatenate(([lo], ts[(ts > lo) & (ts < hi)], [hi]))
            left, right = grid[:-1], grid[1:]
            half = 0.5 * (right - left)
            mid = 0.5 * (right + left)
            nodes = (mid[:, None] + half[:, None] * _GL_X[None, :]).ravel()
            weights = (half[:, None] * _GL_W[None, :]).ravel()
            vals = so(nodes)
            if comps is not None:
                vals = vals[comps]
            total = total + vals @ weights
        return float(total) if np.ndim(total) == 0 else total


def release_times(T, tt):
    """Times kT (k>=1) strictly inside (0, tt): where Mod[t,T]==0 && t>0 fires.
    A release falling exactly on t=tt has no effect on the returned solution
    (there is no time after it), so it is not applied."""
    if T is None or T <= 0:
        return []
    k = 1
    out = []
    while k * T < tt:
        out.append(k * T)
        k += 1
    return out


def n_releases(T, tt):
    """Number of releases actually applied on (0, tt): n_rel(T) = ceil(tt/T) - 1
    (plan, Fase 0, cambio 3). floor(tt/T) over-counts by one whenever T divides
    tt, because the release at t = tt is never applied (see release_times)."""
    if T is None or T <= 0:
        return 0
    n = int(np.ceil(tt / T)) - 1
    assert n == len(release_times(T, tt)), (T, tt)
    return n


def solve_impulsive(rhs, y0, t_span, jump_times=(), jump=None,
                    rtol=None, atol=None, method=None):
    """Integrates y' = rhs(t,y) over t_span, applying y -> jump(y) at each of
    jump_times (sorted, strictly inside t_span)."""
    rtol = RTOL if rtol is None else rtol
    atol = ATOL if atol is None else atol
    method = METHOD if method is None else method
    t0, t1 = t_span
    y = np.array(y0, dtype=float)
    cuts = [tj for tj in sorted(jump_times) if t0 < tj < t1]
    pieces = []
    ta = t0
    for tb in cuts + [t1]:
        s = solve_ivp(rhs, (ta, tb), y, method=method, rtol=rtol, atol=atol, dense_output=True)
        if not s.success:
            raise RuntimeError(f"Integration failed on [{ta},{tb}]: {s.message}")
        pieces.append((ta, tb, s.sol))
        y = s.y[:, -1].copy()
        if tb < t1 and jump is not None:
            y = jump(y)
        ta = tb
    return PiecewiseSolution(pieces, t0, t1)


def subdivide(a, b, n):
    """Subdivide[a,b,n]: n+1 equally spaced points."""
    return np.linspace(a, b, int(n) + 1)
