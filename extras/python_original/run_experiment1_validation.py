"""Sanity check: with Dd=Dt=0, every node of the network model must evolve
exactly like the single-node model, since the diffusive coupling is switched
off. Validates that solve_network_data reduces correctly to the base ODE system
before it is used for experiments."""

import numpy as np

from common import PARAMS, P1, P2, P3, P4
from core import solve_impulsive, subdivide
from network_experiments import make_topology, node_if, solve_network_data, x1, x3

tt = 100
n = 3
ics = [[50] * n, [10] * n, [15] * n, [30] * n, tt]


# Reference: single-node solution, no releases (delta=0), written out explicitly
def ref_rhs(t, y):
    x1s, x2s, x3s, x4s = y
    return np.array([
        x1s * (P1[0] + P1[1] * x2s + P1[2] * x3s - (P1[3] + P1[4] * (P1[1] * x2s + P1[2] * x3s)) * x1s),
        x2s * (P2[0] + P2[1] * x1s + P2[2] * x4s - (P2[3] + P2[4] * (P2[1] * x1s + P2[2] * x4s)) * x2s),
        x3s * (P3[0] + P3[1] * x1s - (P3[2] + P3[3] * P3[1] * x1s) * x3s),
        x4s * (P4[0] + P4[1] * x2s - (P4[2] + P4[3] * P4[1] * x2s) * x4s)])


ref_sol = solve_impulsive(ref_rhs, [50, 10, 15, 30], (0.0, float(tt)))

# Network run: 3 nodes on a Path graph, Dd=Dt=0, no releases
graph = make_topology("Path", n)
net_data = solve_network_data(graph, 0, 0, 0, 15, PARAMS, ics)

tvals = subdivide(0.0, tt, 200)
ref_vals = ref_sol(tvals)
max_err_x1 = max(np.max(np.abs(node_if(net_data, x1, i)(tvals) - ref_vals[0])) for i in range(n))
max_err_x3 = max(np.max(np.abs(node_if(net_data, x3, i)(tvals) - ref_vals[2])) for i in range(n))

# All nodes should also be pairwise identical to each other
max_cross = max(np.max(np.abs(node_if(net_data, x1, 0)(tvals) - node_if(net_data, x1, i)(tvals)))
                for i in range(1, n))

report = {
    "MaxErrorX1_vs_SingleNode": float(max_err_x1),
    "MaxErrorX3_vs_SingleNode": float(max_err_x3),
    "MaxCrossNodeError": float(max_cross),
    "Pass": bool(max_err_x1 < 1e-6 and max_err_x3 < 1e-6 and max_cross < 1e-9),
}

print("Validation report:", report)
if report["Pass"]:
    print("PASS: network model reduces to single-node model when Dd=Dt=0.")
else:
    print("FAIL: network model does NOT reduce to single-node model. Check solve_network_data.")
