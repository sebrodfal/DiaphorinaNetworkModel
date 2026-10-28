"""sigma / (m_t <k>) for the eight graphs of Experiment 10 (section 4, equation (8)).
Writes Results/stability/sigma_ratio.csv."""
import os

import networkx as nx
import numpy as np
import pandas as pd

G = {
    "path": nx.path_graph(8), "cycle": nx.cycle_graph(8),
    "random s9": nx.gnp_random_graph(8, 0.4, seed=9),
    "random s3": nx.gnp_random_graph(8, 0.4, seed=3),
    "random s4": nx.gnp_random_graph(8, 0.4, seed=4),
    "scale-free s2": nx.barabasi_albert_graph(8, 2, seed=2),
    "star": nx.star_graph(7), "complete": nx.complete_graph(8),
}
rows = []
for name, g in G.items():
    A = nx.to_numpy_array(g, nodelist=sorted(g))
    K = np.diag(A.sum(1))
    kmean = A.sum(1).mean()
    for mt in (0.05, 0.1, 0.2, 0.4):
        sigma = np.linalg.eigvalsh(K - (1 - mt) * A).min()
        rows.append({"graph": name, "k_s": int(A[0].sum()), "mean_k": kmean,
                     "m_t": mt, "sigma": sigma, "ratio": sigma / (mt * kmean)})
df = pd.DataFrame(rows)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Results", "stability")
os.makedirs(OUT, exist_ok=True)
df.to_csv(os.path.join(OUT, "sigma_ratio.csv"), index=False)
print(df.pivot(index="graph", columns="m_t", values="ratio").round(3))
