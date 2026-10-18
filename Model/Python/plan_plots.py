"""Common figure style for the plan figures (phases 1-4)."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

R3_COLOR = {-0.2: SERIES[0], -0.3: SERIES[1], -0.45: SERIES[2]}
DEGREE_MARKER = {1: "o", 2: "s", 3: "^", 4: "D", 5: "v", 6: "P", 7: "*"}
NET_ORDER = ["Path", "Cycle", "Random_s9", "Random_s3", "Random_s4", "ScaleFree_s2", "Star", "Complete"]
NET_LABEL = {"Path": "Trayectoria (k=1)", "Cycle": "Ciclo (k=2)", "Random_s9": "Aleatoria s9 (k=3)",
             "Random_s3": "Aleatoria s3 (k=4)", "Random_s4": "Aleatoria s4 (k=5)",
             "ScaleFree_s2": "Libre de escala s2 (k=6)", "Star": "Estrella (k=7)",
             "Complete": "Completo (k=7)"}
NET_COLOR = {n: SERIES[i] for i, n in enumerate(NET_ORDER)}

SEQ = LinearSegmentedColormap.from_list(
    "seq_blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
DIV = LinearSegmentedColormap.from_list(
    "div_blue_red", ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f3a3a2", "#e34948", "#a32524"])


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 150, "font.size": 9.5,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.8, "legend.frameon": False, "axes.titlesize": 10,
        "savefig.bbox": "tight", "figure.facecolor": "white",
    })


style()
