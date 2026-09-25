# Model/Python

Versión en Python de `Model/Mathematica`.

| Mathematica | Python |
|---|---|
| `Functions.m` | `functions.py` |
| `NetworkExperiments.wl` | `network_experiments.py` |
| `NetworkStability.wl` | `network_stability.py` |
| `RunExperimentK_*.wls` | `run_experimentK_*.py` |
| `PlotFromCSV_ExperimentK.wls` | `plot_from_csv_experimentK.py` |
| `Model.nb` | `Model.ipynb` |

`core.py` contiene el campo vectorial local, su Jacobiano y el integrador impulsivo compartido. `common.py` contiene los parámetros de `tabla_latex.txt` y la exportación a CSV.

## Instalación y uso

```bash
cd Model/Python
pip install -r requirements.txt
python run_experiment2_outbreak_containment.py
```

Los resultados se escriben en `Model/Python/Results/`. Orden para reproducir todo desde cero: 1 → 2 → 3 → 5 → 4 → 6 → 7. El 4 lee el CSV resumen del 3.

## Convenciones

Los nodos se indexan desde 0 (el nodo 1 de Mathematica es el nodo 0). En la estrategia `TargetedAtSource` la liberación se aplica a `[0]`.

Las liberaciones se resuelven integrando por tramos entre instantes `kT` y aplicando el salto entre tramos. Evaluada justo en un instante de liberación, la solución devuelve el estado posterior al salto. Una liberación que cae exactamente en `t = tt` no se aplica porque no hay tiempo posterior en el que actúe.

Integrador: DOP853 con `rtol = 1e-10`, `atol = 1e-12` (modificables en `core.py`).

Las topologías Random, SmallWorld y ScaleFree se generan con el RNG de numpy/networkx. Para una misma semilla son instancias distintas de las de Mathematica. Lo mismo vale para las realizaciones de jitter.
