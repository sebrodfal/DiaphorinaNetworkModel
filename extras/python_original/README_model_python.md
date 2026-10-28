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

El experimento 6 se corrigió y regeneró. Ahora cuenta las liberaciones como `ceil(tt/T) - 1`, las que de verdad se aplican (ver `docs/CORRECCIONES.md`).

## Convenciones

Los nodos se indexan desde 0 (el nodo 1 de Mathematica es el nodo 0). En la estrategia `TargetedAtSource` la liberación se aplica a `[0]`.

Las liberaciones se resuelven integrando por tramos entre instantes `kT` y aplicando el salto entre tramos. Evaluada justo en un instante de liberación, la solución devuelve el estado posterior al salto. Una liberación que cae exactamente en `t = tt` no se aplica porque no hay tiempo posterior en el que actúe.

Integrador: DOP853 con `rtol = 1e-10`, `atol = 1e-12` (modificables en `core.py`).

Las topologías Random, SmallWorld y ScaleFree se generan con el RNG de numpy/networkx. Para una misma semilla son instancias distintas de las de Mathematica. Lo mismo vale para las realizaciones de jitter.

## Extensión: mortalidad en tránsito y plan de trabajo (fases 0 a 4)

`network_rhs` y `solve_network_data` aceptan `md` y `mt` (por defecto 0), que son la probabilidad de morir en el trayecto entre nodos del psílido y del parasitoide:

    x1' = f1(x) - Dd (L + md A) x1,      x3' = f3(x) - Dt (L + mt A) x3

`md` y `mt` pueden ser escalares o matrices n×n (m_ij, mortalidad por arista). Con 0 el modelo es idéntico al original. `local=False` apaga la dinámica local y se usa en la prueba de balance de masa. `common.outbreak_ics(n, tt, "A"|"B")` da las dos condiciones iniciales y `core.n_releases(T, tt)` el número de liberaciones aplicadas.

| Archivo | Uso |
|---|---|
| `plan_common.py` | redes focales, Π, Π_t, `run_combo` (3 estrategias + eficiencias) y `run_specs` (ejecutor paralelo reanudable) |
| `plan_plots.py` | estilo y paleta de las figuras |
| `test_phase0.py` | 7 pruebas de la fase 0 |
| `run_phase1_relative_dispersion.py`, `analyze_phase1.py` | fase 1 |
| `run_phase2_goldilocks.py`, `analyze_phase2.py` | fase 2 (`--check-dd0` para la verificación con D_d = 0) |
| `run_phase3_transit_mortality.py`, `analyze_phase3.py` | fases 3a (`--pilot`) y 3b (`--phase3b`) y la parte de AUC de 3c |
| `run_phase3c_stability.py` | 3c, estabilidad del estado sincrónico con m_t > 0 en el ciclo y el grafo completo |
| `run_phase4_edge_mortality.py` | fase 4, grafos geométricos con m_ij = 1 − exp(−d_ij/ℓ) |
| `run_extra_phase3a_extended.py` | fase 3d, la 3a con D_t de 1,78 a 31,6 (se combina con `phase3a_grid.csv`) |
| `run_extra_no_resident.py` | fase 5, liberación sin parasitoide residente (x₃(0) = 0, `CI = B0`) |
| `analyze_extra.py` | tablas, figuras y `summary_extra.json` de 3d y fase 5, y la prueba de colapso con respuesta relativa |
| `make_manuscript_figures.py` | figuras del manuscrito en PDF; acepta nombres (`fig9 figS7`) para regenerar solo algunas |

Las salidas van a `Results/plan/phaseK/`, las de 3d y fase 5 a `Results/plan/extra/` y las figuras del manuscrito a `Results/plan/manuscript/`. La guía completa está en el `README.md` de la raíz y los resultados en `docs/RESULTADOS.md`.
