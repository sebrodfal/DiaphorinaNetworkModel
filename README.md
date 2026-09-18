# DiaphorinaNetworkModel

Extensión de red (metapoblacional) de un modelo de control biológico de *Diaphorina citri* (vector del HLB / "greening" de los cítricos) mediante liberaciones periódicas del parasitoide *Tamarixia radiata*.

El modelo de un solo nodo ya está calibrado y en revisión (preprint/artículo base: **ECOMOD-26-1479**). Ese trabajo menciona explícitamente la extensión a "metapoblaciones espaciales / huertos conectados / migración" como trabajo futuro — este repositorio es esa extensión: en vez de un huerto aislado, se modela una red de huertos (nodos) acoplados por dispersión.

## El modelo

Cada nodo corre el mismo sistema de cuatro compartimentos del paper de un nodo:

| Variable | Significado |
|---|---|
| `x1` | *Diaphorina citri* (la plaga/vector) |
| `x2` | Brotes nuevos de la planta |
| `x3` | *Tamarixia radiata* (el parasitoide liberado) |
| `x4` | Vigor de la planta |

Cada `xi` sigue una ecuación logística acoplada a las otras variables locales del mismo nodo (interacción plaga↔parasitoide, plaga↔brotes↔vigor). Cada `T` unidades de tiempo se libera una dosis `δ` de *Tamarixia* (`x3 → x3 + δ`), simulando liberaciones periódicas de control biológico.

**Extensión de red:** los nodos se acoplan mediante la matriz Laplaciana del grafo de la red, con difusión solo en `x1` (Diaphorina, que vuela entre huertos) y `x3` (Tamarixia, que también se dispersa). Brotes (`x2`) y vigor (`x4`) permanecen estrictamente locales a cada nodo — son propiedades de la planta, no se mueven. Las liberaciones de parasitoide pueden aplicarse de forma **uniforme** (todos los nodos) o **focalizada** (solo en un nodo, p. ej. el foco de un brote de plaga).

Con difusión apagada (`Dd = Dt = 0`) cada nodo se reduce exactamente al modelo de un solo nodo del paper en revisión — esto está validado explícitamente (Experimento 1) para que ambos trabajos sean directamente comparables y usen los mismos parámetros calibrados (`Model/Mathematica/tabla_latex.txt`).

## Estructura del repositorio

```
Model/Mathematica/
  Functions.m                          Modelo base (un nodo) + versión de red para uso interactivo/plots,
                                        más las funciones de análisis del paper de un nodo (sweeps de
                                        liberación, comparación de escenarios, chequeo de positividad,
                                        respuesta de Holling II, exponente de invasión, robustez a jitter)
  NetworkExperiments.wl                Infraestructura de red para correr experimentos sin abrir un notebook:
                                        solver headless (SolveNetworkData), generadores de topología
                                        (MakeTopology: Path/Cycle/Star/Complete/Random/SmallWorld/ScaleFree)
                                        y métricas por experimento (AUC, sincronía, CV espacial, etc.)
  NetworkStability.wl                  Análisis de estabilidad de la sincronización de la red (tipo
                                        Master Stability Function): descompone la estabilidad transversal
                                        por autovalor del Laplaciano, con una versión clásica (Floquet,
                                        requiere órbita periódica exacta) y una generalizada (exponente de
                                        Lyapunov de tiempo finito vía renormalización de Benettin, no
                                        requiere periodicidad)
  RunExperiment1_Validation.wls        Exp. 1: la red con Dd=Dt=0 reproduce exactamente el modelo de un nodo
  RunExperiment2_OutbreakContainment.wls  Exp. 2: contención de un brote puntual × 7 topologías × grid de
                                        difusión × 3 estrategias de liberación (ninguna/uniforme/focalizada)
  RunExperiment3_Replicates.wls        Exp. 3: repite el Exp. 2 con varias semillas de grafo para las
                                        topologías estocásticas (media ± desviación, no una sola instancia)
  RunExperiment4_StabilityValidation.wls  Exp. 4: valida y aplica el análisis de estabilidad de sincronización
                                        de NetworkStability.wl a las 7 topologías reales de los Exp. 2/3
  RunExperiment5_BifurcationMap.wls    Exp. 5: mapa de bifurcación (δ,T) del modelo de un nodo a largo plazo —
                                        busca si existe una órbita periódica genuina (prerrequisito del Exp. 4)
  RunExperiment6_DoseFrequency.wls     Exp. 6: barrido dosis-frecuencia (δ,T) por topología, costo-eficiencia
                                        del control dentro del horizonte corto ya validado por el paper base
  RunExperiment7_Sensitivity.wls       Exp. 7: sensibilidad local (±10%) de los 18 parámetros del modelo
                                        sobre la densidad de plaga en la red
  PlotFromCSV_Experiment3.wls          Genera figuras a partir del CSV del Exp. 3
  PlotFromCSV_Experiment6.wls          Genera figuras a partir del CSV del Exp. 6
  tabla_latex.txt                      Tabla de parámetros calibrados, compartida con el paper de un nodo
  Results/                             CSV y figuras de salida de cada experimento

Reports/
  reporte_extension_red.html           Estado del trabajo, hallazgos y decisiones pendientes de autorización
```

## Cómo correr un experimento

Los experimentos corren en modo headless con `wolframscript` (no generan `.nb`), y escriben sus resultados en `Model/Mathematica/Results/` como CSV (y a veces PNG):

```powershell
cd Model/Mathematica
wolframscript -file RunExperiment2_OutbreakContainment.wls
```

Cada script se basta a sí mismo: carga `Functions.m` y/o `NetworkExperiments.wl`/`NetworkStability.wl` según lo que necesite. El orden lógico para reproducir todo desde cero es 1 → 2 → 3 → 5 → 4 → 6 → 7 (el 5 debe correr antes que el 4 porque el 4 depende de lo que el 5 encontró sobre periodicidad).

## Estado actual y hallazgos

El resumen ejecutivo — qué está hecho, qué se encontró, y qué decisiones de alcance faltan por resolver — vive en [Reports/reporte_extension_red.html](Reports/reporte_extension_red.html). Dos hallazgos abiertos, en breve:

- **Aperiodicidad (Exp. 5):** a los parámetros del paper base, el modelo de un nodo no converge a una órbita periódica limpia en ninguna de las 64 combinaciones (δ,T) probadas — las corridas largas son aperiódicas/caóticas. No invalida los Experimentos 2/3/6/7 (corren en el horizonte corto ya validado, `tt=100`), pero sí limita qué tan lejos puede llegar un análisis de estabilidad clásico tipo Floquet.
- **Estabilidad de sincronización (Exp. 4):** la topología de red predice razonablemente bien (y se valida) qué tan sincronizados quedan los nodos, con una excepción sin resolver en la topología "Aleatoria" (ver el reporte para el detalle).
