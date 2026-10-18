# Dispersión relativa, hipótesis Goldilocks y mortalidad en tránsito

Modelo en red de *Diaphorina citri* (psílido, x₁), brotes (x₂), *Tamarixia radiata* (parasitoide, x₃) y vigor de planta (x₄), con liberaciones periódicas del parasitoide. El repositorio contiene:

- la versión original en Mathematica (`Model/Mathematica`),
- la versión en Python del modelo (`Model/Python`, experimentos 1 a 7 del manuscrito),
- la implementación completa del plan de trabajo [`docs/plan_de_trabajo.pdf`](docs/plan_de_trabajo.pdf), fases 0 a 4,
- dos análisis añadidos después de la revisión del 17 de octubre de 2026: la fase 3 con una malla de D_t extendida (3d) y la fase 5, liberación sin parasitoide residente.

| Documento | Contenido |
|---|---|
| Este README | Qué se hizo, resultados principales, cómo reproducir y cómo usar los resultados |
| [`docs/RESULTADOS.md`](docs/RESULTADOS.md) | Explicación fase por fase con todas las figuras y tablas |
| [`docs/CORRECCIONES.md`](docs/CORRECCIONES.md) | Errores del código e incoherencias del plan encontrados y cómo se resolvieron |
| [`Model/Python/README.md`](Model/Python/README.md) | Referencia técnica del código |

## Resultados en una página

**Fase 0 (aprobada).** El integrador admite mortalidad en tránsito m_d y m_t, dos condiciones iniciales (CI-A y CI-B) y cuenta bien el número de liberaciones. Con m = 0 reproduce las 336 corridas del experimento 2 con un error relativo máximo de 4·10⁻¹⁶, y pasa la prueba de balance de masa, la de D = 0 y la de equivalencia entre la estrella y el grafo completo.

**Fase 1 (resultado negativo).** Ninguna variable adimensional colapsa la eficiencia de la liberación dirigida E_T. El mejor R² con validación cruzada es 0,41 en CI-A y 0,18 en CI-B, cuando el criterio pedía 0,9. La razón es que |r₃| no actúa solo como una escala de tiempo. Una mortalidad natural mayor del parasitoide duplica el AUC de plaga sin liberación y por eso **aumenta** E_T, mientras que Π = D_tτ_t/(D_dτ_d) mueve r₃ en sentido contrario.

**Fase 2 (sin Goldilocks sin costo).** Sin mortalidad en tránsito, E_T crece con D_t hasta el borde de la malla (D_t = 1) en las 42 curvas de CI-A. En CI-B aparecen nueve mesetas con un máximo interior, pero su prominencia no pasa de 1,9 % y ninguna se repite en los dos horizontes. La caída de 5,30 a 4,91 que se veía en el experimento 2 con D_d = 0 es real. Es un valle cuyo mínimo cae entre D_t ≈ 0,015 y 0,1 según la red. Más allá, la eficiencia supera de nuevo 5,30, desde D_t ≈ 0,06 en la estrella y desde D_t ≈ 0,5 en la trayectoria.

**Fase 3 (hay cresta, pero hay que interpretarla bien).** Con m_t > 0 aparece un D_t óptimo muy marcado, con prominencia mediana del 73 %, que baja como D_t* ∝ m_t^−0,91. La compuerta formal hacia la fase 4 se cumple. Hay que tener en cuenta dos cosas al leerla.
1. El pico de E_T **no** es un mínimo del AUC de plaga con liberación. Aparece porque en la corrida sin liberación el parasitoide residente colapsa (su AUC cae a entre 2 y 40 % del valor que tiene sin mortalidad, con mediana del 20 %) mientras que las liberaciones lo mantienen. Es un óptimo de eficiencia relativa a un contrafactual que se deteriora.
2. La regla D_t*·m_t·k_s ≈ |r₃| se cumple solo en orden de magnitud. El cociente tiene una mediana de 1,33 en CI-B, pero va de 0,49 a 5,66 y depende de la red y de r₃. Con el **grado medio** ⟨k⟩ en lugar del grado del foco la dispersión se reduce a la mitad (de 0,60 a 2,05), y el valor absoluto de la pérdida en la cresta, D_t*·m_t·⟨k⟩, se parece más a b₃₁·x̄₁ − |r₃| que a |r₃|. Aquí x̄₁ es la densidad de plaga sin parasitoide.

**Fase 3d (malla extendida, añadida).** La malla original de D_t terminaba en 1 y ahí quedaban 255 de las 384 curvas con m_t > 0. Con la malla extendida hasta 31,6, **las 384 tienen un máximo interior** y D_t* ∝ m_t^−0,92 (4 puntos por condición, 96 condiciones). Con ⟨k⟩, la pérdida en la cresta D_t*·m_t·⟨k⟩ es casi constante dentro de cada r₃ (desviación de 0,04 a 0,06 en log₁₀): vale 0,475 con r₃ = −0,3 y 0,325 con r₃ = −0,45 en CI-B.

Otros resultados de la fase 3. Con m_t > 0 la estrella y el grafo completo dejan de ser equivalentes, con hasta 42 % menos AUC en la estrella. La correlación entre el grado del foco y el AUC con liberación dirigida (−0,99 en CI-A con m = 0) pasa a +0,18 con m_t = 0,4. El estado sincrónico del ciclo y del grafo completo sigue siendo estable para todos los m_t. La mortalidad del psílido m_d baja D_t*. La mediana del desplazamiento es de 2 a 4 % con D_d = 0,01, de 3 a 11 % con D_d = 0,03 y de 18 a 45 % con D_d = 0,1, para m_d = 0,1 y 0,2.

**Fase 4 (opcional, se hizo porque la compuerta pasó).** En 30 grafos geométricos con mortalidad que depende de la distancia, todas las curvas tienen un máximo interior. La regla de la fase 3 con k_s y la mortalidad media de las aristas del foco predice D_t* con un error mediano del 22 % (correlación 0,86). La versión con ⟨k⟩ y la mortalidad media de toda la red lo predice con un error del 10 % (correlación 0,99).

**Fase 5 (añadida, el Goldilocks de control).** Si el parasitoide se libera donde todavía no está establecido (x₃(0) = 0), la corrida sin liberación no depende de D_t y el óptimo de D_t es un **mínimo real de plaga**. Aparece en 236 de 240 curvas y supera en más de 5 % a ambos extremos en 110. El D_t óptimo baja con m_t (pendiente −0,73) y es asimétrico: dispersarse demasiado deja entre 44 y 47 % más plaga, y dispersarse poco entre 3 y 7 % (hasta 18 % en la trayectoria). Este óptimo de control está a un D_t unas 30 veces menor que la cresta de eficiencia de la fase 3.

**Corrección al experimento 6.** Contaba con `floor(t_max/T)` una liberación en t = t_max que nunca se aplica. Tras corregirlo, el calendario de mejor costo-eficiencia pasa de T = 15 a T = 20 en las ocho combinaciones de topología y estrategia.

## Estructura

```
.
├── README.md                      este archivo
├── docs/
│   ├── plan_de_trabajo.pdf        plan original
│   ├── RESULTADOS.md              resultados fase por fase (incluye 3d y fase 5)
│   ├── CORRECCIONES.md            errores e incoherencias resueltos (incluye la revisión del 17-oct)
│   └── historial/                 informe de estado de septiembre de 2026 (antes del plan)
└── Model/
    ├── Mathematica/               versión original (Functions.m, NetworkExperiments.wl, ...)
    └── Python/
        ├── core.py, common.py, functions.py, network_experiments.py, network_stability.py
        ├── run_experiment{1..7}_*.py, plot_from_csv_experiment{3,6}.py, Model.ipynb
        ├── plan_common.py             definiciones del plan (redes, Π, Π_t, eficiencias, ejecutor)
        ├── plan_plots.py              estilo y paleta de las figuras del plan
        ├── test_phase0.py             pruebas de la fase 0 (7 pruebas)
        ├── run_phase1_relative_dispersion.py   + analyze_phase1.py
        ├── run_phase2_goldilocks.py            + analyze_phase2.py
        ├── run_phase3_transit_mortality.py     + analyze_phase3.py
        ├── run_phase3c_stability.py
        ├── run_phase4_edge_mortality.py
        ├── run_extra_phase3a_extended.py       fase 3d  ┐
        ├── run_extra_no_resident.py            fase 5   ├ + analyze_extra.py
        ├── make_manuscript_figures.py          figuras del manuscrito (PDF, en inglés)
        └── Results/
            ├── experiment*.csv|png    experimentos 1 a 7 (el 6 regenerado con la corrección)
            └── plan/
                ├── phase{0..4}/       CSV crudos, tablas y figuras del plan
                ├── extra/             fase 3d y fase 5 (CSV, tablas, figuras, summary_extra.json)
                └── manuscript/        fig4 a fig9 y figS3 a figS7
```

## Cómo reproducir

Requiere Python ≥ 3.10. Los resultados se generaron con Python 3.13, numpy 2.5.3, scipy 1.18.1, networkx 3.6.1, pandas 3.0.5 y matplotlib 3.11.2. Las redes aleatorias dependen de la versión de networkx. Con otra versión, la prueba 1 de la fase 0 (reproducir el experimento 2) puede fallar en las topologías Random, SmallWorld y ScaleFree aunque el código esté bien. Si pasa, conviene fijar `networkx==3.6.1`.

```bash
cd Model/Python
pip install -r requirements.txt

python test_phase0.py                              # ~4 min, debe terminar en "ALL PHASE-0 TESTS PASS"

python run_phase1_relative_dispersion.py --pilot   # 243 corridas
python run_phase1_relative_dispersion.py           # 10 206 corridas
python analyze_phase1.py

python run_phase2_goldilocks.py                    # 6300 corridas
python run_phase2_goldilocks.py --check-dd0        # 1092 corridas
python analyze_phase2.py

python run_phase3_transit_mortality.py --pilot     # 585 corridas
python run_phase3_transit_mortality.py             # 18 720 corridas
python run_phase3_transit_mortality.py --phase3b   # 4032 corridas
python run_phase3c_stability.py
python analyze_phase3.py

python run_phase4_edge_mortality.py                # 900 corridas, lee la compuerta de la fase 3

python run_extra_phase3a_extended.py               # fase 3d, 8640 corridas
python run_extra_no_resident.py                    # fase 5, 18 000 corridas
python analyze_extra.py

python make_manuscript_figures.py                  # todas las figuras; o p. ej. "fig9 figS7"
```

Los tiempos son de una máquina de 2 núcleos: fase 1 unos 25 min, fase 2 unos 35 min, fase 3a unos 40 min y fase 3b unos 9 min. En una de 14 núcleos, la fase 3d tardó 38 min (los D_t grandes son lentos) y la fase 5, 20 min. Los scripts `run_phase*` usan todos los núcleos y son **reanudables**. Escriben una fila por combinación a medida que terminan, y si se interrumpen basta con volver a ejecutarlos para que sigan donde quedaron. Para rehacer una fase desde cero hay que borrar su CSV.

Los barridos usan `rtol = 1e-8` y `atol = 1e-10`. La prueba 6 de la fase 0 comprueba que esto cambia E_T y E_U en menos de 10⁻⁸ en relación con las tolerancias del manuscrito (`1e-10` y `1e-12`), que se siguen usando en la prueba de reproducción y en los experimentos 1 a 7.

## Cómo usar los resultados

### Archivos crudos

Cada fila de `phase1_grid.csv`, `phase2_grid.csv`, `phase3a_grid.csv` y `phase3b_grid.csv` es **una combinación de parámetros**, en la que se resolvieron las tres estrategias con parámetros idénticos. Las columnas son:

| Columna | Significado |
|---|---|
| `Network`, `k_s` | red (ver tabla de redes) y grado del foco (nodo 0) |
| `CI` | `A`: 200 psílidos en el foco y 20 en los demás nodos; `B`: 200 en el foco y 0 en los demás |
| `Dd`, `Dt`, `r3`, `md`, `mt`, `tmax` | parámetros variados |
| `Delta`, `T`, `n_rel` | dosis, período y liberaciones aplicadas, n_rel = ⌈t_max/T⌉ − 1 |
| `Q_T`, `Q_U` | parasitoides liberados, δ·n_rel·1 y δ·n_rel·8 |
| `AUC_None`, `AUC_Uniform`, `AUC_Targeted` | ∫ Σᵢ x₁ᵢ dt en cada estrategia |
| `AUCfocus_*`, `ParAUC_*` | AUC de plaga en el foco y AUC total del parasitoide |
| `E_T`, `E_U`, `E_ratio` | (AUC_None − AUC_estrategia)/Q y el cociente E_T/E_U |
| `Pi_r1`, `Pi_T0`, `Pi_t` | Π con τ_d = 1/\|r₁\| = 5 o τ_d = T₀ = 7,8, y Π_t = D_t k_s/(\|r₃\| + m_t D_t k_s) |

Los CSV se leen con `pd.read_csv(ruta, keep_default_na=False)`. Sin esa opción, pandas convierte la estrategia `"None"` de los CSV de los experimentos 2 a 7 en `NaN`.

### Redes focales

| Nombre | Generador (semilla) | k_s | ⟨k⟩ |
|---|---|---|---|
| Path | trayectoria, foco en un extremo | 1 | 1,75 |
| Cycle | ciclo | 2 | 2,00 |
| Random_s9 | G(8; 0,4), semilla 9 | 3 | 3,00 |
| Random_s3 | G(8; 0,4), semilla 3 | 4 | 2,50 |
| Random_s4 | G(8; 0,4), semilla 4 | 5 | 3,25 |
| ScaleFree_s2 | Barabási–Albert (m = 2), semilla 2 | 6 | 3,00 |
| Star | estrella, foco en el centro | 7 | 1,75 |
| Complete | grafo completo (solo fases 3 y piloto) | 7 | 7,00 |

Las cuatro redes aleatorias son instancias del experimento 3, con el mismo generador y la misma semilla. En cada caso se tomó la primera semilla entre 1 y 10 cuyo nodo 0 tiene el grado buscado.

### Qué se puede afirmar y qué no

Afirmaciones respaldadas:

- La eficiencia de la liberación dirigida no depende de una sola combinación adimensional de D_t, D_d y r₃ (tabla 1.1). Ni D_t/D_d ni Π sirven para generalizar a otros pares plaga–enemigo.
- Sin costo de dispersión, en este modelo y con estos horizontes, un parasitoide más móvil nunca es peor para la liberación dirigida, salvo mesetas de menos del 2 %.
- Con mortalidad en tránsito existe un D_t que maximiza E_T y escala como 1/(m_t·⟨k⟩). Con la malla extendida la cresta es interior en las 384 condiciones con m_t > 0, con pendiente −0,92. La compuerta de la fase 3 se cumple y la fase 4 confirma la forma de la regla en redes heterogéneas.
- Cuando el parasitoide no está establecido antes de liberar (fase 5), existe un D_t que minimiza la plaga y baja con m_t. Un parasitoide demasiado móvil es mucho peor (entre 44 y 47 % más plaga) que uno poco móvil (entre 3 y 7 %, hasta 18 % en la trayectoria).
- La proposición de equivalencia estrella–completo y la correlación negativa entre grado del foco y AUC solo valen con m = 0.

Hay que reportarlas con salvedades:

- La cresta de E_T mide cuánto evita la liberación **respecto de no liberar**, con un parasitoide residente que colapsa en esa corrida de referencia. En la cresta, el AUC de plaga con liberación es mucho mayor que con D_t bajo (figura 3.5). Por eso no es un óptimo de control. Con residente, el D_t que deja menos plaga es casi el más bajo en las redes muy conectadas, pero no en la trayectoria ni en el ciclo con CI-B, donde hay un mínimo interior de hasta 15 %. El óptimo de control propiamente dicho es el de la fase 5.
- La constante de la regla depende de r₃. Con la malla extendida, D_t*·m_t·⟨k⟩ vale 0,475 con r₃ = −0,3 y 0,325 con r₃ = −0,45 (CI-B). La fórmula |r₃| del plan no recoge esa dependencia y la explicación propuesta con b₃₁x̄₁ − |r₃| es una lectura posterior a los datos, no una derivación.
- El óptimo de control de la fase 5 tiene un beneficio moderado frente a un D_t bajo en las redes muy conectadas (1 a 7 % en el grafo completo y la estrella). Su pendiente con m_t (−0,73) no es exactamente −1.
- Las curvas con t_max = 300 tienen oscilaciones en D_t (figura 2.1) que vienen de la dinámica oscilatoria del modelo, no de error numérico.
- La correlación −0,98 entre grado del foco y AUC se obtiene en CI-A, con D_d = D_t ≈ 0,03 y m = 0. Con D = 0,01 o 0,1 baja a −0,42 y −0,63, y en CI-B es positiva (entre +0,40 y +0,88).

### Para el manuscrito

Los tres resultados que el plan señalaba como decisivos son:

| Entregable | Archivo | Resultado |
|---|---|---|
| Tabla 1.1, prueba de colapso | `Results/plan/phase1/table1_1_cv_r2.csv` | negativo, ninguna candidata colapsa |
| Figura 3.1, cresta en (D_t, m_t) | `Results/plan/phase3/fig3_1_maps_Dt_mt.png` | positivo, la cresta existe y baja con m_t |
| Figura 3.2, regla D_t*·m_t·k_s ≈ \|r₃\| | `Results/plan/phase3/fig3_2_rule.png` | parcial, vale en orden de magnitud y mejora con ⟨k⟩ |

Añadidos después de la revisión:

| Entregable | Archivo | Resultado |
|---|---|---|
| Regla con la malla extendida (3d) | `Results/plan/extra/figA_ext_rule.png`, `manuscript/figS7_crest_rule_extended_grid.pdf` | positivo: D_t* ∝ m_t^−0,92 en las 96 condiciones; constante con ⟨k⟩ muy estable dentro de cada r₃ |
| Óptimo de control sin residente (fase 5) | `Results/plan/extra/figB_no_resident.png`, `manuscript/fig9_control_optimum_no_resident.pdf` | positivo: Goldilocks de control, asimétrico, que baja con m_t |

Antes de enviar el manuscrito hay que decidir dos cosas: si las figuras 3.1 y 3.2 (y fig7) se reemplazan por la versión con la malla extendida, y si la fase 5 entra como resultado principal o como complemento.

## Licencia y cita

El repositorio no trae licencia. Hay que añadir una antes de hacerlo público.
