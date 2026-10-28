# Cambios pendientes en el código y el repositorio

Lista de lo que hay que hacer en los programas y en el repositorio para que el manuscrito sea reproducible y coincida con el texto. Va ordenada por prioridad.

## 1. Estabilidad del estado sincrónico (sección 2(h), tabla 3, sección 3(k))

Las cifras actuales del manuscrito no salen del código del repositorio. Las produce `stability/msf_long.py`, que se entrega con este documento.

- **Añadir al repo** `stability/msf_long.py` y `stability/sigma_ratio.py`, junto con sus salidas:
  - `results/stability/modes.csv`: exponente y error estándar de cada modo μ.
  - `results/stability/table3_stability.csv`: exponente máximo por topología, μ del modo máximo, factor e^{100Λ}, tiempo de e-plegado y sincronía.
  - `results/stability/summary.csv`: ρ de Spearman con su p, rango Monte Carlo al 95 %, la misma correlación con λ₂ y el exponente del modo 1.
  - `results/stability/mortality_regular.csv`: ciclo y grafo completo con m_t ∈ {0; 0,05; 0,1; 0,2; 0,4}.
  - `sigma_ratio.csv`: σ/(m_t⟨k⟩) para los ocho grafos del experimento 10. Respalda la frase de la discusión "entre 0,91 y 1 para m_t ≤ 0,1".
- **Comando:** `python msf_long.py --tend 30000 --out results/stability`. Tarda unos 15 minutos. Con `--tend 1200` reproduce exactamente la tabla original (ρ = −0,52, trayectoria con +2,10·10⁻³), y eso sirve de prueba de regresión.
- **Cambios de método respecto del script original:**
  1. El horizonte pasa de t = 1200 (56 segmentos) a t = 30 000 (1400 segmentos), y se sigue descartando el 30 % inicial.
  2. El criterio de convergencia (diferencia entre mitades < 0,02) se sustituye por el error estándar de 20 medias de bloque.
  3. El exponente máximo se busca en **todos** los modos transversales, no solo en λ₂. En la aleatoria y el mundo pequeño el máximo está en μ ≈ 3,3, y por eso sus valores casi coinciden (−2,185 y −2,193·10⁻³). No es un error de copia.
  4. La correlación con la sincronía se acompaña de un rango que propaga el error de los exponentes por Monte Carlo.
- **Pendiente:** integrar `msf_long.py` en el código propio del repo, o comprobar que el script original con `t_end = 30000` da lo mismo. Así la frase de Métodos sobre la comprobación cruzada se refiere a código del repo.

## 2. Mínimo de la plaga con mortalidad en tránsito (sección 3(i))

- El análisis que calcula la posición del mínimo del AUC con liberación dirigida reporta solo la mediana (0,9 %). Hay que guardar también el máximo y desglosarlo por grafo y por condición inicial.
- El texto ahora dice "hasta 15 % en la trayectoria y el ciclo con CI-B". Esa cifra tiene que salir de una tabla del repo, por ejemplo `phase3/table_pest_minimum.csv`.

## 3. Óptimos en el borde de la malla (sección 3(i))

- El texto dice ahora que en la trayectoria y la estrella el óptimo cae fuera de la malla "for the smaller values of m_t".
- Conviene que la tabla de óptimos del experimento 10a incluya una columna `edge_max` (máximo en D_t = 1) por grafo, m_t, D_d, r₃ y condición inicial, para poder dar la cifra exacta.

## 4. Prominencia en el experimento 11

- El texto define la prominencia respecto del extremo superior de la malla, que es D_t = 3 en el experimento 11. Hay que verificar que el código no tenga `E_T(1)` fijo y use el último punto de la malla.

## 5. Material suplementario

- El texto cita las figuras S3–S6, pero el ESM solo tiene S1 y S2. Hay que generar:
  - **S3:** mapas de E_T por red del experimento 8.
  - **S4:** todas las curvas E_T(D_t) del experimento 9, con ambos horizontes y ambas condiciones iniciales.
  - **S5:** planos (D_t, m_t) de todos los grafos del experimento 10a.
  - **S6:** desplazamiento de D_t* con m_d del experimento 10b, en todos los grafos.
- El título y los autores del ESM salen en rojo porque todavía usan `\todo`. Hay que sustituirlos por los datos reales.

## 6. Carpeta `analysis/`

- Hoy solo existe en Overleaf, y genera fig2, fig3, S1 y S2. Hay que copiarla al repositorio.
- "Data accessibility" promete que el repo incluye los scripts que regeneran todas las figuras. Eso también cubre los de las figuras 4–8 y S3–S6.

## 7. Limpieza del proyecto

- Borrar `conclusions.tex` de la raíz. `main.tex` usa `sections/conclusions.tex`.

## 8. Malla extendida de la fase 3 y fase 5 (fig9, figS7)

- El manuscrito usa la malla original: 127 máximos y pendientes log–log calculadas con pocos puntos por curva.
- Si se adopta la malla extendida, hay que recalcular y actualizar en el texto:
  - el número de máximos reales y su prominencia mediana;
  - la mediana y el rango de ρ_s y ρ_k, y la desviación estándar de log₁₀ρ;
  - la pendiente mediana y su rango intercuartílico;
  - las constantes c y c_k que usa el experimento 11.
- Hay que decidir si la fase 5 entra en el artículo. Si entra, necesita un experimento en la tabla 2, su sección de Métodos y Resultados, y las figuras 9 y S7.

## 9. Datos que faltan en `main.tex`

- El DOI del repositorio, la declaración CRediT, la financiación y los agradecimientos.
