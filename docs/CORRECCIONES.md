# Errores e incoherencias encontrados y cómo se resolvieron

## En el código original

**1. El experimento 6 sobrecontaba las liberaciones.** `run_experiment6_dose_frequency.py` calculaba el total liberado con `floor(t_max/T)`. Cuando T divide a t_max (T = 5, 10 y 20 con t_max = 100), eso cuenta una liberación en t = t_max que el integrador nunca aplica, porque no hay tiempo posterior en el que actúe (ver `core.release_times`). El costo quedaba inflado en 5, 11 y 25 % para esos T y la costo-eficiencia subestimada en la misma proporción.

- *Corrección:* se usa `core.n_releases(T, t_max) = ceil(t_max/T) − 1`, que además comprueba que coincide con las liberaciones aplicadas. Se regeneraron `experiment6_dose_frequency.csv`, `experiment6_best_dose_frequency.csv` y `experiment6_summary.png`. Los AUC no cambian (diferencia < 5·10⁻¹⁶). Solo cambian `TotalReleased`, `CostEfficiency` y, por tanto, el mejor calendario.
- *Efecto:* el mejor (δ, T) por costo-eficiencia pasa de δ = 10, T = 15 a δ = 10, T = 20 en las 8 combinaciones de topología y estrategia. Por ejemplo, en la trayectoria con liberación dirigida la mejor costo-eficiencia sube de 5,76 a 6,45.

**2. `functions.solve_system_jitter`** también usa `floor(t_max/T)`. Ahí no es un error. Con jitter, la liberación k = t_max/T puede caer antes de t_max y entonces sí se aplica (las que caen fuera de (0, t_max) se filtran). Se dejó como estaba.

**3. Trampa al leer los CSV.** pandas convierte por defecto la cadena `"None"` (estrategia sin liberación) en `NaN`. No es un error del código, pero rompe cualquier análisis que filtre por `Strategy == "None"`. Todos los scripts nuevos leen con `keep_default_na=False` y el README lo advierte.

## En el plan de trabajo

**4. `analysis/make_figures.py` no existe en el repositorio.** El plan dice que ese script calcula el grado del foco de cada instancia, pero no venía en el .zip. Se reemplazó por `plan_common.FOCAL_NETWORKS` y `plan_common.focus_degree`. Las cuatro redes de grado 3 a 6 se eligieron entre las instancias del experimento 3 (mismo generador, misma semilla) tomando la primera semilla de 1 a 10 con el grado buscado: Random semillas 9, 3 y 4 para los grados 3, 4 y 5. **Ninguna instancia Random tiene el foco con grado 6**, así que para el grado 6 se usó ScaleFree semilla 2. Todas son conexas.

**5. El grafo completo no está entre las siete redes, pero el plan lo usa.** Aparece en los pilotos de las fases 1 y 3a, en la prueba 4 de la fase 0 y en la fase 3c (comparación estrella–completo y estabilidad). Las siete redes usan la estrella para el grado 7. Así se resolvió:
- En las fases 1 y 2 (m = 0) la estrella y el completo son dinámicamente idénticos (prueba 4, diferencia < 10⁻¹⁴), así que basta con la estrella. El piloto con el completo es por tanto un piloto con la estrella.
- En la fase 3 se agregó el completo como octava red, porque con m_t > 0 deja de ser equivalente. Son 18 720 corridas en lugar de 16 380. Las tablas permiten filtrar las 7 redes del plan.

**6. Π con τ_d = 1/|r₁| y Π con τ_d = T₀ no son candidatas distintas.** Ni r₁ ni T₀ varían en el diseño de la fase 1, así que las dos versiones difieren en un factor constante, 7,8/5 = 1,56. En escala logarítmica eso es un corrimiento, y cualquier ajuste suave da exactamente el mismo R². La tabla 1.1 las reporta por separado, pero son la misma prueba. Para distinguirlas habría que variar r₁.

**7. T₀ ≈ 7,8 depende de la amplitud.** El período intrínseco del nodo aislado sin liberaciones no es único. El valor linealizado en el equilibrio de coexistencia es 2π/0,928 = 6,77. El equilibrio es un foco débilmente inestable (parte real +0,0008), y el intervalo entre picos de x₁ crece con la amplitud, desde 6,8 al inicio hasta 7,6 a 8,4 en régimen. El valor 7,8 corresponde a una trayectoria concreta. No afecta a los resultados por el punto 6, pero conviene decirlo si se usa T₀ en el texto.

**8. Ponderar por flujo solo en las aristas del foco equivale a la media simple.** La fase 4 propone, si la predicción falla, usar "la media ponderada por el número de parasitoides que cruza cada arista". Desde el foco el flujo hacia cada vecino es el mismo, D_t·x₃,foco, así que ponderar solo esas aristas devuelve la media simple. Se interpretó como una ponderación sobre **todas** las aristas de la red, con peso A_ij·∫x₃ᵢ dt en la corrida dirigida en D_t*.

**9. La regla D_t*·m_t·k_s ≈ |r₃| no tiene la dependencia correcta en r₃.** No es un error del plan, que la propone como hipótesis a contrastar, pero el resultado hay que registrarlo. La pérdida en la cresta, D_t*·m_t·⟨k⟩, es **menor** con |r₃| mayor (0,46 con r₃ = −0,3 y 0,35 con r₃ = −0,45). Por eso la fase 4 usa la constante de la fase 3 con el mismo r₃, porque mezclar los dos valores de r₃ sesgaba la predicción en un factor 2.

**10. "Prominencia" mide solo la caída hacia la derecha.** El plan la define como la caída relativa hasta D_t = 1. Una curva monótona decreciente desde D_t = 10⁻³ tendría su máximo en el borde izquierdo y no contaría como interior, así que la definición es coherente. Aun así, las tablas incluyen también `prom_left`, la caída hasta D_t = 10⁻³, para que una meseta no se confunda con un pico.

**11. El conteo de corridas de la fase 3a incluye las tres estrategias.** El total del plan, 5 × 13 × 3 × 7 × 2 × 2 × 3 = 16 380, cuenta las tres estrategias aunque el diseño solo nombre "dirigida y uniforme". La corrida sin liberación es necesaria para calcular E, así que el conteo es correcto. Se deja anotado porque en la fase 2 la lista de diseño tampoco la menciona.

## Decisiones de implementación que no estaban fijadas en el plan

- **Ajuste de la prueba de colapso.** Polinomio cúbico en log₁₀ de la candidata (cúbico bivariado para el par). Fuera del rango de entrenamiento se prolonga como constante para no premiar ni castigar la extrapolación de un polinomio. El criterio usa el **menor** de los dos R² validados (r₃ fuera y red fuera).
- **Posición del máximo.** Se toma el máximo de la malla y se refina con una parábola en log₁₀ D_t que pasa por él y sus dos vecinos. Solo se refina si el máximo es interior.
- **Tolerancias de los barridos.** Se usó rtol = 10⁻⁸ en lugar de 10⁻¹⁰, validado en la prueba 6 de la fase 0.
- **Malla de 3b.** Se tomaron los 7 valores de D_t de la malla de 3a que caen en [0,03; 1], que es donde está la cresta para todo m_t > 0, con m_d ∈ {0,1; 0,2} y m_t ∈ {0,05; 0,1; 0,2; 0,4}.
- **Fase 4.** Se usó un radio de conexión de 0,5, ℓ ∈ {0,5; 1; 2} (mortalidades medias 0,47, 0,28 y 0,15, en el mismo rango que la fase 3), D_d = 0,03 y D_t en [0,03; 3], porque la cresta puede pasar de 1 cuando m es pequeño.

## Revisión del 17 de octubre de 2026

Se revisó todo el código del plan (integrador con mortalidad, ejecutor, análisis de las fases 1 a 4, fase 3c, fase 4, pruebas y figuras del manuscrito). Se recalcularon desde los CSV todas las cifras de `RESULTADOS.md`. **No se encontraron errores de código** y todas las cifras coinciden. Sí se encontraron dos problemas de interpretación y uno de diseño:

**12. Afirmación demasiado general sobre el AUC con liberación.** El README decía que, con m_t > 0, "el D_t que deja menos plaga es el más bajo". Eso es falso en la trayectoria y el ciclo con CI-B, donde el AUC con liberación dirigida tiene un mínimo interior de hasta 15 % (31 curvas pasan del 5 %). La figura 3.5 muestra solo el grafo completo, que es justo donde la curva es plana. *Corrección:* se cambió el texto del README y de `RESULTADOS.md`, y se añadió la fase 5.

**13. La malla de D_t de la fase 3 era demasiado corta.** Terminaba en D_t = 1, donde quedaba el máximo de 255 de las 384 curvas con m_t > 0. Por eso 13 de las 19 pendientes de la compuerta usaban dos puntos y la constante de la regla estaba sesgada hacia m_t grande y r₃ = −0,45. *Corrección:* se extendió la malla hasta 31,6 sin modificar `phase3a_grid.csv` (sección 3d de `RESULTADOS.md`). El resultado se confirma y queda más sólido (pendiente −0,92 con 4 puntos en 96 condiciones).

**14. La cresta de E_T depende del parasitoide residente.** No es un error, pero condiciona la lectura. Como todos los nodos empiezan con x₃ = 15, el contrafactual sin liberación se deteriora con D_t cuando m_t > 0. *Corrección:* se añadió la fase 5 (x₃(0) = 0), donde el óptimo de D_t es un mínimo real de plaga.

Los resultados nuevos se escribieron en `Model/Python/Results/plan/extra/` con scripts propios (`run_extra_*.py`, `analyze_extra.py`). No se sobrescribió ninguna salida anterior.

## Incidente durante la ejecución

Mientras corría la fase 3b, el archivo `phase3b_grid.csv` se movió de carpeta y luego se devolvió a su lugar. El proceso siguió escribiendo en el mismo archivo (mismo sistema de archivos, mismo descriptor) y el CSV final tiene las 1344 filas esperadas sin duplicados. Se verificó con `wc -l` y con las claves de combinación.
