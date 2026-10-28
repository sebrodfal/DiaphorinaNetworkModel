# Extras

Todo lo que no hace falta para reproducir el manuscrito. El código y los resultados del artículo están en `code/` y se describen en el `README.md` de la raíz.

| Carpeta | Contenido |
|---|---|
| `mathematica/` | Versión original del modelo en Mathematica, de la que se tradujo `code/`. `tabla_latex.txt` es la fuente de los parámetros de `common.py`. |
| `python_original/` | Partes de la traducción a Python que el artículo no usa: el experimento 1 (validación de la traducción), los gráficos `plot_from_csv_experiment{3,6}.py`, el cuaderno `Model.ipynb` y el README técnico original (`README_model_python.md`). |
| `extended_analyses/` | Análisis añadidos después de la revisión del 17 de octubre de 2026 que todavía no están en el manuscrito: la fase 3d (malla de D_t extendida hasta 31,6) y la fase 5 (liberación sin parasitoide residente), con sus CSV, tablas y las figuras 9 y S7. |
| `docs/` | Plan de trabajo (`plan_de_trabajo.pdf`), resultados fase por fase (`RESULTADOS.md`), errores encontrados y corregidos (`CORRECCIONES.md`), el README de trabajo anterior a esta reorganización (`README_trabajo.md`) y el informe de septiembre de 2026 (`historial/`). |
| `manuscrito/` | Copia del proyecto de Overleaf al 28 de octubre de 2026 (`overleaf_copia_2026-10-28/`) y las secciones corregidas por Luciano ese mismo día, con su lista de cambios pendientes (`cambios_luciano_2026-10-28/`). La versión vigente del manuscrito está en Overleaf. |

## Cómo correr los análisis extendidos

Los scripts de `extended_analyses/` importan el modelo desde `code/` y leen y escriben en `extended_analyses/Results/`.

```bash
cd extras/extended_analyses
python run_extra_phase3a_extended.py   # fase 3d, 8640 corridas (38 min en 14 núcleos)
python run_extra_no_resident.py        # fase 5, 18 000 corridas (20 min en 14 núcleos)
python analyze_extra.py                # tablas, figA, figB y summary_extra.json
python make_extra_figures.py           # figures/fig9 y figures/figS7
```

Los scripts de `python_original/` se dejan como estaban. Para correrlos hay que poner `code/` en el `PYTHONPATH`, y escriben en `code/Results/`.

Las rutas que aparecen en el texto de `docs/README_trabajo.md` y `docs/CORRECCIONES.md` son las de la estructura anterior: `Model/Python/` es ahora `code/`, y `Model/Python/Results/plan/extra/` es ahora `extras/extended_analyses/Results/`.
