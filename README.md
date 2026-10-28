# Dispersal costs and network structure in augmentative biological control of *Diaphorina citri*

Code and simulation outputs for the manuscript *Dispersal costs determine when network structure matters for augmentative biological control of* Diaphorina citri (S. Rodríguez-Falcón, H. Arias-Castro, J. Galeano and L. Stucchi).

The model couples four variables in each patch of a network: the psyllid *D. citri* (x₁), host-plant flush (x₂), the parasitoid *Tamarixia radiata* (x₃) and plant vigour (x₄). Parasitoids are released periodically (δ individuals every T time units), either in every patch (uniform) or only in the outbreak patch (targeted). Psyllids and parasitoids disperse along the edges of the graph and may die in transit:

    x1' = f1(x) - Dd (L + md A) x1,      x3' = f3(x) - Dt (L + mt A) x3,

where L is the graph Laplacian, A the adjacency matrix and md, mt are transit-mortality probabilities (scalars, or an n×n matrix m_ij in Experiment 11). Flush and vigour do not disperse.

## Repository layout

```
.
├── README.md
├── requirements.txt
└── code/
    ├── core.py                 local vector field, Jacobian and impulsive integrator
    ├── common.py               parameters (table 1), initial conditions, CSV output
    ├── functions.py            helpers shared with Experiment 5
    ├── network_experiments.py  network model with transit mortality, topologies, metrics
    ├── network_stability.py    master stability function (Benettin exponents)
    ├── plan_common.py          focal graphs, efficiencies and the parallel, resumable runner
    ├── plan_plots.py           palette and style of the diagnostic figures
    ├── test_phase0.py          checks of the transit-mortality code (Methods)
    ├── run_experiment{2..7}_*.py           Experiments 2-7
    ├── run_phase{1..4}_*.py, run_phase3c_stability.py   Experiments 8-11
    ├── analyze_phase{1,2,3}.py             tables and diagnostic figures of Experiments 8-10
    ├── make_manuscript_figures.py          figures 4-8 and S3-S6
    ├── analysis/               figures 2, 3, S1 and S2, and the numbers of Experiments 2-7
    ├── stability/              long-run transverse Lyapunov exponents (table 3) and sigma/(m_t<k>)
    ├── figures/                figures 2, 3, S1 and S2 (PDF)
    └── Results/
        ├── experiment{2..7}_*.csv|png      Experiments 2-7
        ├── plan/phase0/                    checks of the transit-mortality code
        ├── plan/phase{1..4}/               Experiments 8-11: raw grids, tables, diagnostic figures
        ├── plan/manuscript/                figures 4-8 and S3-S6 (PDF)
        └── stability/                      table 3 and the stability numbers of sections 3 and 4
```

The code of Experiments 8-11 keeps the names of the work plan it was written for: phase *k* is Experiment *k* + 7. The diagnostic figures in `Results/plan/phase*/` have Spanish labels and are not part of the manuscript.

## Experiments

Experiments are numbered as in the manuscript (table 2). Experiment 1 was a validation of the Python port of the original Mathematica code and is not reported.

| Experiment | Content | Script | Output in `code/Results/` |
|---|---|---|---|
| 2 | 7 topologies × D_d × D_t × 3 strategies | `run_experiment2_outbreak_containment.py` | `experiment2_*` |
| 3 | 34 stochastic graphs, 10 seeds | `run_experiment3_replicates.py` | `experiment3_*` |
| 4 | stability of the synchronous state | `run_experiment4_stability_validation.py` (t = 1200), `stability/msf_long.py` (t = 30 000, reported values) | `experiment4_*`, `stability/` |
| 5 | single patch, dose and period grid | `run_experiment5_bifurcation_map.py` | `experiment5_*` |
| 6 | dose and frequency | `run_experiment6_dose_frequency.py` | `experiment6_*` |
| 7 | one-at-a-time sensitivity | `run_experiment7_sensitivity.py` | `experiment7_*` |
| 8 | relative dispersal, D_d × D_t | `run_phase1_relative_dispersion.py`, `analyze_phase1.py` | `plan/phase1/` |
| 9 | optimum of D_t without transit cost | `run_phase2_goldilocks.py` (`--check-dd0` for D_d = 0), `analyze_phase2.py` | `plan/phase2/` |
| 10a, 10b | transit mortality of the parasitoid, and of both species | `run_phase3_transit_mortality.py` (`--phase3b` for 10b), `analyze_phase3.py` | `plan/phase3/` |
| 11 | edge-specific transit mortality on random geometric graphs | `run_phase4_edge_mortality.py` | `plan/phase4/` |
| checks | the four checks of the transit-mortality code described in Methods, and three numerical checks | `test_phase0.py` | `plan/phase0/phase0_tests.json` |

`run_phase3c_stability.py` checks numerically the reduction of the synchronous state with transit mortality on regular graphs, and `stability/sigma_ratio.py` computes σ/(m_t⟨k⟩) for the eight graphs of Experiment 10 (section 4).

## Figures and tables

| Item | Script (run from `code/`) | File |
|---|---|---|
| Table 1 | parameters in `common.py` | |
| Figure 1 | drawn in the manuscript (TikZ) | |
| Figures 2, 3, S2 | `analysis/make_figures.py` | `figures/` |
| Figure S1 | `analysis/intrinsic_dynamics.py` (after `make_figures.py`) | `figures/` |
| Table 3, rank correlations of the stability analysis, exponents with transit mortality | `stability/msf_long.py` | `Results/stability/` |
| σ/(m_t⟨k⟩), section 4 | `stability/sigma_ratio.py` | `Results/stability/sigma_ratio.csv` |
| Figures 4-8, S3-S6 | `make_manuscript_figures.py` | `Results/plan/manuscript/` |

`analysis/make_figures.py` also prints the numbers of Experiments 2-7 quoted in the Results. The Spearman correlation it prints for the stability analysis comes from the short run of Experiment 4 (t = 1200); the values reported in the manuscript are those of `stability/msf_long.py`.

## Reproducing the results

The stored results were generated with Python 3.13, numpy 2.5.3, scipy 1.18.1, networkx 3.6.1, pandas 3.0.5 and matplotlib 3.11.2. The random, small-world and scale-free graphs depend on the networkx version, so `requirements.txt` fixes networkx.

```bash
pip install -r requirements.txt
cd code

# Experiments 2-7 (in this order: Experiment 4 reads the summary of Experiment 3)
python run_experiment2_outbreak_containment.py
python run_experiment3_replicates.py
python run_experiment5_bifurcation_map.py
python run_experiment4_stability_validation.py
python run_experiment6_dose_frequency.py
python run_experiment7_sensitivity.py

python test_phase0.py                              # ~4 min, ends with "ALL PHASE-0 TESTS PASS"

python run_phase1_relative_dispersion.py           # Experiment 8, 10 206 runs
python analyze_phase1.py
python run_phase2_goldilocks.py                    # Experiment 9, 6300 runs
python run_phase2_goldilocks.py --check-dd0        # 1092 runs
python analyze_phase2.py
python run_phase3_transit_mortality.py             # Experiment 10a, 18 720 runs
python run_phase3_transit_mortality.py --phase3b   # Experiment 10b, 4032 runs
python run_phase3c_stability.py
python analyze_phase3.py
python run_phase4_edge_mortality.py                # Experiment 11, 900 runs

python stability/msf_long.py --tend 30000          # ~15 min
python stability/sigma_ratio.py

python analysis/make_figures.py
python analysis/intrinsic_dynamics.py
python make_manuscript_figures.py
```

Notes:

- The `run_experiment*.py` scripts run on import. Run them as scripts, never import them.
- The `run_phase*.py` scripts use all CPU cores and are resumable: they append one row per parameter combination, and if interrupted they continue where they stopped. Delete the CSV to rerun an experiment from scratch. `--pilot` runs the pilots of Experiments 8 and 10 (`phase1_pilot.csv`, `phase3a_pilot.csv`).
- On a 2-core machine Experiment 8 takes about 25 min, Experiment 9 about 35 min and Experiment 10a about 40 min.
- Experiments 2-7 use DOP853 with `rtol = 1e-10`, `atol = 1e-12`. The sweeps of Experiments 8-11 use `rtol = 1e-8`, `atol = 1e-10`; check 6 of `test_phase0.py` shows that this changes the efficiencies by about 10⁻⁸ (relative).
- `python stability/msf_long.py --tend 1200 --out Results/stability_short` reproduces the original exponents of Experiment 4 and serves as a regression test.
- Read the CSV files of Experiments 2-7 with `pd.read_csv(path, keep_default_na=False)`; otherwise pandas turns the strategy label `"None"` into `NaN`.

## Output of Experiments 8-11

Each row of `phase1_grid.csv`, `phase2_grid.csv`, `phase3a_grid.csv` and `phase3b_grid.csv` is one parameter combination, solved with the three strategies.

| Column | Meaning |
|---|---|
| `Network`, `k_s` | focal graph and degree of the outbreak patch (node 0) |
| `CI` | initial condition. `A`: 200 psyllids in the outbreak patch and 20 elsewhere (IC-A); `B`: 200 in the outbreak patch and 0 elsewhere (IC-B) |
| `Dd`, `Dt`, `r3`, `md`, `mt`, `tmax` | parameters varied |
| `Delta`, `T`, `n_rel` | dose, period and number of releases applied, n_rel = ⌈t_max/T⌉ − 1 |
| `Q_T`, `Q_U` | parasitoids released, δ·n_rel and 8·δ·n_rel |
| `AUC_None`, `AUC_Uniform`, `AUC_Targeted` | network pest load ∫ Σᵢ x₁ᵢ dt for each strategy |
| `AUCfocus_*`, `ParAUC_*` | pest load in the outbreak patch and parasitoid load of the network |
| `E_T`, `E_U`, `E_ratio` | efficiencies (AUC_None − AUC)/Q and their ratio E_T/E_U |
| `Pi_r1`, `Pi_T0`, `Pi_t` | dimensionless dispersal ratios tested in Experiment 8 |

Focal graphs of Experiments 8-10 (n = 8, outbreak at node 0):

| Name | Generator (seed) | k_s | ⟨k⟩ |
|---|---|---|---|
| Path | path, outbreak at one end | 1 | 1.75 |
| Cycle | cycle | 2 | 2.00 |
| Random_s9 | G(8, 0.4), seed 9 | 3 | 3.00 |
| Random_s3 | G(8, 0.4), seed 3 | 4 | 2.50 |
| Random_s4 | G(8, 0.4), seed 4 | 5 | 3.25 |
| ScaleFree_s2 | Barabási–Albert (m = 2), seed 2 | 6 | 3.00 |
| Star | star, outbreak at the hub | 7 | 1.75 |
| Complete | complete graph (Experiment 10 only) | 7 | 7.00 |

The four stochastic graphs are instances of Experiment 3: for each degree, the first seed between 1 and 10 whose node 0 has that degree.

## Licence and citation

To be added before the repository is made public.
