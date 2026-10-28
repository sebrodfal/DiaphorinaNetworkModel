# Overleaf project

- `main.tex` compiles with pdfLaTeX + BibTeX (Overleaf default). `esm.tex` is the electronic supplementary material and compiles separately.
- Filled sections: `sections/methods.tex` and `sections/results.tex`. Red text marks placeholders.
- Figures are regenerated from the model outputs with

      python analysis/make_figures.py  <path>/Model/Python/Results
      python analysis/intrinsic_dynamics.py  <path>/Model/Python

  The first script also prints every number quoted in the Results.
