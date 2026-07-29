# Demo 01: T-Rex+GVS on Hastie's Equicorrelated Blocks (Python)

## Purpose

The following demo presents the T-Rex+GVS selector according to [[1]](#references) and [[2]](#references).
We evaluate grouped-signal recovery under strong equicorrelation, using the classical three-group design of
 Zou & Hastie (2005) [[3]](#references).
 In our demos all variables within an active group are themselves active, so the task is recovering all
 of them, but on an individual-variable basis, rather than using a group-level selection method.
 See [What is actually measured](../README.md#what-is-actually-measured-in-these-demos).
 For details on grouped variable selection and FDR control see [[4]](#references).

This is the Python port of the C++ demo
[demo_trex_gvs_01_mc_sim_hastie_en_blocks](../../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_01_mc_sim_hastie_en_blocks/),
using `trex_selector_neo`'s `TRexGVSSelector`.

---

## Data Generation Parameters (`dgp_hastie_snr`)

We consider the linear model:

```math
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon},
\qquad \boldsymbol{\epsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma_{\varepsilon}^2 \boldsymbol{I}_n)
```

- $\boldsymbol{y} \in \mathbb{R}^n$ is the response vector.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ is the design matrix.
- $\boldsymbol{\beta} \in \mathbb{R}^p$ is the coefficient vector, with $s$ nonzero entries.
- $\boldsymbol{\epsilon}$ is the noise vector, i.i.d. standard normal.
- $\sigma_{\varepsilon}^2$ is the noise variance, calibrated to achieve a target linear signal-to-noise ratio (SNR).
- $n = 200$, $p = 500$, $s = 150$.

The design matrix $\boldsymbol{X}$ is generated from a latent-factor model with three active groups:
```math
X_{ij} = Z_{i,\,g(j)} + \sigma_x\, \xi_{ij}, \qquad \xi_{ij} \sim \mathcal{N}(0,1).
```

- Three latent factors $Z_{\cdot,k} \sim \mathcal{N}(0,1)$ for $k \in \lbrace0,1,2\rbrace$.
- **Group 0** (columns 0–49), **Group 1** (columns 50–99), **Group 2** (columns 100–149) — all active, $\beta_j = 3$.
- **Background** (columns 150–499): i.i.d. $\mathcal{N}(0,1)$, each a singleton group.
- Within-group correlation $\rho = 1/(1+\sigma_x^2)$; the default $\sigma_x = \sqrt{0.01}$ gives $\rho \approx 0.99$.
- SNR calibrated via $\sigma_{\varepsilon}^2 = \mathrm{Var}(\boldsymbol{X}\beta)/\mathrm{SNR}$.

![Design-matrix correlation — Hastie equicorrelated groups, natural order (no clustering)](../dgp_corr_plots/gvs_hastie.png)

Empirical Pearson correlation of an exemplary design drawn from this DGP (`dgp_hastie_snr`, $n=200$, $p=500$,
seed 2026, default $\sigma_x = \sqrt{0.01}$), in **natural variable order — no hierarchical clustering**. The three
dark squares on the diagonal are the active groups (columns 0–49 / 50–99 / 100–149) at $\rho \approx 0.99$; the 350
background columns are mutually uncorrelated. One figure per DGP lives in
[`../dgp_corr_plots/`](../dgp_corr_plots/); regenerate the whole suite with
[`../generate_corr_plots.sh`](../generate_corr_plots.sh) (shared plotter
[`../../corr_plt_utils.py`](../../corr_plt_utils.py)).

---

## Control Parameters

```text
K = 20                       # Random experiments per T-loop iteration
tFDR = 0.1                   # Target FDR
corr_max = 0.98              # HAC auto-clustering correlation threshold
hc_linkage = Single          # Single-linkage HAC
lambda2_method = CV_1SE_CCD  # Elastic-net penalty selection
MC = 200                     # Monte Carlo repetitions per grid point
```

---

## Methods Compared

Three T-Rex+GVS solver variants:

- **TENET** — terminating elastic net pathwise solver (`GVSType.EN`, `ENSolverType.TENET`) derived on the idea of a
  pathwise solver introduced in [[3]](#references).
- **TENET_AUG** — terminating elastic net solver based on row-augmented data to solve via a terminating lasso
   (`GVSType.EN`, `ENSolverType.TENET_AUG`) as proposed in [[1]](#references).
- **TIENET_AUG** — terminating elastic net solver based on row-augmented data to solve via a terminating lasso
   the informed elastic net (`GVSType.IEN`) proposed in [[2]](#references).

---

## Three Parts

1. **2-D SNR $\times$ $\rho$ sweep** — TENET / TENET_AUG / TIENET_AUG over
   $\mathrm{SNR} \in \lbrace0.2, 0.5, 1, 2, 5\rbrace$ and
   $\rho \in \lbrace0.30, 0.50, 0.70, 0.90, 0.95, 0.99\rbrace$ (derived via
   $\sigma_x = \sqrt{(1-\rho)/\rho}$), 200 MC trials per cell.
2. $\lambda_2$**-method comparison** — `CV_1SE_SVD` vs. `CV_1SE_CCD` penalty selection at the
   fixed operating point $\mathrm{SNR}=2$, $\rho=0.7$ (TENET and TIENET_AUG).
3. **HAC linkage comparison** — Single / Complete / Average / WPGMA linkage at the same fixed
   operating point (TENET and TIENET_AUG).

Parts 2–3 are supplementary robustness checks over the $\lambda_2$-selection and linkage implementation choices.

---

## Output Files

Written to `simulation_results/data/`:

- `gvs_Hastie_2d.txt` / `.csv` — 2-D SNR $\times$ $\rho$ grid (Part 1).
- `gvs_Hastie_lambda2_method.txt` / `.csv` — $\lambda_2$-selection comparison (Part 2).
- `gvs_Hastie_hc_linkage.txt` / `.csv` — HAC-linkage comparison (Part 3).

The file names and CSV schemas match the C++ demo, and the suite plotter
[`../trex_gvs_plt_utils.py`](../trex_gvs_plt_utils.py) is the same module the C++ suite uses.

Figures (PNG + PDF) go to `simulation_results/plots/`, produced by `./generate_plots.sh`.

---

## Running the Demo

```bash
python demo_trex_gvs_01_mc_sim_hastie_en_blocks.py       # full run (hours)
python demo_trex_gvs_01_mc_sim_hastie_en_blocks.py 8     # with 8 worker processes
./generate_plots.sh                                      # render the figures below from the saved CSVs
```

Monte Carlo trials are parallelized with `ProcessPoolExecutor`; the optional
first CLI argument sets the number of worker processes (default
`min(6, cpu_count)`). For a quick smoke test, temporarily reduce `num_MC` and
the sweep grids at the top of the script.

---

## Simulation Results

The numbers and figures below are this Python run's own results (200 MC trials per grid point). They agree with
the [C++ reference](../../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_01_mc_sim_hastie_en_blocks/README.md#simulation-results)
up to Monte-Carlo noise — the two ports draw from independent RNG streams (NumPy here, the C++ engine there).

### Part 1 — 2-D SNR $\times$ $\rho$ sweep

- FDR stays at or below the $\mathrm{tFDR} = 0.1$ target in every cell (no cyan violation outlines):
  **TENET** and **TENET_AUG** — numerically identical throughout this demo — realize $0.043\text{–}0.091$, and
  **TIENET_AUG** never exceeds $0.030$.
- **TENET / TENET_AUG** reach $\mathrm{TPR} \approx 1$ once $\mathrm{SNR} \gtrsim 1$ for $\rho \ge 0.5$; the
  weakest column, $\rho = 0.30$, still climbs $0.88 \to 0.96 \to 0.99$ over $\mathrm{SNR} \in \lbrace1, 2, 5\rbrace$.
- **TIENET_AUG** is more conservative (very low FDR, TPR that grows more slowly with SNR — e.g. $0.58$ at
  $\mathrm{SNR}=1$, $\rho = 0.3$ vs. $0.88$ for TENET) and its power **collapses in the extreme $\rho = 0.99$
  column** (TPR $\approx 0.02$, flat in SNR), where the within-group columns become nearly collinear. This
  conservatism is not specific to this design — the same pattern recurs across the other block demos in this
  folder.
- Along the $\rho$ axis, performance is otherwise stable across correlation levels, since all three methods use
  HAC-based grouping to aggregate evidence across strongly correlated columns; only $\rho = 0.99$ stresses
  TIENET_AUG.

TPR (top) and FDR (bottom) heatmaps, one column per solver; FDR cells above the
$\mathrm{tFDR} = 0.1$ target would be outlined in cyan.

![Hastie SNR × rho sweep](simulation_results/plots/gvs_Hastie_2d.png)

---

### Part 2 — $\lambda_2$ selection method

- At the fixed operating point both backends behave the same where it matters: TENET realizes
  FDR $0.071$ (`CV_1SE_SVD`) vs. $0.070$ (`CV_1SE_CCD`) at TPR $1.000$, and TIENET_AUG reaches TPR
  $1.000$ vs. $0.996$.
- The only visible difference is TIENET_AUG's FDR ($0.004$ SVD vs. $0.017$ CCD) — both far below the
  $0.1$ target, so the pipeline is insensitive to the penalty-selection backend.

Grouped bar charts at the fixed operating point $\mathrm{SNR} = 2$, $\rho = 0.7$; error bars show
$\pm 1$ Monte Carlo standard deviation.

![lambda2 selection comparison](simulation_results/plots/gvs_Hastie_lambda2_method.png)

---

### Part 3 — HAC linkage

- Single / Complete / Average / WPGMA linkage produce **numerically identical** results here (TENET:
  FDR $0.070$, TPR $1.000$; TIENET_AUG: FDR $0.017$, TPR $0.996$) — with groups this strongly
  equicorrelated, every linkage rule cuts the dendrogram the same way.
- Taken together, Parts 2 and 3 confirm that the Part 1 conclusions are not sensitive to the
  $\lambda_2$-selection backend or the linkage rule for this data scenario.

Grouped bar charts at the same fixed operating point as Part 2; error bars show $\pm 1$ Monte
Carlo standard deviation.

![HAC linkage comparison](simulation_results/plots/gvs_Hastie_hc_linkage.png)

---

## References

1. Machkour, J., Muma, M., & Palomar, D. P., "False Discovery Rate Control for Grouped Variable Selection
   in High-Dimensional Linear Models using the T-Knock Filter.", European Signal Processing Conference (EUSIPCO), 2022,
    pp. 892–896, EURASIP.
    [DOI-Link](https://doi.org/10.23919/EUSIPCO55093.2022.9909883)
2. Machkour, J., Muma, M., & Palomar, D. P., "The Informed Elastic Net for Fast Grouped Variable Selection and
   FDR Control in Genomics Research.", Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP),
    2023, pp. 466–470, IEEE.
    [DOI-Link](https://doi.org/10.1109/CAMSAP58249.2023.10403489)
3. Zou, H., & Hastie, T. (2005). "Regularization and variable selection via the elastic net." *Journal of the Royal
   Statistical Society: Series B (Statistical Methodology)*, 67(2), pp. 301–320.
   [DOI-Link](https://doi.org/10.1111/j.1467-9868.2005.00503.x)
4. Dai, R., & Barber, R. F., "The knockoff filter for FDR control in group-sparse and multitask regression."
   Proceedings of the 33rd International Conference on Machine Learning (ICML), vol. 48, pp. 1851–1859, 2016.
   [DOI-Link](http://proceedings.mlr.press/v48/daia16.pdf)
