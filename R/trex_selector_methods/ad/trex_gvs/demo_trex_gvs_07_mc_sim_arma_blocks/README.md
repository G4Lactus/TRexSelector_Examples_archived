# Demo 07: T-Rex+GVS on Heterogeneous ARMA Blocks

## Purpose

Stress-test the Terminating Random Experiments (T-Rex) selector with grouped variable selection (GVS), T-Rex+GVS,
when different blocks follow **different time-series structures** rather than a single uniform correlation model.
This is the most heterogeneous of the four-block scenarios. Each block has its own autoregressive moving-average
(ARMA) specification, which tests whether the selector can still recover grouped signal when "the group correlation
model" is not one fixed thing. Selection and evaluation are per variable (see
[What is actually measured](../../../../../cpp/trex_selector_methods/ad/trex_gvs/README.md#what-is-actually-measured-in-these-demos)
in the C++ suite README).

This is the R port of the C++ demo
[demo_trex_gvs_07_mc_sim_arma_blocks](../../../../../cpp/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_07_mc_sim_arma_blocks/README.md),
using the `TRexGVSSelector` class of `TRexSelectorNeo`. It keeps the data-generating process (DGP), control
parameters, sweep, console summaries and result-file naming of the C++ demo. Indices are 1-based (R convention).

---

## Data Generation Parameters (`dgp_arma_blocks_snr`)

The data follow the linear model

```math
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon},
\qquad \boldsymbol{\epsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma_{\varepsilon}^2 \boldsymbol{I}_n)
```

- $\boldsymbol{y} \in \mathbb{R}^n$ is the response vector.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ is the design matrix.
- $\boldsymbol{\beta} \in \mathbb{R}^p$ is the coefficient vector, with $s$ nonzero entries.
- $\boldsymbol{\epsilon}$ is the noise vector, independent and identically distributed (i.i.d.) normal.
- $\sigma_{\varepsilon}^2$ is the noise variance, calibrated to a target signal-to-noise ratio (SNR) as
  $\sigma_{\varepsilon}^2 = \widehat{\mathrm{Var}}(\boldsymbol{X}\boldsymbol{\beta}) / \mathrm{SNR}$.
- $n = 200$, $p = 500$, $s = 150$ (the sum of the three active blocks).

The design matrix keeps the four-block geometry of Demos 03, 05 and 06 (sizes $\lbrace20, 50, 80, 65\rbrace$, 215
grouped columns in total), but **each block follows a distinct ARMA process**. In every row, the columns of a block
are one realization of its process, so for a stationary process the correlation between two columns of a block is
the autocorrelation at their column distance. The helper `.sim_arma()` generates each realization by direct
recursion after a burn-in of 100 steps, as the C++ generator does (it does not call `arima.sim` and does not enforce
stationarity), and every block column is then standardized to unit empirical variance.

- **Block 1** (20 vars, active, $\beta_j = 3$): first-order autoregressive (AR(1)), $\phi_1 = \texttt{ar\_coef}$.
- **Block 2** (50 vars, active): moving-average (MA(3)), $\theta = (0.5, 0.3, 0.1)$.
- **Block 3** (80 vars, active): ARMA(2,1), $\phi = (\texttt{ar\_coef}, \texttt{ar\_coef}/5)$, $\theta = 0.5$.
- **Block 4** (65 vars, inactive trap): AR(1), $\phi_1 = \texttt{ar\_coef}$, structurally identical to Block 1.
- Each trial shuffles the block order and spreads the 285 remaining columns as gaps around the blocks, with at least
  one gap column between two blocks. The gap columns are i.i.d. $\mathcal{N}(0,1)$.

Unlike Demo 06, the sweep parameter here is `ar_coef` (shared by the AR(1) and ARMA(2,1) blocks), **not** $\rho$. The
MA(3) block has no AR component at all, so it does not respond to the sweep. The two AR coefficients of the ARMA(2,1)
block sum to $1.2 \cdot \texttt{ar\_coef}$, so that block is stationary only for `ar_coef` $< 5/6 \approx 0.833$ (see
the comment on `.sim_arma()`). At `ar_coef` $= 0.9$ and $0.95$ its recursion grows geometrically, and after the burn-in
the columns of the block are close to multiples of one another, that is, nearly collinear.

![Design-matrix correlation, block-diagonal ARMA, natural order (no clustering)](../dgp_corr_plots/gvs_arma_blocks.png)

Empirical Pearson correlation of an exemplary design drawn from this DGP (`dgp_arma_blocks_snr`, $n=200$, $p=500$,
`ar_coef` $=0.8$, seed 2026), in **natural variable order without hierarchical clustering**. The heterogeneity is
visible directly. The ARMA(2,1) block (80 vars) carries a broad, slowly decaying band, the two AR(1) blocks (20 and 65
vars, active and trap) look alike, and the MA(3) block (50 vars) is nearly a bare diagonal, since its correlation
vanishes beyond lag 3. One figure per DGP lives in [`../dgp_corr_plots/`](../dgp_corr_plots/). Regenerate the whole
suite with [`../generate_corr_plots.sh`](../generate_corr_plots.sh) (shared plotter
[`../../../corr_plt_utils.R`](../../../corr_plt_utils.R)).

---

## Control Parameters

```text
K = 20                          # Random experiments per T-loop iteration
tFDR = 0.1                      # Target FDR
corr_max = 0.98                 # HAC auto-clustering correlation threshold
hc_linkage = "Single"           # Single-linkage HAC
lambda2_method = "CV_1SE_CCD"   # Elastic-net penalty selection
num_MC = 200                    # Monte Carlo repetitions per grid point
seed = 2026                     # Base seed, trial k uses 2026 + k - 1
```

---

## Methods Compared

Three T-Rex+GVS solver variants are compared. **TENET** runs the elastic net (EN) [[3]](#references),
**TENET_AUG** its row-augmented form [[1]](#references), and **TIENET_AUG** the informed elastic net (IEN)
[[2]](#references). All three use single-linkage hierarchical agglomerative clustering (HAC) grouping and
`CV_1SE_CCD` for the $\lambda_2$ selection, as in the C++ and Python demos. The script lists them in `GVS_METHODS`.

- **TENET**: `GVS_type = "EN"`, `en_solver = "TENET"`
- **TENET_AUG**: `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **TIENET_AUG**: `GVS_type = "IEN"`, with `en_solver` left at the library default, which the IEN path ignores

Each Monte Carlo (MC) trial builds the selector through the shared adapter `.gvs_select()` in
[`../../../simulation_utils.R`](../../../simulation_utils.R), here for TENET:

```r
library(TRexSelectorNeo)

ctrl <- trex_control(
  solver = "TLARS", K = 20L, max_dummy_multiplier = 10, use_max_T_stop = TRUE,
  lloop_strategy = "HCONCAT", tloop_stagnation_stop = TRUE,
  tloop_max_stagnant_steps = 5, use_memory_mapping = FALSE
)
sel <- TRexGVSSelector$new(
  X, y, tFDR = 0.1, seed = trial_seed, verbose = FALSE,
  gvs_control = trex_gvs_control(
    gvs_type = "EN", corr_max = 0.98, hc_linkage = "Single", en_solver = "TENET",
    lambda2_method = "CV_1SE_CCD", cv_seed = trial_seed
  ),
  control = ctrl
)
sel$select()
fdp <- trex_fdp(sel$selected_indices, true_support)
tpp <- trex_tpp(sel$selected_indices, true_support)
```

TENET_AUG changes only `en_solver`. TIENET_AUG sets `gvs_type = "IEN"` and the base solver `"TIENET_AUG"`, so its
rows run the augmented IEN solver as in C++ and Python. The trial seed $2026 + k - 1$ of trial $k$ seeds the DGP, the
selector and the cross-validation folds.

---

## The Sweep

A single **2-D SNR $\times$ `ar_coef` grid**, over
$\mathrm{SNR} \in \lbrace0.2, 0.5, 1, 2, 5\rbrace$ and
`ar_coef` $\in \lbrace0.30, 0.50, 0.70, 0.80, 0.90, 0.95\rbrace$
(the AR coefficient shared by the AR(1) and ARMA(2,1) blocks, labeled `ac=` in the result tables), with 200 MC trials
per cell. Because the MA(3) block carries no AR component, only three of the four blocks respond to the `ar_coef`
axis.

---

## Output Files

The run writes `gvs_arma_blocks_2d.txt` and `gvs_arma_blocks_2d.csv` to `simulation_results/data/`. They hold the mean
false discovery proportion (FDP) and true positive proportion (TPP) of each method over the SNR $\times$ `ar_coef`
grid, and the comma-separated values (CSV) file adds their standard deviations. The file names and the CSV schema
match the C++ and Python demos.

Figures in Portable Network Graphics (PNG) and Portable Document Format (PDF) go to `simulation_results/plots/`,
with the same file names as in the Python suite.

---

## Running the Demo

```bash
demo=R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_07_mc_sim_arma_blocks/demo_trex_gvs_07_mc_sim_arma_blocks.R

Rscript "$demo"      # full run, TREX_NUM_WORKERS workers (default 8)
Rscript "$demo" 8    # 8 worker cores
```

Monte Carlo trials run in parallel on a parallel socket (PSOCK) cluster of spawned R worker processes
([`../../../mc_runtime.R`](../../../mc_runtime.R)). The worker count is the first command-line argument, else
`TREX_NUM_WORKERS`, else 8, and `0` uses every core. Every finished grid cell is appended to a cache in
`simulation_results/cache/`, so a stopped run resumes where it left off, and a completed run removes the cache. The
script also sources [`../../../support_generators.R`](../../../support_generators.R),
[`../../../simulation_utils.R`](../../../simulation_utils.R), [`../trex_gvs_dgps.R`](../trex_gvs_dgps.R) and
[`../trex_gvs_sim_common.R`](../trex_gvs_sim_common.R). The recorded run (library commit `44cf1ab`, 8 workers) took
about 2.7 hours. For a quick smoke test, reduce `num_MC` in `CFG` and the grids in `PRESET` at the top of the script.

Afterwards, regenerate the figure from the CSV file with [`generate_plots.sh`](generate_plots.sh), which drives the
suite's ggplot2 plotting module [`../trex_gvs_plt_utils.R`](../trex_gvs_plt_utils.R) through `Rscript` (set `RSCRIPT`
to use another R installation).

---

## Simulation Results

The numbers and figure below are this R run's own results (200 MC trials per cell). The false discovery rate (FDR) and
the true positive rate (TPR) of a cell are the means of FDP and TPP over its trials. R draws its own random numbers,
so the results agree with the
[C++ reference](../../../../../cpp/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_07_mc_sim_arma_blocks/README.md#simulation-results)
and the
[Python port](../../../../../Python/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_07_mc_sim_arma_blocks/README.md#simulation-results)
up to Monte Carlo noise, not exactly. A cell-by-cell comparison with the C++ run, with $z$ the difference divided by
the combined Monte Carlo standard error of the two 200-trial means, finds 4 cells with $|z| > 3$ among the 180 cells
of this demo (FDR and TPR of every solver at every grid point). They are two grid points counted twice, because
TENET and TENET_AUG return identical numbers. Both are TENET FDR values that R puts above C++, $0.086$ against
$0.071$ at $\mathrm{SNR} = 0.5$ and `ar_coef` $= 0.95$ ($z = 3.02$), and $0.064$ against $0.050$ at
$\mathrm{SNR} = 2$ and `ar_coef` $= 0.5$ ($z = 3.14$). Python lies between the two in both cells ($0.077$ and
$0.055$), and all three stay below the target. Across the 30 grid points R's TENET FDR lies above the C++ value in 29,
by $0.006$ on average. Every grid point reuses the same trial seeds, so the cells of one run are correlated and a
shared sign is weaker evidence than it would be for independent cells. Over the whole T-Rex+GVS suite these are the
only 4 of 1,698 cells (0.24 %) with $|z| > 3$, a rate within Monte Carlo noise.

- The aggregate TPR **levels off** instead of approaching 1. At $\mathrm{SNR} = 5$, TENET's TPR rises over `ar_coef`
  $\in \lbrace0.3, 0.5, 0.7\rbrace$ ($0.11 \to 0.36 \to 0.59$) and then stays at $0.58\text{–}0.60$. From `ar_coef`
  $= 0.8$ on it lies at $0.55\text{–}0.60$ for every $\mathrm{SNR} \ge 0.5$ (and $0.48\text{–}0.55$ at
  $\mathrm{SNR} = 0.2$), roughly 83 to 91 of the 150 active variables. That is fewer than the 100 variables of the AR(1)
  and ARMA(2,1) blocks together. A plausible reading is that the MA(3) block, whose correlation vanishes beyond lag 3,
  adds few selections, but the demo records no per-block rates.
- **TENET** and **TENET_AUG** return identical numbers (means and standard deviations agree in every cell) and keep
  the FDR below the $\mathrm{tFDR} = 0.1$ target in every cell (no cyan outlines). Their FDR grows with `ar_coef`,
  from at most $0.047$ at `ar_coef` $= 0.3$ to $0.081\text{–}0.086$ at $0.9$ and $0.95$, with the maximum $0.086$ at
  $\mathrm{SNR} = 2$, `ar_coef` $= 0.95$. The trap block (AR(1), same `ar_coef` as the active AR(1) block) is
  structurally similar to an active block and is the main test of false-group exclusion here. The FDR stays below
  target despite it, although the demo does not record which false discoveries fall into the trap.
- **TIENET_AUG** is more conservative, with an FDR of at most $0.046$. Its TPR peaks at $0.51$ (`ar_coef` $= 0.7$,
  $\mathrm{SNR} = 5$), reaches up to $0.45$ at `ar_coef` $= 0.8$, and collapses to $0.006\text{–}0.007$ from
  `ar_coef` $= 0.9$ on, while TENET keeps its $0.55\text{–}0.60$ there. In the six collapse cells from
  $\mathrm{SNR} = 1$ on, the TPR is exactly $1/150$ and the FDR exactly $0$, both with standard deviation $0$, so
  every trial selected a single variable, an active one. The collapse sets in where the ARMA(2,1) block turns
  non-stationary and nearly collinear (see the DGP section). A plausible reading is that TIENET_AUG meets within-block
  correlation above `corr_max = 0.98` as in the $\rho = 0.99$ columns of Demos 05 and 06, where it collapses as
  well, although the ARMA(2,1) block is far more extreme (pairwise correlation of at least $0.999999$).
- Against Demo 06 (uniform AR(1) blocks, R run), the ARMA design has a lower TENET ceiling, a top TPR of $0.605$ here
  against $0.961$ there. At a coefficient of $0.9$ and $\mathrm{SNR} \ge 0.5$, Demo 06 reaches $0.63\text{–}0.95$
  against $0.55\text{–}0.58$ here. At coefficients up to $0.7$, in contrast, the ARMA design gives the higher TPR at all
  15 shared grid points ($0.36$ against $0.12$ at $0.5$ and $\mathrm{SNR} = 5$, for example). A plausible reading is
  that its 80-variable ARMA(2,1) block is more strongly correlated at short lags than an AR(1) block with the same
  coefficient (lag-1 autocorrelation $0.62$, $0.77$ and $0.91$ against $0.3$, $0.5$ and $0.7$), because its second AR
  term and its MA term add correlation.

TPR (top) and FDR (bottom) heatmaps over the SNR $\times$ `ar_coef` grid, one column per solver
(row axis = SNR, column axis = AR coefficient). FDR cells above the $\mathrm{tFDR} = 0.1$ target
would be outlined in cyan (none are here).

![ARMA blocks SNR × AR-coefficient sweep](simulation_results/plots/gvs_arma_blocks_2d.png)

---

## References

1. Machkour, J., Muma, M., & Palomar, D. P., "False Discovery Rate Control for Grouped Variable Selection
   in High-Dimensional Linear Models using the T-Knock Filter.", European Signal Processing Conference (EUSIPCO), 2022,
   pp. 892–896, European Association for Signal Processing (EURASIP).
   [DOI-Link](https://doi.org/10.23919/EUSIPCO55093.2022.9909883)
2. Machkour, J., Muma, M., & Palomar, D. P., "The Informed Elastic Net for Fast Grouped Variable Selection and
   FDR Control in Genomics Research.", Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP),
   2023, pp. 466–470, Institute of Electrical and Electronics Engineers (IEEE).
   [DOI-Link](https://doi.org/10.1109/CAMSAP58249.2023.10403489)
3. Zou, H., & Hastie, T. (2005). "Regularization and variable selection via the elastic net." *Journal of the Royal
   Statistical Society: Series B (Statistical Methodology)*, 67(2), pp. 301–320.
   [DOI-Link](https://doi.org/10.1111/j.1467-9868.2005.00503.x)
