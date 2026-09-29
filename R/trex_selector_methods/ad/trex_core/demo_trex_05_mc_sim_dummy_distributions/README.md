# Demo 05: Dummy Distribution Comparison

## Purpose

Compare the T-Rex **dummy distributions** — the distribution from which the injected
dummy (null) variables are drawn — across three base solvers
(**TLARS** — equiangular LARS path; **TOMP** and **TAFS** with `rho_afs = 0.3` —
greedy forward selection), with the L-loop strategy fixed at **STANDARD**. Since dummies are centered and
L2-normalized before entering the solver, their *scale* is immaterial; the demo probes
whether the *shape* (tail weight, discreteness, skewness, sparsity, cross-column
dependence) affects FDR, TPR, and the calibrated dummy multiplier `L` / stopping time `T`.

## Data-generating process

`dgp_gauss_snr` — i.i.d. Gaussian design, `y = X beta + eps` calibrated to SNR.
High-dimensional config: `n = 300`, `p = 1000`, `s = 10`, `tFDR = 0.1`, fixed
coefficients `beta_j = 1`. The support is a random size-`s` subset redrawn per trial
via `make_support_random(s, p, seed = trial_seed)` (the `block_support = TRUE`
contiguous-block variant is available but its call is commented out). SNR grid
`{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}` (7 values, mirrors the current C++ grid),
`num_MC = 200` per solver x distribution x SNR.

The *design* is the suite's shared one and stays fixed across the sweep —
independent columns, empirical column-correlation essentially the identity — so
the only thing varying below is the distribution of the injected dummies:

![Design-matrix correlation of the shared i.i.d. Gaussian DGP (natural order, no clustering)](../dgp_corr_plots/trex_gauss.png)

The solver sweep (`SOLVERS`) is the outer loop; one `.txt` + `.csv` file pair is
written per solver, with the lowercased solver name in the stem
(`..._{tlars,tomp,tafs}_random_support`).

## Distributions compared

`DISTRIBUTIONS` — **12 rows**, each built with `trex_dummy_distribution()` and passed
to `trex_control(dummy_distribution = ...)`:

- Normal — N(0, 1); the reference choice.
- Uniform — U(-sqrt(3), sqrt(3)); bounded, light tails.
- Rademacher — {-1, +1} equiprobable; discrete two-point.
- StudentT_df3 / StudentT_df5 — heavy tails (unit-variance scaled).
- Laplace — double exponential; heavier-than-normal tails.
- Gumbel — extreme-value; skewed (centered at 0).
- Triangle — bounded, symmetric; (-sqrt(6), 0, sqrt(6)).
- Logistic — tail weight between normal and Laplace.
- Mammen — asymmetric golden-ratio two-point distribution.
- SparseRad_s0.1 — constrained sparse Rademacher, sparsity `s = 0.1`.
- UnifSphere_d5 — uniform on the unit 5-sphere; consecutive 5-column groups are
  dependent (norm constraint), a deliberate probe of the exchangeability axis.

**Note on UniformSphere**: the library requests dummies in multiples of `p`, so the
sphere dimension must divide `p`; the default `dim = 3` would throw at `p = 1000`,
hence `dim = 5`.

Shared control per solver x distribution: `K = 20`, `max_dummy_multiplier = 10`,
`use_max_T_stop = TRUE`, `opt_threshold = 0.75`, `lloop_strategy = "STANDARD"`,
`tloop_max_stagnant_steps = 5`, `use_openmp = FALSE`, `use_memory_mapping = FALSE`,
`dummy_distribution` swept; `tloop_stagnation_stop` is left unset (AUTO — the
per-solver-family default). Reports averaged FDR, TPR, Avg L, Avg T per distribution x
SNR (`track_L_T = TRUE`). MC trials run in parallel over a hardcoded `num_cores <- 6L`.

## Expected pattern

The C++ reference run shows all i.i.d. unit-variance rows FDR-equivalent at or below
`tFDR` (isolated exceedances at very low SNR are small-`num_MC` noise), statistically
indistinguishable TPR, and no systematic deviation even for the dependent
UnifSphere_d5 rows.

## Results

This R run reproduces that pattern. `./generate_plots.sh` renders a cross-solver
comparison grid plus one figure set per base solver through the shared ggplot2
plotter [`../trex_plt_utils.R`](../trex_plt_utils.R), into
`simulation_results/plots/` as PNG + PDF (this port renders no interactive HTML).

### Cross-solver comparison grid

The headline figure — a 2x3 layout, rows = metric (FDR / TPR), columns = base
solver (TLARS / TOMP / TAFS), with all twelve dummy distributions drawn in every
panel. Colours are keyed to the distribution, so a distribution can be tracked
across solvers; the three FDR panels share one y-scale and the three TPR panels
share `[0, 1]`:

![FDR/TPR per base solver per dummy distribution, 2x3 comparison grid](simulation_results/plots/demo_trex_05_dummy_distributions_grid.png)

Within each column the twelve curves are essentially one thick line. At SNR 1.0
— at or next to the FDR peak for every solver — the entire twelve-distribution
FDR spread is 0.0475 (Triangle) to 0.0575 (StudentT_df5) under TLARS, 0.0256
(Rademacher) to 0.0325 (Gumbel) under TOMP, and 0.0311 (Gumbel) to 0.0429
(Laplace) under TAFS: every value well under the target `tFDR = 0.1`, and the
within-solver spread is of the same order as the gap between solvers. TPR at the same
point is equally tight (TLARS 0.723–0.741, TOMP 0.859–0.875, TAFS 0.865–0.885),
and every distribution reaches 1.0000 by SNR 5.0 under all three solvers. The
*level* is solver-dependent — TLARS carries the familiar FDR hump around SNR
0.6–1.0 while the greedy TOMP / TAFS run flatter and lower — but the choice of
dummy distribution moves nothing systematically, not even the deliberately
dependent **UnifSphere_d5** row (TLARS 0.0510, TOMP 0.0277, TAFS 0.0351 at SNR
1.0, all mid-pack). Given that dummies are centered and L2-normalized before
entering the solver, this is the demo's conclusion: the dummy distribution is
effectively a free parameter, and the Normal default is justified.

Note that the ranking of distributions is *not* stable across solvers — the
TLARS minimum (Triangle) is not the TOMP minimum (Rademacher) or the TAFS
minimum (Gumbel), and the TAFS maximum (Laplace) is mid-pack elsewhere. With
`num_MC = 200` the residual spread is Monte Carlo noise, not an ordering.

### Per-solver figure sets

The usual overview (all twelve distributions on one FDR | TPR pair, log-scaled
x-axis) plus a de-cluttered grouped 2x2 view that splits the distributions across
two columns. The CSV's "solver" column holds the distribution name, so the
distributions are the legend.

**TLARS**:

![TLARS: FDR and TPR versus SNR for all 12 dummy distributions](simulation_results/plots/demo_trex_05_dummy_distributions_results_n300_p1000_tlars_random_support_fdr_tpr_vs_snr.png)

![TLARS: FDR and TPR versus SNR, distributions split across a 2x2 grid](simulation_results/plots/demo_trex_05_dummy_distributions_results_n300_p1000_tlars_random_support_fdr_tpr_vs_snr_grouped.png)

**TOMP**:

![TOMP: FDR and TPR versus SNR for all 12 dummy distributions](simulation_results/plots/demo_trex_05_dummy_distributions_results_n300_p1000_tomp_random_support_fdr_tpr_vs_snr.png)

![TOMP: FDR and TPR versus SNR, distributions split across a 2x2 grid](simulation_results/plots/demo_trex_05_dummy_distributions_results_n300_p1000_tomp_random_support_fdr_tpr_vs_snr_grouped.png)

**TAFS** (`rho_afs = 0.3`):

![TAFS: FDR and TPR versus SNR for all 12 dummy distributions](simulation_results/plots/demo_trex_05_dummy_distributions_results_n300_p1000_tafs_random_support_fdr_tpr_vs_snr.png)

![TAFS: FDR and TPR versus SNR, distributions split across a 2x2 grid](simulation_results/plots/demo_trex_05_dummy_distributions_results_n300_p1000_tafs_random_support_fdr_tpr_vs_snr_grouped.png)

The grouped views are the ones to read when the overview curves overlap: with the
twelve rows split across two columns on a shared y-scale, "no systematic
deviation" is easier to confirm than on a single crowded panel.

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_core/demo_trex_05_mc_sim_dummy_distributions/demo_trex_05_mc_sim_dummy_distributions.R
```

Results are written to this demo's own `simulation_results/data/` folder as a `.txt` table
and a tidy `.csv` (the "solver" column holds the distribution name), stems
`demo_trex_05_dummy_distributions_results_n300_p1000_{tlars,tomp,tafs}_random_support`. The script sources
`../../support_generators.R`, `../../simulation_utils.R`, and `../trex_sim_utils.R`.
