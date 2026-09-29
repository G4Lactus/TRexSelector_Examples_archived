# Demo 06: T-Rex+DA on Multi-Level Prior-Groups Data

## Purpose

Test Dependency-Aware T-Rex (DA-TRex) on the structurally richest dependency model in this suite, a three-level
hierarchical latent-factor structure with a non-exchangeable Toeplitz leaf layer, using explicit prior group labels
rather than an automatically discovered dendrogram cut. This isolates the difficulty of the multi-level
dependency-aware correction itself from any clustering-quality artifact.

## Data-generating process

`dgp_groups_toeplitz_leaf` draws from a multi-level latent-factor model where coarser levels use shared factors and the
finest level uses a Toeplitz-correlated leaf block (`Cov(u_r, u_s) = phi^|r-s|`). There are three nested grouping
levels, with groups of 10, 50, and 250 (fine to coarse). Config:
`n=300, p=1000, s=10, amplitude=3.0, tFDR=0.2, K=20, num_MC=200, seed=2026`, with cumulative correlation levels
`rho_levels = {0.55, 0.25, 0.10}`, Toeplitz leaf decay `phi_leaf = 0.50`, and Random support.

## DA method

In the prior-groups path, `trex_da_control(prior_groups = groups)` supplies the 3-level hierarchy explicitly. The core
routes on `prior_groups` being non-empty, so the `da_method` value is irrelevant here, and this is not the hierarchical
agglomerative clustering (HAC) / binary-tree (BT) discovery path. The labels act as merge CONSTRAINTS. The selector
sub-clusters within the labels' finest common refinement and calibrates over a data-driven rho grid ending at the
conservative rho = 1 anchor. Solver rows from `da_sim_common.R`: `TREX-DA: TLARS` / `TREX-DA: TAFS` / `TREX-DA: TOMP`.

## What it runs

- Part 1: signal-to-noise ratio (SNR) sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` (n=300, p=1,000, s=10, Random support
  redrawn per trial), 200 Monte Carlo trials per point.

Result stem (`simulation_results/data/`): `da_trex_mc_da_groups_toeplitz_leaf`
(`.txt` + `.csv`).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_06_mc_sim_groups/demo_trex_da_06_mc_sim_groups.R      # full run
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_06_mc_sim_groups/demo_trex_da_06_mc_sim_groups.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` as `.txt` and `.csv` files. The
script sources the consolidated `../trex_da_dgps.R`, the shared
`../da_sim_common.R`, and `../../../support_generators.R`.

## Interpretation

- The three-level nested dependency structure is genuinely harder for DA-TRex to control than the single-level
  structures in the other demos, so the realized false discovery rate (FDR) can overshoot the 0.2 target more than
  elsewhere in this suite.
- Because prior groups are supplied directly rather than HAC-discovered, any residual FDR inflation reflects genuine
  difficulty in the underlying multi-level correction, not a clustering-quality confound — which differentiates Demo 06
  from the BT-based Demos 02-05.
- The demo currently exercises only the SNR sweep. Extending it with rho-level or group-size sweeps (analogous to
  Demos 02-05's Q/M sweeps) would help clarify which aspect of the nested structure drives the elevated FDR.
