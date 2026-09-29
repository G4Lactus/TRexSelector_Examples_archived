# Demo 02: T-Rex+DA+BT on Block AR(1) + White-Noise Data

## Purpose

Extend Demo 03's block-diagonal first-order autoregressive (AR(1)) design by appending a large block of independent and
identically distributed (i.i.d.) white-noise columns, testing whether `BT`-based Dependency-Aware T-Rex (DA-TRex) can
still isolate the true AR(1) blocks when they are a small fraction of a much larger, otherwise-uncorrelated design.
Holding the total dimension fixed lets M and Q vary without changing p_total.

## Data-generating process

`dgp_ar1_block_white_snr` builds `M` AR(1) blocks of size `Q` (`p_ar = M * Q`), plus `p_white` appended i.i.d. N(0,1)
columns, with total dimension fixed at `p_total = p_ar + p_white = 1000`. Active variables (`s = M`, one per AR(1)
block) lie only in the AR(1) part. Base config: `n=300, p_total=1000, M=5, Q=5 (p_ar=25, p_white=975, s=5)`,
`amplitude=1.0, rho=0.7, SNR=2.0, tFDR=0.1, K=20, num_MC=200, seed=2026`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "BT", hc_linkage = ...), control = trex_control(...))`
with the default solver set from `da_sim_common.R`: `TREX-DA+BT: TLARS` / `TREX-DA+BT: TAFS` / `TREX-DA+BT: TOMP`.
Each part loops over linkage in `{Single, Complete, Average}`.

## What it runs

- Part 1: signal-to-noise ratio (SNR) sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}` with `p_ar=25, p_white=975` fixed,
  outer loop over linkage.
- Part 2: rho sweep `seq(0.0, 0.9, by = 0.1)`, outer loop over linkage.
- Part 3: Q sweep `{5, 10, ..., 50}` with `p_ar = M*Q` varying over 25..250, `p_white = 1000 - p_ar`, and `s = M = 5`.
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}`, in which `p_ar`, `p_white`, and `s = M` all vary.

Result stems (`simulation_results/data/`):
`da_trex_mc_da_ar1_blocks_plus_white_{snr,rho,Q,M}_{Single,Complete,Average}`
(each as `.txt` + `.csv`).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_02_mc_sim_ar1_blocks_plus_white/demo_trex_da_02_mc_sim_ar1_blocks_plus_white.R      # full run
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_02_mc_sim_ar1_blocks_plus_white/demo_trex_da_02_mc_sim_ar1_blocks_plus_white.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` as `.txt` and `.csv` files. The
script sources the consolidated `../trex_da_dgps.R`, the shared
`../da_sim_common.R`, and `../../../support_generators.R`.

## Interpretation

- The false discovery rate (FDR) is typically very tightly controlled here (well below the lower 0.1 target) and lower
  than in Demo 03's block-only design, because the many uncorrelated white-noise columns dilute the design and leave
  much less correlated, false-positive-prone structure to latch onto.
- The true positive rate (TPR) reaches a strong plateau similar to Demo 03, suggesting the added white-noise columns do
  not meaningfully hurt detection of the true AR(1) blocks, consistent with Terminating Random Experiments (T-Rex)'s
  robustness to null variables.
- Compare directly against Demo 03 (no white noise, same block core) to isolate the effect of a larger, mostly-null p.
  The Q/M/linkage sweeps remain the key diagnostics for block size, block count, and linkage choice.
