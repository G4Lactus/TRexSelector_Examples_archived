# Demo 05: T-Rex+DA+BT on Heavy-Tailed Block + White-Noise Data

## Purpose

Combine Demo 02's white-noise dilution with Demo 04's heavy-tailed robustness study. Heavy-tailed first-order
autoregressive (AR(1)) Toeplitz blocks sit in a larger design padded with independent and identically distributed
(i.i.d.) Student-t white-noise columns, with total dimension held fixed. It decomposes how heavy tails and white-noise
dilution each affect false discovery rate (FDR) control and power for the `BT` method.

## Data-generating process

`dgp_ht_block_white_snr` builds `M` blocks of size `Q` with heavy-tailed (t(nu)) AR(1)-Toeplitz within-block
covariance (one chi-squared mixing variable per (row, block)), plus `p_white` appended i.i.d. t(nu) white-noise
columns, with the total dimension fixed at `p_total = 500`. Predictors X are always heavy-tailed, while the response
noise is Gaussian (`s1_Gauss`) or t(nu) (`s2_Heavy`) by scenario. Active variables (`s = M`, one per AR(1) block) lie
in the AR(1) part only. Base config: `n=150, p_total=500, M=5, Q=5 (p_ar=25, p_white=475, s=5)`,
`amplitude=1.0, rho=0.8, nu=3.0, SNR=2.0, tFDR=0.2, K=20, num_MC=200, seed=2026`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "BT", hc_linkage = ...), control = trex_control(...))`
with the default solver set from `da_sim_common.R`: `TREX-DA+BT: TLARS` / `TREX-DA+BT: TAFS` / `TREX-DA+BT: TOMP`.
Outer loops over noise scenario `{s1_Gauss = Gaussian noise, s2_Heavy = heavy noise}` x linkage
`{Single, Complete, Average}`.

## What it runs

- Part 1: signal-to-noise ratio (SNR) sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}` with `p_ar=25, p_white=475`.
- Part 2: rho sweep `seq(0.0, 0.9, by = 0.1)`.
- Part 3: Q sweep `{5, 10, ..., 50}` (`p_ar = M*Q`, `p_white = p_total - p_ar`).
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}`.
- Part 5: target FDR sweep `seq(0.05, 0.50, by = 0.05)`.

Result stems (`simulation_results/data/`):
`da_trex_mc_da_ht_blocks_plus_white_{snr,rho,Q,M,tFDR}_{s1_Gauss,s2_Heavy}_{Single,Complete,Average}`
(each as `.txt` + `.csv`).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_05_mc_sim_ht_blocks_plus_white/demo_trex_da_05_mc_sim_ht_blocks_plus_white.R      # full run
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_05_mc_sim_ht_blocks_plus_white/demo_trex_da_05_mc_sim_ht_blocks_plus_white.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` as `.txt` and `.csv` files. The
script sources the consolidated `../trex_da_dgps.R`, the shared
`../da_sim_common.R`, and `../../../support_generators.R`.

## Interpretation

- Unlike Demo 04 (heavy-tailed blocks without white-noise dilution), FDR here tends to decrease as SNR increases and
  end up quite low, closer to Demo 02's white-noise-diluted pattern — the large white-noise padding helps FDR control
  even under heavy tails.
- The true positive rate (TPR) typically reaches higher levels than Demo 04's plateau at the same nominal SNR, again
  consistent with the white-noise dilution making the true AR(1) blocks comparatively easier to isolate than in
  Demo 04's smaller all-block design.
- Read this demo as the combined-effects check and compare it against Demo 03 (Gaussian, no white noise), Demo 02
  (Gaussian, with white noise), and Demo 04 (heavy-tailed, no white noise) to decompose how heavy tails and white-noise
  dilution each contribute. The target-FDR sweep (Part 5) checks target tracking on this diluted design.
