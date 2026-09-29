# Demo 07: T-Rex+DA+NN on Banded (Nearest-Neighbor) Data

## Purpose

Study the `NN` DA method, which targets a banded MA(kappa) covariance where only
the kappa nearest-neighboring columns are correlated and correlation is exactly
zero beyond that band. Parts 1-2 cover the correctly-specified case (NN correction
on NN data); Part 3 applies the same NN correction to AR(1) data — a deliberate
misspecification stress test (formerly its own demo, 03b, in the legacy R
numbering).

## Data-generating process

`dgp_nn_snr` builds `X[i,j] = sum_{l=0..kappa} theta_l * eta[i,j+l]` with
`theta_l ~ rho^l` normalized to unit variance, so `Sigma[j,k] = 0` once
`|j-k| > kappa`. Part 3 uses `dgp_ar1_snr` (the same AR(1) Toeplitz design as
Demo 01) with the NN correction. Monte Carlo config: `n=300, p=1000, s=10,
amplitude=3.0, kappa=3, rho=0.7, SNR=2.0, tFDR=0.2, K=20, num_MC=200, seed=2026,
max_gap=20`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "NN"),
control = trex_control(...))` with the default solver set from
`da_sim_common.R`: `TREX-DA+NN: TLARS` / `TREX-DA+NN: TAFS` /
`TREX-DA+NN: TOMP`. Part 1 also runs the base (no-DA) `TREX: TLARS/TAFS/TOMP`
comparison rows.

## What it runs

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` on NN data (kappa=3, rho=0.7), CappedSpread(max_gap=20) and Random support, DA + base T-Rex rows.
- Part 2: 2D kappa x rho sweep on NN data; `kappa in {1,2,3,5,7,10,15}`, `rho in {0.1,0.3,0.5,0.7,0.8,0.9}` at SNR=2.0, CappedSpread and Random sub-sections.
- Part 3: 2D SNR x rho sweep on AR(1) data with DA+NN; `SNR in {0.1,0.2,0.5,1.0,2.0,5.0}`, `rho in seq(0.1, 0.9, by = 0.1)`, CappedSpread and Random sub-sections.

Result stems (`simulation_results/data/`): `da_trex_mc_da_nn_snr_{capped,random}`,
`da_trex_mc_da_nn_kappa_rho`, `da_trex_mc_da_nn_ar_snr_rho`
(each as `.txt` + `.csv`).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_07_mc_sim_nn/demo_trex_da_07_mc_sim_nn.R      # full run
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_07_mc_sim_nn/demo_trex_da_07_mc_sim_nn.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` as TXT and CSV files; the
script sources the consolidated `../trex_da_dgps.R`, the shared
`../da_sim_common.R`, and `../../support_generators.R`.

## Interpretation

- Because NN correlation has a hard cutoff at distance kappa, DA-NN's correction is most effective when active variables are spaced further apart than kappa, analogous to Demo 01's gap-vs-kappa story but with a sharp rather than decaying profile.
- The 2D kappa x rho sweep (Part 2) is the most informative view: as kappa grows the effective correlation neighborhood widens, and DA-NN needs progressively wider correction windows with a corresponding TPR cost.
- Part 3 directly tests robustness to model misspecification: NN assumes correlation vanishes beyond kappa columns, while AR(1) correlation decays but never truly reaches zero, so some residual correlation should leak through as extra false discoveries relative to Demo 01's matched AR1-on-AR1 results. If DA-NN performs nearly as well as DA-AR1 there, the method is reasonably robust to this mismatch; a large gap would indicate sensitivity to matching the true dependency structure.
