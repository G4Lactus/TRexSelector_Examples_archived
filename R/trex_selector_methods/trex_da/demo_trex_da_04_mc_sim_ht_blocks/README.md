# Demo 04: T-Rex+DA+BT on Heavy-Tailed Block Data

## Purpose

Stress-test `BT`-based DA-TRex under heavy-tailed (Student-t, nu=3) predictors on a
block-diagonal AR(1)-Toeplitz design. It studies both statistical robustness
(Gaussian vs. heavy-tailed response noise) and sensitivity across five sweep
dimensions plus a dedicated linkage-method sweep.

## Data-generating process

`dgp_block_toeplitz_hvt_snr` builds `M` blocks of size `Q` (`p = M * Q`) with
within-block AR(1)-Toeplitz covariance; predictors X are always drawn multivariate
Student-t (nu=3, one chi-squared mixing variable per (row, block)), and the
response noise is Gaussian (`s1_Gauss`) or Student-t (`s2_Heavy`) depending on the
scenario. Support is one active variable per block (`s = M`). Base config: `n=150,
M=5, Q=5 (p=25, s=5), amplitude=1.0, rho=0.8, nu=3.0, SNR=2.0, tFDR=0.2, K=20,
num_MC=200, seed=2026`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "BT",
hc_linkage = ...), control = trex_control(...))` with the default solver set
from `da_sim_common.R`: `TREX-DA+BT: TLARS` / `TREX-DA+BT: TAFS` /
`TREX-DA+BT: TOMP`. Outer loops over noise scenario
`{s1_Gauss = Gaussian noise, s2_Heavy = heavy noise}` x linkage
`{Single, Complete, Average}`.

## What it runs

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- Part 2: rho sweep `seq(0.0, 0.9, by = 0.1)`.
- Part 3: Q sweep `{5, 10, ..., 50}` (`p = M*Q`).
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}`.
- Part 5: tFDR sweep `seq(0.05, 0.50, by = 0.05)`.
- Part 6: linkage sweep, codes `{1=Single, 2=Complete, 3=Average}`, at base parameters.

Result stems (`simulation_results/data/`):
`da_trex_mc_da_ht_blocks_{snr,rho,Q,M,tFDR}_{s1_Gauss,s2_Heavy}_{Single,Complete,Average}`
and `da_trex_mc_da_ht_blocks_linkage_{s1_Gauss,s2_Heavy}` (each as `.txt` + `.csv`).

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_04_mc_sim_ht_blocks/demo_trex_da_04_mc_sim_ht_blocks.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_04_mc_sim_ht_blocks/demo_trex_da_04_mc_sim_ht_blocks.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` as TXT and CSV files; the
script sources the consolidated `../trex_da_dgps.R`, the shared
`../da_sim_common.R`, and `../../support_generators.R`.

## Interpretation

- Under heavy-tailed data, FDR tends to hover around or somewhat above the tFDR=0.2 target — less tightly controlled than the Gaussian block demos (02-03) — reflecting the added estimation difficulty from occasional large t(3) outliers.
- TPR typically plateaus lower than Demo 03's Gaussian block-AR(1) plateau, consistent with heavy tails making the true signal harder to separate from noise even at high nominal SNR.
- Compare the Gaussian-noise vs. heavy-noise scenarios directly to isolate the cost of heavy tails; the dedicated linkage sweep (Part 6) matters more here since correlation estimates feeding HAC clustering are noisier, and the tFDR sweep (Part 5) checks whether realized FDR tracks the target across levels.
