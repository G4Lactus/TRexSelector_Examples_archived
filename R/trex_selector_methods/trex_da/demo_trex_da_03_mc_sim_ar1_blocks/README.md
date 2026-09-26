# Demo 03: T-Rex+DA+BT on Block-Diagonal AR(1) Data

## Purpose

Introduce the `BT` (binary-tree / dendrogram HAC clustering) DA method on a clean
block-diagonal AR(1) design, and study sensitivity to the HAC linkage method
(Single / Complete / Average) across four sweep dimensions: SNR, rho, block size Q,
and number of blocks M.

## Data-generating process

`dgp_ar1_block_snr` partitions `p` predictors into `M` independent AR(1) blocks of
size `Q` (`p = M * Q`); within each block correlation follows an AR(1) Toeplitz
recursion, with zero correlation between blocks. Support is one active variable per
block (`make_support_one_per_block`), so `s = M`. Base config: `n=150, M=5,
Q=5 (p=25, s=5), amplitude=1.0, rho=0.7, SNR=2.0, tFDR=0.2, K=20, num_MC=200,
seed=2026`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "BT",
hc_linkage = ...), control = trex_control(...))` with the default solver set
from `da_sim_common.R`: `TREX-DA+BT: TLARS` / `TREX-DA+BT: TAFS` /
`TREX-DA+BT: TOMP`. Each part loops over linkage in
`{Single, Complete, Average}`.

## What it runs

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`, outer loop over linkage.
- Part 2: rho sweep `seq(0.0, 0.9, by = 0.1)`, outer loop over linkage.
- Part 3: Q sweep `{5, 10, ..., 50}`; `p = M*Q` varies from 25 to 250, `s = M = 5` fixed.
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}`; `p = M*Q` and `s = M` both vary.

Result stems (`simulation_results/data/`):
`da_trex_mc_da_ar1_blocks_{snr,rho,Q,M}_{Single,Complete,Average}`
(each as `.txt` + `.csv`).

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_03_mc_sim_ar1_blocks/demo_trex_da_03_mc_sim_ar1_blocks.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_03_mc_sim_ar1_blocks/demo_trex_da_03_mc_sim_ar1_blocks.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` as TXT and CSV files; the
script sources the consolidated `../trex_da_dgps.R`, the shared
`../da_sim_common.R`, and `../../support_generators.R`.

## Interpretation

- FDR stays broadly near the tFDR=0.2 target while TPR climbs to a plateau (rather than reaching ~1.0 as in the unblocked AR(1) case), reflecting that recovering every one of the M active blocks consistently is harder than recovering isolated active variables.
- The Q and M sweeps are the most novel diagnostics: as block size Q or block count M grows the problem becomes harder, so check whether TPR degrades gracefully or drops sharply beyond some threshold.
- The linkage loop answers a practical question for BT users: does Single/Complete/Average linkage materially change FDR/TPR here, or is performance robust to that choice?
