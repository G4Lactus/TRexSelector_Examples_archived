# DA-TRex: Dependency-Aware T-Rex — R Demonstration Suite

> **Status**: NEW — migrated to the **TRexSelectorNeo** R6 application programming interface (API)
> (`TRexDASelector$new(..., da_control = trex_da_control(...))$select()`, with `trex_fdp()` / `trex_tpp()`). No
> longer depends on the `TRexSelector` package from the Comprehensive R Archive Network (CRAN).

## Purpose

R demos and Monte Carlo (MC) simulations for the **Dependency-Aware T-Rex (DA-TRex)** selector across correlated
designs: first-order autoregressive (AR(1)) Toeplitz, equicorrelation, nearest-neighbor (NN) with a banded
moving-average (MA) structure, block-diagonal AR(1) (with and without white-noise augmentation), and heavy-tailed
Student-t blocks.

The folder layout mirrors the C++ / Python suites, with one subfolder per demo, plus a single consolidated
data-generating process (DGP) file and the shared MC / writer infrastructure.

```txt
trex_da/
  ├── README.md
  ├── da_sim_common.R                  <- solver set, MC runners, sweeps, writers
  ├── trex_da_dgps.R                   <- all scenario DGPs (one file)
  ├── demo_trex_da_01_mc_sim_ar1/
  │   ├── demo_trex_da_01_mc_sim_ar1.R
  │   ├── README.md
  │   └── simulation_results/data/     <- da_trex_mc_*.{txt,csv}
  ├── demo_trex_da_02_mc_sim_ar1_blocks_plus_white/
  │   └── ...
  └── demo_trex_da_08_mc_sim_equi_and_bt/
```

Dependencies: `TRexSelectorNeo`, `parallel` (MC loops). Each demo sources the
consolidated [trex_da_dgps.R](trex_da_dgps.R) (all scenario generators in one
file, mirroring the C++ `dgp_generators.hpp` / Python `dgp_generators.py`), the
shared [da_sim_common.R](da_sim_common.R) (default solver set, parallel MC
runners, parameter sweeps, and the `.txt`/`.csv` result writers, mirroring the Python
`da_sim_common.py` / C++ `trex_da_sim_common.hpp`), and
[../../support_generators.R](../../support_generators.R).

---

## Demos

| File | Description | cpp counterpart |
|---|---|---|
| [demo_trex_da_01_mc_sim_ar1](demo_trex_da_01_mc_sim_ar1/) | AR(1) Toeplitz DGP: signal-to-noise ratio (SNR) sweep and rho sweep (CappedSpread and Random support, dependency-aware (DA) + base T-Rex rows), plus a 2D gap x rho sweep exploring the kappa boundary | `cpp/trex_selector_methods/ad/trex_da/demo_trex_da_01_mc_sim_ar1` |
| [demo_trex_da_02_mc_sim_ar1_blocks_plus_white](demo_trex_da_02_mc_sim_ar1_blocks_plus_white/) | AR(1)-block + white-noise DGP (fixed p_total = 1,000, p_white = p_total − M*Q): SNR, rho, Q, and M sweeps, each with an outer linkage loop | `cpp/trex_selector_methods/ad/trex_da/demo_trex_da_02_mc_sim_ar1_blocks_plus_white` |
| [demo_trex_da_03_mc_sim_ar1_blocks](demo_trex_da_03_mc_sim_ar1_blocks/) | Pure AR(1)-block DGP (p = M*Q): SNR, rho, Q, and M sweeps, each with an outer linkage loop | `cpp/trex_selector_methods/ad/trex_da/demo_trex_da_03_mc_sim_ar1_blocks` |
| [demo_trex_da_04_mc_sim_ht_blocks](demo_trex_da_04_mc_sim_ht_blocks/) | Heavy-tailed Student-t(3) AR(1)-block data, Gaussian vs Student-t noise scenarios: SNR, rho, Q, M, target false discovery rate (tFDR), and linkage sweeps | `cpp/trex_selector_methods/ad/trex_da/demo_trex_da_04_mc_sim_ht_blocks` |
| [demo_trex_da_05_mc_sim_ht_blocks_plus_white](demo_trex_da_05_mc_sim_ht_blocks_plus_white/) | Heavy-tailed blocks + heavy-tailed white augmentation, Gaussian vs heavy-tailed noise: SNR, rho, Q, M, and tFDR sweeps | `cpp/trex_selector_methods/ad/trex_da/demo_trex_da_05_mc_sim_ht_blocks_plus_white` |
| [demo_trex_da_06_mc_sim_groups](demo_trex_da_06_mc_sim_groups/) | Prior-groups method: 3-level hierarchy (groups of 10 / 50 / 250) over a multi-level latent-factor DGP with Toeplitz leaf blocks, SNR sweep. The method is selected by supplying a non-empty `prior_groups` to `trex_da_control()` | `cpp/trex_selector_methods/ad/trex_da/demo_trex_da_06_mc_sim_groups` |
| [demo_trex_da_07_mc_sim_nn](demo_trex_da_07_mc_sim_nn/) | Nearest-Neighbours / MA(kappa) banded DGP: SNR sweep (DA + base T-Rex rows), a 2D kappa x rho sweep, plus TREX+DA+NN applied to AR(1) data (method-mismatch test) as a 2D SNR x rho sweep | `cpp/trex_selector_methods/ad/trex_da/demo_trex_da_07_mc_sim_nn` |
| [demo_trex_da_08_mc_sim_equi_and_bt](demo_trex_da_08_mc_sim_equi_and_bt/) | Equicorrelated DGP (single latent factor): one 2D SNR x rho sweep with a rho grid dense near 0, one `.txt`/`.csv` pair per SNR column, DA + base T-Rex rows | `cpp/trex_selector_methods/ad/trex_da/demo_trex_da_08_mc_sim_equi_and_bt` |

Notes:

- Every MC part runs three DA solver rows — `TREX-DA[+variant]: TLARS / TAFS / TOMP` (variant tags `+AR1`, `+BT`,
  `+NN`, `+EQUI`, while the prior-groups demo 06 uses plain `TREX-DA:`) — with per-solver exchangeable tie-breaking
  (`exch_tie_alpha = 0.25` for the greedy TAFS/TOMP) and stagnation control. Demos 01, 07, and 08 additionally run the
  base (no-DA) `TREX: TLARS / TAFS / TOMP` comparison rows.
- Hierarchical agglomerative clustering (HAC) linkage labels are capitalized (`Single`, `Complete`, `Average`) and also
  appear as filename tokens in the block demos' result stems.
- The prior-groups method (demo 06) is driven by passing `prior_groups` to `trex_da_control()`. The C++ core routes on
  `prior_groups` being non-empty, so `da_method` is then irrelevant.

---

## DGPs

All scenario generators live in the single [trex_da_dgps.R](trex_da_dgps.R)
(mirroring the C++ `dgp_generators.hpp` / Python `dgp_generators.py`):

- `dgp_ar1_snr` — AR(1) Toeplitz design, SNR-calibrated noise (demos 01, 07).
- `dgp_ar1_block_white_snr` / `dgp_ar1_block_snr` — AR(1)-block covariance with
  and without white augmentation (demos 02, 03).
- `dgp_block_toeplitz_hvt_snr` / `dgp_ht_block_white_snr` — heavy-tailed
  Student-t AR(1) blocks via Gaussian scale mixture, with and without t(nu)
  white augmentation (demos 04, 05).
- `dgp_groups_toeplitz_leaf` — multi-level latent-factor design with Toeplitz
  leaf blocks (prior-groups demo 06).
- `dgp_nn_snr` — banded MA(kappa) covariance via causal convolution of a shared
  innovation sheet (demo 07).
- `dgp_equi_snr` — pure equicorrelation via a single latent factor (demo 08).

All DGPs calibrate the noise variance to a target SNR = Var(signal)/Var(noise)
and take the true support as 1-based R indices. (`dgp_bt_snr` and
`dgp_nn_block_snr` are demo-less leftovers kept for reference.)

---

## Results

Each demo writes its result tables to its own `simulation_results/data/` folder as `da_trex_mc_<scenario_tag>.{txt,csv}`
— the same stems and tidy comma-separated values (CSV) schemas as the C++ / Python suites.

---

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_01_mc_sim_ar1/demo_trex_da_01_mc_sim_ar1.R      # full run
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_01_mc_sim_ar1/demo_trex_da_01_mc_sim_ar1.R 8    # 8 worker cores
```

The optional first command-line argument sets the number of parallel worker cores
(default `TREX_NUM_WORKERS`, else 8). Scripts resolve their own directory, so they run from any working
directory. Individual parts are gated behind `RUN_PART_*` flags at the top of
each file. Each demo folder has its own `README.md`, and
[run_overnight.sh](run_overnight.sh) regenerates the whole suite's result data.
