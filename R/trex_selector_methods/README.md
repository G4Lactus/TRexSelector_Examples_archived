# T-Rex Selector Methods — R Demos

R demonstration suites for the Terminating Random Experiments (T-Rex) selector variants, mirroring
[cpp/trex_selector_methods/](../../cpp/trex_selector_methods/). Each variant
lives in its own sub-directory, and Monte Carlo (MC) demos write results to a
source-side `simulation_results/` folder next to the demo file.

The folders follow the layout of the library, which sorts the selectors by the way they handle the dummies. [ad/](ad/)
holds the variants with augmented dummies (AD), which append an explicit dummy matrix to the predictors, and all current
demos. [sd/](sd/README.md) is reserved for the sparse-dummy (SD) selector and [vd/](vd/README.md) for the virtual-dummy
(VD) selector.

---

## Contents

| Folder | Status | What it covers |
|---|---|---|
| [ad/trex_core/](ad/trex_core/README.md) | NEW — TRexSelectorNeo R6 application programming interface (API) | Classical T-Rex: single runs, MC studies (fixed/variable support), l-loop strategies, memory-mapped workflows, scalability benchmark (demos 01–07, folder-per-demo) |
| [ad/trex_da/](ad/trex_da/README.md) | NEW — TRexSelectorNeo R6 API | Dependency-Aware T-Rex (DA-TRex) via `TRexDASelector` over first-order autoregressive (AR(1)), equicorrelated, nearest-neighbor, block-diagonal, and heavy-tailed designs (demos 01–08, folder-per-demo) |
| [ad/trex_gvs/](ad/trex_gvs/README.md) | NEW — TRexSelectorNeo R6 API | Grouped Variable Selection (GVS) T-Rex via `TRexGVSSelector` with the elastic net (EN), the elastic net in augmented-lasso form (EN+AUG), and the informed elastic net (IEN) over Hastie, scattered, mixed, negative-trap, equi-block, t-distributed, AR(1), and autoregressive moving-average (ARMA) designs (demos 01–08, folder-per-demo) |
| [ad/trex_spca/](ad/trex_spca/README.md) | NEW — TRexSelectorNeo R6 API | T-Rex sparse principal component analysis (PCA) with `TRexSPCASelector` (2 solvers × 2 modes) vs. ordinary/oracle PCA on a sparse 3-factor model (demos 01-02, folder-per-demo) |
| [ad/trex_screening/](ad/trex_screening/README.md) | NEW — TRexSelectorNeo R6 API | Screen-TRex (`TRexScreeningSelector` / `TRexBiobankScreeningSelector`): in-memory and memory-mapped screening, correlated designs, biobank-scale, and solver comparison (6 demos) |

The R-side cross-check probes for trex_spca (the rdump generator and its C++/R
comparison probes) moved to the TRexSelector library test suite
(`TRexSelector/cpp/tests/validation/trex_selector_methods/ad/trex_spca/`).

---

## Shared helpers

| File | Status | Description |
|---|---|---|
| [support_generators.R](support_generators.R) | Shared | Support-index construction helpers, sourced by demos across all variants: `make_support_capped_spread()` (deterministic evenly-spaced support with capped gap), `make_support_random()` (uniformly random support, internally using `set.seed(seed + 500000L)`), `make_support_bt_one_per_block()` (one active predictor per Binary-Tree (BT) block). All return 1-based R indices and mirror the C++ helpers. |
| [simulation_utils.R](simulation_utils.R) | NEW — migrated to the TRexSelectorNeo R6 API (`TRexGVSSelector` + `trex_fdp()` / `trex_tpp()`) | Shared MC engine for the trex_gvs demos plus generic output helpers: `plot_cormat()` (Plotly heatmap), table printers, and the GVS wrappers `.run_mc_hastie()` … `.run_mc_arma_blocks()` comparing EN / EN+AUG / IEN, plus the `.gvs_select()` adapter. (The DA-TRex MC engine lives in [ad/trex_da/da_sim_common.R](ad/trex_da/da_sim_common.R), sourced by the trex_da demos directly.) |

The NEW-API trex demos additionally use their own utility file,
[ad/trex_core/trex_sim_utils.R](ad/trex_core/trex_sim_utils.R) (see the
[trex_core README](ad/trex_core/README.md)).

---

## Running

All demos run from any working directory, e.g.

```bash
Rscript R/trex_selector_methods/ad/trex_core/demo_trex_02_mc_sim_fixed_support/demo_trex_02_mc_sim_fixed_support.R
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_01_mc_sim_ar1/demo_trex_da_01_mc_sim_ar1.R
```

The gvs MC suite additionally requires `plotly` (correlation heatmaps), and the
da / gvs MC suites use `parallel` for the MC loops. The trex_spca cross-check programs (which
additionally require `glmnet`) moved to the TRexSelector library test suite
(`TRexSelector/cpp/tests/validation/trex_selector_methods/ad/trex_spca/`).
