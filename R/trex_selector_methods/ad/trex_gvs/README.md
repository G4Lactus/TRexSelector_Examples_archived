# T-Rex+GVS: Grouped Variable Selection — R Demonstration Suite

> **Status**: migrated to the **TRexSelectorNeo** R6 application programming interface (API)
> (`TRexGVSSelector$new(..., gvs_control = trex_gvs_control(...))$select()`,
> with `trex_fdp()` / `trex_tpp()`). No longer depends on the `TRexSelector`
> package from the Comprehensive R Archive Network (CRAN).

## Purpose

R demos and Monte Carlo (MC) simulations for the **Terminating Random Experiments (T-Rex) selector with Grouped
Variable Selection (T-Rex+GVS)** on designs where the active variables cluster into groups. The folder layout mirrors
the C++/Python suites one-to-one, with one subfolder per demo plus a single consolidated data-generating process (DGP)
file and a shared simulation-common file.

```txt
trex_gvs/
  ├── README.md
  ├── trex_gvs_dgps.R                 <- all scenario DGPs (one file)
  ├── trex_gvs_sim_common.R           <- shared sweep driver + result writers
  ├── demo_trex_gvs_01_mc_sim_hastie_en_blocks/
  │   ├── demo_trex_gvs_01_mc_sim_hastie_en_blocks.R
  │   ├── README.md
  │   └── simulation_results/
  ├── demo_trex_gvs_02_mc_sim_scattered_blocks/
  │   └── ...
  └── demo_trex_gvs_08_mc_sim_block_bench/
```

Demos 01–07 each run a single Monte Carlo **2-D signal-to-noise ratio (SNR) × correlation sweep** comparing the three
dummy-generation methods **TENET / TENET_AUG / TIENET_AUG** (demo 01 adds two fixed-operating-point method
comparisons: `lambda2_method` and `hc_linkage`). The typical fixed settings are n = 200, p = 500, a target false
discovery rate (FDR) of $\mathrm{tFDR} = 0.1$, K = 20 and 200 MC runs.

Dependencies: `TRexSelectorNeo`, `parallel`. Each demo sources the
consolidated [trex_gvs_dgps.R](trex_gvs_dgps.R) (all scenario generators in one
file, mirroring the C++ `trex_gvs_dgps.hpp`), the shared
[trex_gvs_sim_common.R](trex_gvs_sim_common.R) (2-D sweep driver, matrix/table
renderers, and the `.txt`/`.csv` result writers, mirroring the Python
`gvs_sim_common.py`), and — except demo 08 — the shared
[../../support_generators.R](../../support_generators.R) /
[../../simulation_utils.R](../../simulation_utils.R) (which provides the `.run_mc_*()`
TENET/TENET_AUG/TIENET_AUG wrappers and the `.gvs_select()` selector adapter).
Demo 08's block runner and metrics live in the demo file. Each demo writes its
result tables to `simulation_results/data/` as `.txt` + `.csv`.

---

## Demos

| # | Scenario | Sweeps | DGP (in `trex_gvs_dgps.R`) | cpp counterpart |
|---|---|---|---|---|
| [01](demo_trex_gvs_01_mc_sim_hastie_en_blocks/) | Hastie: 3 equicorrelated active groups of 50 (all active, s = 150) | 2-D SNR × rho, **lambda2_method**, **hc_linkage** | `dgp_hastie_snr` | `demo_trex_gvs_01_mc_sim_hastie_en_blocks` |
| [02](demo_trex_gvs_02_mc_sim_scattered_blocks/) | Scattered-grouped: same 3 groups, members permuted across columns | 2-D SNR × rho | `dgp_scattered_grouped_snr` | `demo_trex_gvs_02_mc_sim_scattered_blocks` |
| [03](demo_trex_gvs_03_mc_sim_mixed_blocks/) | Mixed blocks: 3 active blocks (20/50/80) + 1 trap (65), random order, noise gaps | 2-D SNR × rho (5 × 7 union grid) | `dgp_mixed_blocks_snr` | `demo_trex_gvs_03_mc_sim_mixed_blocks` |
| [04](demo_trex_gvs_04_mc_sim_neg_traps/) | Negative traps: active +Z/−Z group + two inactive +Z/−Z traps (s = 100) | 2-D SNR × rho | `dgp_neg_traps_snr` | `demo_trex_gvs_04_mc_sim_neg_traps` |
| [05](demo_trex_gvs_05_mc_sim_t3_blocks/) | Heavy-tailed blocks: demo-03 geometry with Student-t(3) factors/noise | 2-D SNR × rho | `dgp_t3_equi_blocks_snr` | `demo_trex_gvs_05_mc_sim_t3_blocks` |
| [06](demo_trex_gvs_06_mc_sim_ar1_blocks/) | First-order autoregressive (AR(1)) blocks: demo-03 geometry with within-block Cor = rho^\|j−k\| | 2-D SNR × rho | `dgp_ar1_blocks_snr` | `demo_trex_gvs_06_mc_sim_ar1_blocks` |
| [07](demo_trex_gvs_07_mc_sim_arma_blocks/) | Autoregressive moving-average (ARMA) blocks: AR(1) / moving-average (MA(3)) / ARMA(2,1) active blocks + AR(1) trap | 2-D SNR × **ar_coef** | `dgp_arma_blocks_snr` | `demo_trex_gvs_07_mc_sim_arma_blocks` |
| [08](demo_trex_gvs_08_mc_sim_block_bench/) | Block benchmark: 100 blocks × 5 with 10 active, groups discovered by hierarchical agglomerative clustering (HAC) vs oracle groups | 4 methods (M1–M4) × 2 truth scenarios over a (rho × snr) grid | `dgp_block_equicorr` | `demo_trex_gvs_08_mc_sim_block_bench` |

Each demo folder has its own `README.md`. Demo 01 exposes `run_part_*` flags
at the top to enable/disable its three parts.

### Demo 08 notes (block benchmark)

Demo 08 compares four variants — **M1** = elastic net (EN) + HAC-discovered groups, **M2** = EN + oracle groups,
**M3** = informed elastic net (IEN) + HAC, **M4** = IEN + oracle — across two within-block truth scenarios
(**WHOLE** with all 5 vars/block active as the headline scenario, and **INDIVIDUAL** with 1 active var/block as the
cautionary appendix). Oracle groups are supplied via `trex_gvs_control(groups = ...)` with contiguous 1..G labels,
while the HAC variants discover groups from the data (`corr_max = 0.5`, single linkage). It reports the coordinate
true positive rate (TPR) and FDR, the block-level FDR and a scenario-specific block metric (`blk_hit` / `full_blk`),
plus **purity** (computed from the selector's per-variable cluster labels, `$groups`), `T_stop`, and `M_found`
(number of clusters formed, from `$max_clusters`).

---

## DGPs

All scenario generators live in the single [trex_gvs_dgps.R](trex_gvs_dgps.R)
(mirroring the C++ `trex_gvs_dgps.hpp`): `dgp_hastie_snr`,
`dgp_scattered_grouped_snr`, `dgp_mixed_blocks_snr`, `dgp_neg_traps_snr`,
`dgp_t3_equi_blocks_snr`, `dgp_ar1_blocks_snr`, `dgp_arma_blocks_snr`
(+ the `.sim_arma()` helper), and `dgp_block_equicorr`. Earlier speculative,
demo-less DGPs (factor-model, hapgen, neg-corr, random/unequal blocks, the
`dgp_gvs` collection) were removed in the consolidation, and so were the
fixed-noise `generate_*()` variants of the scenarios, which no demo used.

---

## Recorded results

Result tables are written to each demo's `simulation_results/data/` folder
(`<stem>.txt` + `<stem>.csv`). Manually-recorded result summaries also live
there, e.g.
[demo_trex_gvs_01_mc_sim_hastie_en_blocks/simulation_results/data/demo_trex_gvs_01_hastie_models.md](demo_trex_gvs_01_mc_sim_hastie_en_blocks/simulation_results/data/demo_trex_gvs_01_hastie_models.md).

---

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_01_mc_sim_hastie_en_blocks/demo_trex_gvs_01_mc_sim_hastie_en_blocks.R      # all parts
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_08_mc_sim_block_bench/demo_trex_gvs_08_mc_sim_block_bench.R 8              # 8 worker cores
```

Scripts resolve their own directory, so they run from any working directory.
The optional first argument sets the number of parallel worker cores (default `TREX_NUM_WORKERS`, else 8).
These are full Monte Carlo simulations (200 trials per grid point on n = 200,
p = 500). For a quick smoke run, reduce `num_MC` / grid sizes.
