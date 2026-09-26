# Demo 08: T-Rex+GVS Block Benchmark — HAC-Discovered vs Oracle Groups

## Purpose

Benchmark how much accurate **group discovery** (via HAC clustering) costs compared to knowing the true groups in advance ("oracle" groups). Structurally the most different demo in the suite: instead of an SNR/rho sweep for a single method set, it runs a **4-method x 2-scenario** comparison over a `(rho x snr)` grid.

## Data-generating process

`dgp_block_equicorr` — `G = 100` equal-size blocks of `block_size = 5` (`p = 500`), with `n_active_blocks = 10` chosen at random; each active block gets a random sign coefficient of magnitude `b = 3.0`. Columns follow `X_j = sqrt(rho) * Z_g + sqrt(1 - rho) * eps`, Gaussian. `n = 200`. Two truth scenarios control which variables inside an active block are active: **WHOLE** (all 5/block, `s = 50`; the headline scenario) and **INDIVIDUAL** (1/block, `s = 10`; the cautionary appendix).

Grid: `snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` with scenario-specific rho grids — `rho_grid_whole = {0.2, 0.35, 0.5, 0.65, 0.8, 0.9}` (crosses the `corr_max = 0.5` HAC discoverability boundary) and `rho_grid_individual = {0.35, 0.9}` — at `num_MC = 200` per cell. Note `corr_max = 0.5` here (vs 0.98 in demos 01–07). Seeds are staggered per grid cell (`cell_base_seed = base_seed + cell_index * num_MC`), so each `(rho, snr)` cell draws from a distinct seed band (indexed within each scenario's own rho grid).

## Methods compared

Four variants = EN/IEN x HAC-discovered/oracle groups:

- **M1 (EN-C)** — `GVS_type = "EN"`, HAC-clustered groups. EN variants use `en_solver = "TENET_AUG"`.
- **M2 (EN-O)** — `GVS_type = "EN"`, oracle groups (the DGP's true block labels, passed via `trex_gvs_control(groups = ...)` as contiguous 1..G IDs — GVS uses `groups`, not the DA-only `prior_groups`).
- **M3 (IE-C)** — `GVS_type = "IEN"`, HAC-clustered groups (solver `TIENET_AUG`).
- **M4 (IE-O)** — `GVS_type = "IEN"`, oracle groups (solver `TIENET_AUG`).

Shared controls: `tFDR = 0.1`, `K = 20`, `hc_linkage = "Single"`, `lambda2_method = "CV_1SE_CCD"`. The central comparisons are M1-vs-M2 and M3-vs-M4.

## Parts

- **Part 1** — Full benchmark over the two truth scenarios WHOLE (headline) then INDIVIDUAL (appendix), each printed as one grid table per rho value.

Metrics: coordinate TPR/FDR, `block_FDR`, a scenario-specific block metric (`blk_hit` for INDIVIDUAL, `full_blk` for WHOLE), `purity` (fraction of true blocks whose members share one cluster label, computed from the selector's `$groups`; 1.0 by construction for the oracle methods), `T_stop`, and `M_found` (number of clusters formed, from `$max_clusters`). The CSV additionally records exact-support, block-hit, full-block, null-block-activation, and lambda2 aggregates.

## Running

```bash
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_08_mc_sim_block_bench/demo_trex_gvs_08_mc_sim_block_bench.R      # full benchmark
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_08_mc_sim_block_bench/demo_trex_gvs_08_mc_sim_block_bench.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` (`gvs_block_bench_whole` and `gvs_block_bench_individual`, each as `.txt` + `.csv`); the script sources the shared `../trex_gvs_dgps.R` and `../trex_gvs_sim_common.R`.

## Interpretation

- Central comparison is M1-vs-M2 (and M3-vs-M4): if `corr_max` is well matched to the true within-block `rho`, HAC-clustered groups should track oracle-group performance closely; a persistent gap indicates the clustering threshold is splitting or merging true blocks.
- INDIVIDUAL is expected to be hardest (only 1 of 5 correlated variables per block is truly active, so a discovered group can "win" on a correlated null neighbor); WHOLE should be easiest.
- Comparing low vs high `rho` shows how within-block correlation strength helps group discovery but can hurt localizing signal within a block in INDIVIDUAL.
- EN vs IEN at matched group-source shows whether the iterative elastic net changes the HAC-vs-oracle gap itself or only the overall detection level.
