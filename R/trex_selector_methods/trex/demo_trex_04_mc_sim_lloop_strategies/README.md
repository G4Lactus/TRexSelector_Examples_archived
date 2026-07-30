# Demo 04: L-Loop Strategy Comparison

## Purpose

Compare the T-Rex **L-loop strategies** — the mechanisms by which the dummy-variable
structure is generated/expanded across L-loop iterations — repeated across three
base solvers (**TLARS**, **TOMP**, **TAFS** with `rho_afs = 0.3`). Shows how the
strategy affects FDR, TPR, and the calibrated dummy multiplier `L` / stopping
time `T`.

## Data-generating process

`dgp_gauss_snr` — i.i.d. Gaussian design, `y = X beta + eps` calibrated to SNR.
High-dimensional config: `n = 300`, `p = 1000`, `s = 10`, `tFDR = 0.1`, fixed
coefficients `beta_j = 1`. The support is a random size-`s` subset redrawn per trial
via `make_support_random(s, p, seed = trial_seed)` (the `block_support = TRUE`
contiguous-block variant is available but its call is commented out). SNR grid
`seq(0.1, 2.0, by = 0.1)` plus `5.0` (21 values), `num_MC = 200` per strategy x SNR.

## Strategies compared

`L_STRATEGIES` — **8 rows** (matches the 8-row reference study):

- STANDARD — fresh i.i.d. dummy matrix each L-loop iteration.
- HCONCAT — horizontally concatenated dummy columns.
- PERMUTATION — re-use the stored base dummy matrix via deterministic row
  permutations per experiment.
- PERMUTATION_SEEDED — seed-derived base + row permutations, nothing stored
  (bit-identical to PERMUTATION for the same seed).
- SEEDED — seed-derived i.i.d. draws, nothing stored.
- SKIPL_5p / SKIPL_10p / SKIPL_20p — skip the L-loop, fixed
  `L = max_dummy_multiplier * p` at 5p / 10p / 20p.

(The former `DIRECT` / `PERMUTATION_DIRECT` strategy names were renamed to
`SEEDED` / `PERMUTATION_SEEDED`; the old names are no longer accepted. The
C++/Python demos still carry the removed `ONDEMAND` / `PERMUTATION_ONDEMAND`
labels and are stale against the current library — this R file is the corrected
version.)

The solver sweep (`SOLVERS`) is the outer loop; one `.txt` + `.csv` file pair is
written per solver. Shared control per solver x strategy: `K = 20`,
`use_max_T_stop = TRUE`, `opt_threshold = 0.75`, `tloop_max_stagnant_steps = 5`,
`use_openmp = FALSE`, `use_memory_mapping = FALSE`, `lloop_strategy` and
`max_dummy_multiplier` swept; `tloop_stagnation_stop` is left unset (AUTO — the
per-solver-family default). Reports averaged FDR, TPR, Avg L, Avg T per strategy
x SNR (`track_L_T = TRUE`). MC trials run in parallel over a hardcoded
`num_cores <- 6L`.

## Running

```bash
Rscript R/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies/demo_trex_04_mc_sim_lloop_strategies.R
```

Results are written to this demo's own `simulation_results/data/` folder as a `.txt` table
and a tidy `.csv` (the "solver" column holds the strategy name), stems
`demo_trex_04_lloop_strategies_results_n300_p1000_{tlars,tomp,tafs}_random_support`.
The script sources `../../support_generators.R`, `../../simulation_utils.R`, and
`../trex_sim_utils.R`.
