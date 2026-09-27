# Demo: T-MP (Matching Pursuit) Solver

## Purpose

This demo illustrates the **T-MP (Terminating Matching Pursuit)** solver.

T-MP is the simplest pursuit-style greedy method in this group: it repeatedly selects a variable based on residual correlation and updates coefficients incrementally without performing a full orthogonal refit on the active set. The demo shows how the solver is used in five settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data,
5. terminated solution paths of one warm-started object, drawn as a figure.

---

## Key idea

Like the other terminating solvers, T-MP stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set:

```math
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TMP}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
```

In the demo, this is done through calls of the form

```cpp
executeStep(T_stop, /*early_stop=*/true)
```

so, with `early_stop=true`, the solver terminates early once the number of dummy variables in the active set reaches $T_{\text{stop}}$.

---

## Demo structure

### 1. Basic T-MP with early stopping

The first part runs T-MP with **internal normalization** and an intercept on synthetic data, in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part applies external L2 normalization to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-MP with internal normalization disabled.

This demonstrates the externally normalized workflow and allows it to be compared conceptually with the internally normalized workflow from the first part.

### 3. Serialization and warm start

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues from the checkpoint to the same final stopping level.

The data are normalized once in front (`scaler::LpNormScaler` on `X` and `D`, `y` centered) and the reference, checkpoint and resumed solvers run as readers (`.normalize = false, .intercept = false`), so the three solvers share the data without copies and the paths are bit-identical. Three writer solvers on the same data would each center and scale it in place, which moves the data by rounding and made the paths differ (the 2026-09-25 finding). Under the library's in-place contract, a second writer on a held buffer now raises an error instead.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, and equality of the recorded action path.

### 4. Memory-mapped data

The fourth part generates synthetic data directly on disk, retrieves memory-mapped views of $X$, $D$, and $y$, and runs T-MP on those mapped objects.

This shows that the solver can also be used in workflows where large matrices are stored on disk rather than held entirely in ordinary in-memory matrices.

### 5. Terminated solution paths

The fifth part repeats the quick start of the R tlars vignette. With $n = 150$, $p = 300$, five true coefficients
equal to 1 and $L = p = 300$ dummies, one T-MP object runs to $T = 1$ and is then warm-started to $T = 2$ and $T = 5$.
Each call continues the stored path, and a snapshot of the path after each call goes to `simulation_results/data/`.
The setting `snr = 5` gives a noise standard deviation of about 1, as in the vignette.

```text
=== Demo 5: T-MP Terminated Solution Paths ===
============================================================
Configuration:
  n = 150, p = 300
  num_dummies = 300
  T_stop = 5
  True support size = 5
  True support: 0 1 2 3 4
  True coefficients: 1 1 1 1 1
  SNR = 5

T_stop = 1: steps =  8, predictors =  7, dummies = 1, FDP = 0.286, TPP = 1.000
T_stop = 2: steps = 11, predictors =  9, dummies = 2, FDP = 0.444, TPP = 1.000
T_stop = 5: steps = 15, predictors = 10, dummies = 5, FDP = 0.500, TPP = 1.000

Selection:
  Predictors: {0, 1, 4, 2, 3, 72, 183, 106, 194, 37}
  Dummies: {383, 553, 399, 469, 407}
  Path: [+v0, +v1, +v4, +v2, +v3, +v72, +v183, +D83, +v106, +v194, +D253, +D99, +v37, +D169, +D107]

Path data written to:
  .../cpp/tsolvers/demo_ts_09_tmp/simulation_results/data/tmp_terminated_paths_path.csv
  .../cpp/tsolvers/demo_ts_09_tmp/simulation_results/data/tmp_terminated_paths_actions.csv
```

Matching pursuit (MP) adds the variable most correlated with the residual and updates only that variable's
coefficient. Each coefficient is therefore set once, when its variable enters, and then stays flat, and the five true
ones keep their values of 0.87 to 1.06 from $T = 1$ to $T = 5$. Without refits, MP is the cheapest solver of its
family per step. The entry order equals T-OMP up to the first dummy, which gives a false discovery proportion (FDP) of
0.286 at $T = 1$. At $T = 5$ the FDP is 0.500, and the true positive proportion (TPP) stays at 1. T-NCGMP in its
line-search variant follows the same path to the last digit. The figure shows the three snapshots. Original variables
are black, dummies red and dashed, a red vertical line marks each included dummy, and the top axis names the variable
that enters at each step.

![T-MP terminated solution paths for T = 1, 2 and 5](simulation_results/plots/tmp_terminated_paths.png)

---

## Running

From the repository root, after `cmake --preset release && cmake --build --preset release-build` in `cpp/`:

```bash
cpp/build/release/bin/tsolvers/demo_ts_09_tmp/demo_ts_09_tmp      # Parts 1-5
cpp/build/release/bin/tsolvers/demo_ts_09_tmp/demo_ts_09_tmp 5    # Part 5 only
cpp/tsolvers/demo_ts_09_tmp/generate_plots.sh                        # path figure from the Part 5 data
```

The optional arguments select the parts to run. Part 4 writes about 44 GB of temporary files, so `5` refreshes the
path data in well under a second without it. `generate_plots.sh` uses the repository's `.venv` and writes the figure
as PNG and PDF to `simulation_results/plots/`.

The demo prints its diagnostics to the console, including configuration summaries, selected variables, selection-quality measures, serialization checks, and memory-mapped execution status messages.

---

## What to look for

When reading the console output, focus on the following points:

- whether the selected variables recover the known true support reasonably well,
- how the solver behaves under both high-dimensional and low-dimensional settings,
- whether externally normalized and internally normalized workflows both execute successfully,
- whether the serialized solver state can be reloaded and continued without changing the final path,
- whether the memory-mapped workflow completes successfully and cleans up its temporary files,
- how the terminated path grows with $T$ when the same object is warm-started.

Compared with T-OMP, T-MP does not perform a full orthogonal refit on the active set, so it can converge more slowly and may revisit the same variable multiple times. It therefore serves as a useful baseline for comparison against more elaborate pursuit variants.

---

## Notes

- The synthetic data are generated with known support and coefficients, so the selected path can be compared against ground truth.
- The serialization example writes a temporary checkpoint file named `tmp_checkpoint.bin` and removes it afterward.
- The memory-mapped example writes temporary files for $X$, $D$, and $y$, runs T-MP on those mapped views, and removes those files at the end.
- The path data and figures under `simulation_results/` are committed, so the figure is visible without a run.
