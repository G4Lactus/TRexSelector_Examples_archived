# T-Solvers: Demos

## Overview

This folder contains demos for the **terminating solvers** in `tsolvers/`.

These are variable-selection algorithms designed for high-dimensional linear-model problems inside the T-Rex framework. Each solver follows a different path-following, greedy, or pursuit-style strategy, but all share the same terminating principle: instead of computing a full regularization path, they stop early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set.

The demos run each solver in a standalone way, outside the full T-Rex calibration workflow, so that solver behavior, normalization choices, and serialization can be inspected in isolation.

---

## Demo overview

| # | Solver | Folder | Strategy |
|---|--------|--------|----------|
| 01 | **T-LARS** | [demo_ts_01_tlars/](demo_ts_01_tlars/README.md) | Least Angle Regression |
| 02 | **T-LASSO** | [demo_ts_02_tlasso/](demo_ts_02_tlasso/README.md) | Lasso path |
| 03 | **T-Stepwise** | [demo_ts_03_tstepwise/](demo_ts_03_tstepwise/README.md) | Forward stepwise regression |
| 04 | **T-ENET** | [demo_ts_04_tenet/](demo_ts_04_tenet/README.md) | Elastic net |
| 05 | **T-Stagewise** | [demo_ts_05_tstagewise/](demo_ts_05_tstagewise/README.md) | Forward stagewise regression |
| 06 | **T-OMP** | [demo_ts_06_tomp/](demo_ts_06_tomp/README.md) | Orthogonal matching pursuit |
| 07 | **T-GP** | [demo_ts_07_tgp/](demo_ts_07_tgp/README.md) | Gradient pursuit |
| 08 | **T-ACGP** | [demo_ts_08_tacgp/](demo_ts_08_tacgp/README.md) | Approximate conjugate gradient pursuit |
| 09 | **T-MP** | [demo_ts_09_tmp/](demo_ts_09_tmp/README.md) | Matching pursuit |
| 10 | **T-NCGMP** | [demo_ts_10_tncgmp/](demo_ts_10_tncgmp/README.md) | Norm-corrective generalized matching pursuit |
| 11 | **T-OLS** | [demo_ts_11_tools/](demo_ts_11_tools/README.md) | Orthogonal least squares |
| 12 | **T-AFS** | [demo_ts_12_tafs/](demo_ts_12_tafs/README.md) | Adaptive forward stepwise |

The performance demos compare all solvers side by side, for example by false discovery rate (FDR) and true
positive rate (TPR), with figures drawn by `ts_plt_utils.py`:

| # | Folder | Question |
|---|--------|----------|
| P1 | [demo_ts_perf_01_fdr_tpr_vs_T/](demo_ts_perf_01_fdr_tpr_vs_T/README.md) | FDR and TPR of the terminated paths vs. $T$ (the simulation of the R tlars vignette) |

---

## Key concepts

### Terminating principle

A terminating solver runs its selection path with `executeStep(T_stop, early_stop=true)`, where $T_{\text{stop}}$ is the **number of dummy variables that must enter the active set before the solver terminates early** (it stops once the count of active dummies reaches $T_{\text{stop}}$), usually well before the full path is computed.
The actual coefficient vector is computed internally, but serves only as tool for the selection process.

### Normalization modes

The demos are organized around two common usage patterns:

- **External normalization**: the caller preprocesses $\mathbf{X}$ before passing it to the solver.
- **Internal normalization**: the solver handles normalization itself and may expose both normalized and raw coefficient views.

### Serialization

Several solver demos also illustrate Cereal-based serialization workflows such as `save(...)` and reload operations, which are relevant for checkpointing and repeated solver use inside larger pipelines.

### Terminated solution paths

Part 5 of every demo repeats the quick start of the R tlars vignette. One solver object with $n = 150$, $p = 300$,
five true coefficients equal to 1 and $L = p = 300$ dummies is warm-started $T = 1 \to 2 \to 5$, and
`ts_perf_utils.hpp` records each snapshot of its coefficient path. `ts_plt_utils.py` then draws the coefficients
against the number of included dummies, and each demo README embeds its figure. On this data several solvers share
one path to the last digit:

| Shared path | Solvers |
| ----------- | ------- |
| T-LARS | T-LARS, T-LASSO, T-Stagewise |
| T-OMP | T-OMP, T-Stepwise, T-AFS ($\rho = 1$), T-NCGMP (fully corrective variant) |
| T-MP | T-MP, T-NCGMP (line-search variant) |
| T-GP (same entry order, coefficients within 0.074) | T-GP, T-ACGP |

T-ENET and T-OLS take paths of their own. The T-AFS figure shows four values of $\rho$ from 0.1 to 1,
and the T-NCGMP figure both variants.

Over full paths an identity holds as far as the rules of the solvers agree. In
[demo_ts_perf_01](demo_ts_perf_01_fdr_tpr_vs_T/) (100 runs up to $T = 100$, precentered readers, seeded
tie-breaking):

| Solvers | Full path |
| ------- | --------- |
| T-OMP, T-AFS ($\rho = 1$), T-NCGMP (fully corrective) | every step, bit for bit |
| T-OMP, T-Stepwise | the same steps; coefficients within $1.1 \cdot 10^{-12}$ |
| T-MP, T-NCGMP (line search) | the same steps; coefficients within $2.1 \cdot 10^{-14}$ (T-NCGMP divides each step by $\lVert x_j \rVert_2^2$, which is 1 up to rounding) |
| T-LARS, T-LASSO | bit for bit up to the first T-LASSO drop (a coefficient crossing zero) |
| T-LARS, T-Stagewise | within $2.7 \cdot 10^{-13}$ up to the first T-Stagewise removal (the NNLS reduction) |

After its first drop T-LASSO takes its own path, and so does T-Stagewise after its first removal. A variable
that T-Stagewise removes keeps its coefficient (frozen), so its selection, read from the coefficient path, can
exceed the active set; `cdiagnostics::print_selection` lists such variables separately.

---

## Building

```bash
cd cpp
cmake --preset release
cmake --build --preset release-build
```

---

## Running

Run any demo from the repository root. For example:

```bash
cpp/build/release/bin/tsolvers/demo_ts_01_tlars/demo_ts_01_tlars      # Parts 1-5
cpp/build/release/bin/tsolvers/demo_ts_01_tlars/demo_ts_01_tlars 5    # Part 5 only
cpp/tsolvers/demo_ts_01_tlars/generate_plots.sh                        # path figure from the Part 5 data
```

Every demo takes the parts to run as optional arguments. Part 4 writes about 44 GB of temporary files, so `5`
refreshes the path data in well under a second without it.

---

## Folder contents

```txt
tsolvers/
  ├── README.md
  ├── CMakeLists.txt
  ├── ts_perf_utils.hpp      # shared C++ helpers (part selection, path recorder, selection)
  ├── ts_perf_solvers.hpp    # the solvers compared by the performance demos
  ├── ts_plt_utils.py        # shared matplotlib figures
  ├── demo_ts_01_tlars/
  ├── demo_ts_02_tlasso/
  ├── demo_ts_03_tstepwise/
  ├── demo_ts_04_tenet/
  ├── demo_ts_05_tstagewise/
  ├── demo_ts_06_tomp/
  ├── demo_ts_07_tgp/
  ├── demo_ts_08_tacgp/
  ├── demo_ts_09_tmp/
  ├── demo_ts_10_tncgmp/
  ├── demo_ts_11_tools/
  ├── demo_ts_12_tafs/
  └── demo_ts_perf_01_fdr_tpr_vs_T/
```
