# TRexSelector_Examples — Python Sources

## Overview

This directory contains the standalone Python examples project. It consumes the
**TRexSelectorNeo** package (1.0.0), the pybind11 binding layer over the C++
TRexSelector backend, imported as

```python
import trex_selector_neo
```

The demos mirror their C++ counterparts in the sibling [cpp/](../cpp/) tree
(and the R demos in [R/](../R/)): same demo names and folder-per-demo layout
(one subfolder per demo, each with its own `README.md`), same data-generating
processes, same control parameters, same result summaries. All demos are plain
scripts — no build step is required. Monte Carlo demos write their output to a
`simulation_results/` folder next to the demo script.

> **Indexing note**: All variable indices in the Python demos are **0-based**
> (matching C++). The R package is 1-based, so R support sets are shifted by
> one relative to the values printed here.

---

## Import Conventions

All demos follow one canonical import scheme so that every symbol's origin is
visible at its use site:

```python
import os
import sys

import numpy as np
import pandas as pd

import trex_selector_neo as trex
from trex_selector_neo.ml_methods import LpNormScaler, NormType
from trex_selector_neo.ml_methods.clustering import LinkageMethod, agglomerative_cluster
from trex_selector_neo.tsolvers.ad.lars_based import TLARS_Solver
from trex_selector_neo.utils import compute_fdp, compute_tpp, numpy_to_memmap

from trex_sim_common import run_mc_trex  # local demo helpers last
```

The rules:

1. **Selector layer via the `trex` alias.** Everything the package root exports —
   selector classes (`trex.TRexSelector`, `trex.TRexDASelector`, …), control
   parameters (`trex.TRexControlParameter`, …), and enums
   (`trex.SolverTypeForTRex`, `trex.LLoopStrategy`, `trex.DummyDistribution`, …) —
   is accessed through `import trex_selector_neo as trex`, mirroring the C++
   `namespace` qualification and R's `trexselector::`. Root symbols are never
   `from`-imported, and the alias is only introduced in files that use it.
2. **Subpackage symbols via explicit `from`-imports.** Solvers come from
   `trex_selector_neo.tsolvers.ad.{lars_based,omp_based,afs_based}`, machine learning (ML) building
   blocks from `trex_selector_neo.ml_methods` (clustering from
   `trex_selector_neo.ml_methods.clustering`), and metrics/memory-mapping
   helpers from `trex_selector_neo.utils`. Call sites use the bare names.
   Subpackages are never aliased and never accessed as `trex.utils.…`.
3. **Only the public surface.** Demos never import from the raw binding module
   `trex_selector_neo.trex_selector_methods` (it bypasses the ergonomic Python
   wrappers), never from `trex_selector_neo._core`, and never through the
   legacy `ml_methods.clustering.hierarchical` shim.
4. **Import order.** Standard library → third-party (`numpy`, `pandas`,
   `matplotlib`) → `trex_selector_neo` (root alias first, then subpackages
   alphabetically) → local demo helper modules.

---

## Prerequisites and Installation

- Python >= 3.10
- numpy >= 1.26 (installed automatically with TRexSelectorNeo)
- pandas (used by the demos to write comma-separated values (CSV) result files)
- matplotlib and plotly (used by the plotters behind each demo's
  `generate_plots.sh`)

TRexSelectorNeo is built with scikit-build-core + pybind11 from the main
repository. Install it into a virtual environment at the repository root,
because the overnight runners and the `generate_plots.sh` scripts look for
`.venv/bin/python` there:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install <path-to>/TRexSelector/Python/TRexSelectorNeo
pip install pandas matplotlib plotly
pip install -r requirements-dev.txt   # optional: mypy + stubs for the linter
```

The old package name `trex_selector` is obsolete — all demos in this tree
import `trex_selector_neo`.

---

## Running a Demo

Each demo is run directly with the Python interpreter, e.g.:

```bash
python trex_selector_methods/ad/trex_core/demo_trex_01_single_run/demo_trex_01_single_run.py
```

Notes:

- Each demo inserts its own directory into `sys.path`, so sibling helper
  modules resolve regardless of the current working directory.
- The Monte Carlo simulation demos parallelize trials with
  `concurrent.futures.ProcessPoolExecutor` and write results
  (pandas CSV + plain-text tables) to a `simulation_results/` folder
  next to the demo script.

---

## Contents

Each row below is one demo area. The selector extensions in the lower rows cover Grouped Variable Selection (GVS),
Screen-TRex with its Ordinary rule and its Bootstrap rule built on a confidence interval (CI), and sparse principal
component analysis (SPCA).

| Area | Folder | What it covers | C++ counterpart |
|---|---|---|---|
| Memory mapping | [memory_mapping/](memory_mapping/) | Disk-backed matrices via `MemoryMappedMatrix`, `numpy_to_memmap`, out-of-core workflows | [cpp/memory_mapping/](../cpp/memory_mapping/) |
| ML methods — HAC | [ml_methods/hac_clustering/](ml_methods/hac_clustering/) | Hierarchical agglomerative clustering (HAC): sample/variable clustering, locality-sensitive hashing (LSH) linkage, linkage comparisons, out-of-core mmap workflow (`demo_mlm_hac_01`, `demo_mlm_hac_02_mmap`) | [cpp/ml_methods/hac_clustering/](../cpp/ml_methods/hac_clustering/) |
| ML methods — Scalers | [ml_methods/scaler_methods/](ml_methods/scaler_methods/) | Z-score and Lp-norm column scaling, inverse transform, serialization, mmap-backed (`demo_mlm_scaler_01`) | [cpp/ml_methods/scaler_methods/](../cpp/ml_methods/scaler_methods/) |
| ML methods — PCA | [ml_methods/pca/](ml_methods/pca/) | Principal component analysis (PCA) fit, restore round-trip, out-of-sample transform (`demo_mlm_pca_01`) | [cpp/ml_methods/pca/](../cpp/ml_methods/pca/) |
| ML methods — SVD | [ml_methods/svd/](ml_methods/svd/) | Truncated singular value decomposition (SVD), direct and Gram paths (`demo_mlm_svd_01`) | [cpp/ml_methods/svd/](../cpp/ml_methods/svd/) |
| ML methods — Ridge | [ml_methods/ridge_regression/](ml_methods/ridge_regression/) | L2-penalized least squares, primal/dual, lambda sweep (`demo_mlm_ridge_01`) | [cpp/ml_methods/ridge_regression/](../cpp/ml_methods/ridge_regression/) |
| ML methods — Model selection | [ml_methods/model_selection/](ml_methods/model_selection/) | K-fold cross-validation (CV) for ridge via SVD and elastic net via cyclic coordinate descent (CCD), lambda.min / lambda.1se (`demo_mlm_ms_01`, `demo_mlm_ms_02`) | [cpp/ml_methods/model_selection/](../cpp/ml_methods/model_selection/) |
| T-Rex selector | [trex_selector_methods/ad/trex_core/](trex_selector_methods/ad/trex_core/) | Classical Terminating Random Experiments (T-Rex) selector: single runs, Monte Carlo simulations, L-loop strategies, dummy distributions, memory-mapped pipelines, scalability benchmark (Demos 01–08) | [cpp/trex_selector_methods/ad/trex_core/](../cpp/trex_selector_methods/ad/trex_core/) |
| Dependency-Aware T-Rex | [trex_selector_methods/ad/trex_da/](trex_selector_methods/ad/trex_da/) | T-Rex+DA: first-order autoregressive (AR(1)), AR(1)/heavy-tailed blocks (± white noise), prior groups, nearest-neighbor, equicorrelated + Binary-Tree (BT) (Demos 01–08) | [cpp/trex_selector_methods/ad/trex_da/](../cpp/trex_selector_methods/ad/trex_da/) |
| Grouped Variable Selection T-Rex | [trex_selector_methods/ad/trex_gvs/](trex_selector_methods/ad/trex_gvs/) | T-Rex+GVS: Hastie, scattered, mixed/negative-trap, heavy-tailed, AR(1), autoregressive moving-average (ARMA), block benchmark, with the elastic net (EN), the elastic net in augmented-lasso form (EN+AUG) and the informed elastic net (IEN) (Demos 01–08) | [cpp/trex_selector_methods/ad/trex_gvs/](../cpp/trex_selector_methods/ad/trex_gvs/) |
| Screen-TRex | [trex_selector_methods/ad/trex_screening/](trex_selector_methods/ad/trex_screening/) | Ultra-high-dimensional screening: Ordinary/Bootstrap-CI, mmap, correlated designs, biobank routing, solver backends (Demos 01–06) | [cpp/trex_selector_methods/ad/trex_screening/](../cpp/trex_selector_methods/ad/trex_screening/) |
| T-Rex Sparse PCA | [trex_selector_methods/ad/trex_spca/](trex_selector_methods/ad/trex_spca/) | T-Rex SPCA: Monte Carlo (MC) comparison vs. ordinary/oracle PCA on a sparse 3-factor model, plus Fig.-3-style proportion of explained variance (PEV) sweeps and the plugin-vs-oracle false discovery rate (FDR) / true positive rate (TPR) heatmap study (Demos 01–02) | [cpp/trex_selector_methods/ad/trex_spca/](../cpp/trex_selector_methods/ad/trex_spca/) |
| T-Solvers | [tsolvers/](tsolvers/README.md) | 12 standalone terminating-solver demos (`demo_ts_01` … `demo_ts_12`: early stopping, external normalization, serialization, mmap). The validation programs moved to the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/tsolvers/`) | [cpp/tsolvers/](../cpp/tsolvers/) |

The `tsolvers/` and `trex_selector_methods/` folders follow the library layout by dummy type. `ad/` holds the
augmented-dummy (AD) demos, while `sd/` (sparse dummies) and `vd/` (virtual dummies) hold a README until the Python
bindings expose these methods.

---

## Coverage vs. cpp

The Python examples cover the classical T-Rex selector, the Dependency-Aware
T-Rex, the Grouped-Variable-Selection T-Rex, Screen-TRex, T-Rex Sparse PCA,
the standalone terminating solvers (`tsolvers/`, 12 demos), memory mapping,
and the full `ml_methods` demo set (HAC clustering, scaler methods, PCA, SVD,
ridge regression, model selection) — every cpp demo suite is ported, with
matching demo names, folder layout, and result-file schemas (as of the
last sync with cpp).

The scaler/HAC cross-language reference programs that used to live under
`cpp/ml_methods/validation/` moved to the TRexSelector library test suite
(`TRexSelector/cpp/tests/validation/hac_clustering/`). As validation
infrastructure rather than demos, they have no Python port by design.
