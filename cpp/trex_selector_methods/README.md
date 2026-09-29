# T-Rex Selector Methods: Demos

## Overview

This folder contains all **T-Rex selector variants** built on top of the `tsolvers` and `ml_methods` layers of the
C++ project.

Each variant addresses a different variable-selection scenario — the classical baseline, dependent/correlated designs,
grouped-variable structures, ultra-high-dimensional screening, and sparse principal component analysis (PCA) — while
sharing the same core idea of controlling the false discovery rate (FDR) with the Terminating Random Experiments
(T-Rex) framework. Demo programs illustrate usage and typical results for each variant. The validation programs,
which cross-check C++ behavior against R reference implementations, moved to the TRexSelector library test suite
(`TRexSelector/cpp/tests/validation/trex_selector_methods/ad/`).

The folders follow the layout of the library, which sorts the selectors by the way they handle the dummies. [ad/](ad/)
holds the variants with augmented dummies (AD), which append an explicit dummy matrix to the predictors, and all current
demos. [sd/](sd/README.md) is reserved for the sparse-dummy (SD) selector and [vd/](vd/README.md) for the virtual-dummy
(VD) selector.

---

## Category overview

| Category | Folder | Purpose |
|----------|--------|---------|
| **Classical T-Rex** | [ad/trex_core/](ad/trex_core/README.md) | Terminating Random Experiments Selector — baseline FDR-controlled variable selection in low- and high-dimensional Gaussian linear models |
| **Dependency-Aware T-Rex (DA-TRex)** | [ad/trex_da/](ad/trex_da/README.md) | Variable selection under correlated/dependent design structures: first-order autoregressive (AR(1)) Toeplitz, equicorrelated, banded nearest-neighbor, and hierarchical block-covariance designs |
| **Grouped Variable Selection (T-Rex+GVS)** | [ad/trex_gvs/](ad/trex_gvs/README.md) | Group-structured selection via the elastic net (EN), the elastic net in augmented-lasso form (EN-AUG) and the informed elastic net (IEN) based on hierarchical agglomerative clustering (HAC) across equicorrelated, scattered, and block-covariance patterns |
| **Screening (Screen-TRex)** | [ad/trex_screening/](ad/trex_screening/README.md) | Ultra-high-dimensional pre-screening ahead of T-Rex selection, with in-memory, memory-mapped, and biobank-scale variants |
| **Sparse PCA (T-Rex SPCA)** | [ad/trex_spca/](ad/trex_spca/README.md) | Sparse PCA via per-component T-Rex+GVS selection on plug-in principal components (PCs), benchmarked against ordinary PCA and oracle-thresholded baselines: Fig.-2/Fig.-3-style support-recovery and Definition-1 proportion of explained variance (PEV) studies, plus mechanism studies of the plug-in construction (plug-in vs. oracle response, union-support FDR per PC) |

---

## Architecture context

These components form the `trex_selector_methods` layer in the project structure:

```txt
utils -> ml_methods -> tsolvers -> trex_selector_methods -> demos
```

In particular:

- each selector variant dispatches to one or more base solvers from `tsolvers` (e.g. TLARS, TLASSO, TENET) within
  its T-loop,
- variants that rely on preprocessing reuse components from `ml_methods`, for example HAC clustering for group
  construction in `trex_gvs`, or scaling utilities shared across variants,
- `trex_screening` and the memory-mapped demos in `trex` build on the memory-mapping utilities for large design
  matrices,
- the validation programs in the library test suite check C++ selector output against trusted R reference pipelines.

---

## Building

```bash
cd TRexSelector_Examples/cpp
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"
cmake --build build/debug
```

---

## Running

Run any demo executable from the build tree. For example:

```bash
./build/debug/bin/trex_selector_methods/ad/trex_core/demo_trex_01_single_run/demo_trex_01_single_run
```

For each variant's demo lineup and required parameters, see the corresponding subfolder README (once available).

---

## Folder contents

```txt
trex_selector_methods/
  ├── README.md
  ├── CMakeLists.txt
  ├── reformat_result_txt.py    # re-render saved .txt tables from CSVs (maintenance)
  ├── ad/                       # augmented dummies (AD)
  │   ├── trex_core/
  │   ├── trex_screening/
  │   ├── trex_da/
  │   ├── trex_gvs/
  │   └── trex_spca/
  ├── sd/                       # sparse dummies (SD), README only
  └── vd/                       # virtual dummies (VD), README only
```
