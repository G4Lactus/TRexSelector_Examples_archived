# ML Methods: Demos

## Overview

This folder contains the R demos for the `ml_methods` layer of the TRexSelector library, exercised through the
**TRexSelectorNeo** R6 API.

It provides the numerical and preprocessing components that support higher-level T-Rex workflows, including column
scaling, singular value decomposition, principal component analysis, ridge regression, model parameter selection based
on k-fold cross-validation, and hierarchical agglomerative clustering.

Each method family has its own subfolder. Demo programs illustrate usage and interface behavior. Cross-language
correctness checks against R references live in the TRexSelector library's own test suite (`cpp/tests/validation/`),
run on demand via its `run_validations` target, rather than in this examples package.

The demos mirror the C++ suite in [cpp/ml_methods/](../../cpp/ml_methods/): same sub-demo structure, parameters, and
console-output layout, with Eigen matrices replaced by ordinary R matrices. Where the C++ demos run at very large scale
(e.g. OpenMP-parallel generation of matrices with $p \ge 50{,}000$), the R mirrors scale the problem sizes down — the
reductions are noted in each script's header comment.

---

## Category overview

| Category | Folder | Purpose |
| -------- | ------ | ------- |
| **Scaler Methods** | [scaler_methods/](scaler_methods/README.md) | Column normalization utilities, including Z-score and Lp-norm scaling |
| **PCA** | [pca/](pca/README.md) | Principal component analysis |
| **SVD** | [svd/](svd/README.md) | Singular value decomposition for stable matrix computations |
| **Ridge Regression** | [ridge_regression/](ridge_regression/README.md) | L2-penalized least-squares regression |
| **Model Selection** | [model_selection/](model_selection/README.md) | Cross-validation workflows for regularized regression methods |
| **HAC Clustering** | [hac_clustering/](hac_clustering/README.md) | Hierarchical agglomerative clustering, including memory-mapped workflows |

---

## Architecture context

These components form the `ml_methods` layer in the project structure:

```txt
utils -> ml_methods -> tsolvers -> trex_selector_methods -> demos
```

In particular:

- scaler methods support feature normalization before downstream modeling,
- SVD and ridge-related tools support stable regression and dimensionality-reduction workflows,
- HAC clustering supports grouped-variable workflows such as GVS.

---

## Running

Each demo is a standalone script run directly with `Rscript` (no build step required), from any working directory:

```bash
Rscript R/ml_methods/pca/demo_mlm_pca_01/demo_mlm_pca_01.R
```

All demos assume an installed **TRexSelectorNeo** (v0.2.0); see the [R README](../README.md) for installation.

---

## Folder contents

```txt
ml_methods/
  ├── README.md
  ├── scaler_methods/
  ├── pca/
  ├── svd/
  ├── ridge_regression/
  ├── model_selection/
  └── hac_clustering/
```

---

## Notes on the R port

- All demos print to the console only; no `simulation_results/` folders are produced in this tree.
- Index convention: everything user-facing is **1-based** in R (`ml_cut_tree()` labels, `$index_min()` /
  `$index_1se()`), following the R language convention — the binding converts deliberately at the Rcpp boundary. The
  C++ and Python ports are 0-based.
- RNG seeds match the C++ demos, but R's Mersenne-Twister stream differs from `std::mt19937`, so printed numbers agree
  qualitatively, not bitwise.
- Memory-mapped inputs are first-class: `mmap_matrix` handles are accepted directly by `ml_agglomerative_cluster()`, the
  scalers, and `RidgeCV`, and are bound zero-copy through the same Eigen map interface the C++ demos use (the scalers
  require `"readwrite"` handles; `RidgeCV` also accepts `"readonly"`).
- The `scaler_methods/` folder was previously named `standardization/`; it was renamed and its two flat scripts were
  folded into `demo_mlm_scaler_01/` to match the C++ tree.
- The C++ validation programs that used to live under `cpp/ml_methods/validation/` moved to the TRexSelector library
  test suite (`TRexSelector/cpp/tests/validation/hac_clustering/` and `.../ml_methods/model_selection/`); they are
  cross-language checks by design, not demos.
