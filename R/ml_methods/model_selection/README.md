# Model Selection: Demonstration Suite

## Overview

This folder contains R examples for the **model selection** routines of the **TRexSelectorNeo** package (`RidgeCV`,
`ElasticNet`, `ElasticNetCV`).

The demos focus on regularized Gaussian regression models and show how hyperparameters are selected from synthetic data
using K-fold cross-validation. In the TRexSelector project, these workflows are relevant for tuning penalization
strength in ridge- and elastic-net-type estimators.

The main goals of this folder are:

1. to demonstrate model selection for regularized linear regression from R,
2. to illustrate how K-fold cross-validation is used to choose penalty parameters,
3. to compare numerically stable SVD-based ridge fitting with coordinate-descent elastic-net fitting,
4. to connect fitted penalty choices to downstream TRexSelector workflows.

Equivalent model selection demos are available in the C++ and Python packages; see
[cpp/ml_methods/model_selection/](../../../cpp/ml_methods/model_selection/README.md) and
[Python/ml_methods/model_selection/](../../../Python/ml_methods/model_selection/README.md) for the corresponding
walkthroughs.

---

## What this folder covers

The demos study penalized regression models of the form

```math
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma^2 \boldsymbol{I}_n).
```

They focus on:

- **ridge regression**, where coefficient shrinkage is controlled by a quadratic penalty,
- **elastic-net regression**, which combines $\ell_1$ and $\ell_2$ regularization,
- **cross-validation workflows**, used to select tuning parameters such as $\lambda$ and $\alpha$.

---

## Start here

If you are new to this folder, begin with:

1. **`demo_mlm_ms_01_ridge_cv_svd/`**
   A ridge-regression CV demo using the SVD-based solver across several dimensional regimes.

2. **`demo_mlm_ms_02_enet_cv_ccd/`**
   An elastic-net demo covering both regularization-path fitting and K-fold CV via coordinate descent.

---

## Demo overview

| Folder | Purpose |
| ------ | ------- |
| [demo_mlm_ms_01_ridge_cv_svd/](demo_mlm_ms_01_ridge_cv_svd/README.md) | Ridge regression K-fold CV via `RidgeCV`, including a memory-mapped design-matrix example |
| [demo_mlm_ms_02_enet_cv_ccd/](demo_mlm_ms_02_enet_cv_ccd/README.md) | Elastic-net path fitting (`ElasticNet`) and K-fold CV (`ElasticNetCV`) via coordinate descent, including ridge/lasso/elastic-net comparisons |

---

## Running

```bash
Rscript R/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/demo_mlm_ms_01_ridge_cv_svd.R
Rscript R/ml_methods/model_selection/demo_mlm_ms_02_enet_cv_ccd/demo_mlm_ms_02_enet_cv_ccd.R
```

Console output only.

---

## Folder contents

```txt
model_selection/
  ├── README.md
  ├── demo_mlm_ms_01_ridge_cv_svd/
  │   ├── demo_mlm_ms_01_ridge_cv_svd.R
  │   └── README.md
  └── demo_mlm_ms_02_enet_cv_ccd/
      ├── demo_mlm_ms_02_enet_cv_ccd.R
      └── README.md
```

---

## Notes for new users

- Both demos assert their invariants with `stopifnot()` (ordering of `lambda.min` / `lambda.1se`, monotone deviance
  ratio, `$fit_grid()` reproduction), so a violated expectation aborts the script rather than printing quietly.
- Index convention: `$index_min()` / `$index_1se()` are **1-based**, following the R language convention (the
  binding converts deliberately at the Rcpp boundary); the C++ and Python ports are 0-based.
- The glmnet cross-check that used to sit here (`demo_en_glmnet_compare.R` + `rdump_en/`) is local-testing material, not
  a demo. It moved to the TRexSelector library test suite
  (`TRexSelector/cpp/tests/validation/ml_methods/model_selection/`).
