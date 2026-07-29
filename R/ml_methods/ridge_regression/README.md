# Ridge Regression: Demonstration Suite

## Overview

This folder contains R examples for the **ridge regression** solver of the **TRexSelectorNeo** package
(`RidgeSolver`).

Ridge regression adds an $\ell_2$-penalty to least squares in order to stabilize estimation and shrink coefficients
toward zero. The demo in this folder illustrates ridge fitting in both low-dimensional and high-dimensional settings,
together with the effect of varying the penalty parameter.

Equivalent ridge regression demos are available in the C++ and Python packages; see
[cpp/ml_methods/ridge_regression/](../../../cpp/ml_methods/ridge_regression/README.md) and
[Python/ml_methods/ridge_regression/](../../../Python/ml_methods/ridge_regression/README.md) for the corresponding
walkthroughs.

---

## Key concept

Given a design matrix $\boldsymbol{X} \in \mathbb{R}^{n \times p}$, response vector $\boldsymbol{y} \in \mathbb{R}^n$,
and penalty parameter $\lambda \ge 0$, ridge regression estimates

```math
\widehat{\boldsymbol{\beta}}_{\mathrm{ridge}} =
\arg\min_{\boldsymbol{\beta} \in \mathbb{R}^p}
\lVert \boldsymbol{y} - \boldsymbol{X}\boldsymbol{\beta} \rVert_2^2
+
\lambda \lVert \boldsymbol{\beta} \rVert_2^2 .
```

---

## Start here

| Folder | Purpose |
| ------ | ------- |
| [demo_mlm_ridge_01/](demo_mlm_ridge_01/README.md) | Demonstrates ridge regression in primal and high-dimensional regimes, plus a lambda sweep |

---

## Running

```bash
Rscript R/ml_methods/ridge_regression/demo_mlm_ridge_01/demo_mlm_ridge_01.R
```

Console output only.

---

## APIs used

- `RidgeSolver$new()` with `$solve(X, y, lambda)`, returning the coefficient vector of length $p$. The primal
  ($n \ge p$) and dual ($n < p$) paths are dispatched internally.

For the cross-validated choice of $\lambda$, see
[../model_selection/demo_mlm_ms_01_ridge_cv_svd/](../model_selection/demo_mlm_ms_01_ridge_cv_svd/README.md).

---

## Folder contents

```txt
ridge_regression/
  ├── README.md
  └── demo_mlm_ridge_01/
      ├── demo_mlm_ridge_01.R
      └── README.md
```
