# Demo: Ridge Regression K-Fold Cross-Validation (`RidgeCV`)

## Purpose

This demo illustrates K-fold cross-validation for the ridge penalty through the `RidgeCV` R6 class: lambda-grid
construction, the CV-MSE curve, and the `lambda.min` / `lambda.1se` selections across four scenarios.

It mirrors the corresponding C++ example in
[cpp/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/](../../../../cpp/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/),
replacing Eigen matrices with ordinary R matrices.

---

## Statistical setup

The synthetic data follow

```math
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta}_{\mathrm{true}} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma^2 \boldsymbol{I}_n),
```

with $\sigma = 1$ and a sparse ground truth: 10 coefficients equal to $1.0$ at evenly spaced column indices.

For each grid value $\lambda$, K-fold CV reports the mean squared error and its standard error across folds. Two
selections are returned:

- `lambda.min`, the grid value minimizing the CV-MSE,
- `lambda.1se`, the largest $\lambda$ whose CV-MSE lies within one standard error of the minimum — the more
  parsimonious, more heavily regularized choice.

By construction `lambda.min` $\le$ `lambda.1se`; every scenario asserts this ordering.

---

## Algorithm notes

`RidgeCV` runs an SVD per fold with per-fold centering and column $L_2$-normalization. The lambda grid is anchored at
the largest useful penalty and spans a `lambda_ratio` range (100 points over a $1000\times$ range here, matching the C++
defaults); `seed` fixes the fold assignment.

---

## Subproblem 1: Low-dimensional regime

```math
n = 1000, \qquad p = 500, \qquad K = 10 \text{ folds}.
```

The classical $n > p$ setting: the CV curve has a clear interior minimum.

---

## Subproblem 2: High-dimensional regime

```math
n = 300, \qquad p = 500, \qquad K = 10 \text{ folds}.
```

Here $p > n$, so regularization is required for the problem to be well posed.

---

## Subproblem 3: Very high-dimensional regime

```math
n = 300, \qquad p = 1000, \qquad K = 5 \text{ folds}.
```

With $p \gg n$, five folds are used so each fold's training set stays large enough relative to $n$.

---

## Subproblem 4: Memory-mapped design matrix

```math
n = 2000, \qquad p = 500, \qquad K = 10 \text{ folds}
```

(the C++ demo uses $n = 5000$; the size is scaled down here to keep the R demo fast). The design matrix is written to a
memory-mapped file with `mem_matrix_from()`, the writer handle is released, and the file is re-opened read-only
before the fit.

The `mmap_matrix` handle is passed straight to `$fit()`, which binds the disk-backed buffer zero-copy through a const
Eigen map — `"readonly"` mode suffices since the fit does not modify $\boldsymbol{X}$, exactly mirroring the C++ demo's
`ConstMapType` / `Eigen::Ref` path. Internally each fold materializes dense `X_train` / `X_test` sub-matrices in RAM,
so there is no out-of-core streaming path in either port.

---

## Running

```bash
Rscript R/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/demo_mlm_ms_01_ridge_cv_svd.R
```

Console output only; the mmap backing file is removed on exit.

---

## What to look for

When reading the console output, focus on these points:

- whether the CV-MSE curve has an interior minimum, and how sharp it is in each dimensional regime,
- how far `lambda.1se` sits from `lambda.min`, and how that gap changes with $n / p$,
- whether the ordering check (`lambda.min <= lambda.1se`) passes in every scenario,
- whether the memory-mapped scenario reproduces the same qualitative behavior as the in-memory ones.

---

## Technical notes

- APIs used: `MLRidgeCV$new()` with `$fit(X, y, num_folds, n_lambda, lambda_ratio, seed)`, `$cv_min()`, `$cv_1se()`,
  `$index_min()`, `$index_1se()`, `$get_lambdas()`, `$get_cv_errors()`, `$get_cv_std()`.
- `$index_min()` / `$index_1se()` return **1-based** indices into the lambda grid, following the R language
  convention — the binding converts deliberately at the Rcpp boundary (the C++ and Python accessors are 0-based).
- The R lambda grid is reported in ascending order; the printed curve therefore starts at the smallest penalty.
- The ordering checks are enforced with `stopifnot()` in addition to the PASS/FAIL line, so a violation aborts the
  script.
- RNG seeds match the C++ demo, but R's Mersenne-Twister stream differs from `std::mt19937`, so numbers agree
  qualitatively, not bitwise.
