# Demo: Elastic-Net Path Fitting and K-Fold CV (`ElasticNet` / `ElasticNetCV`)

## Purpose

This demo illustrates glmnet-style pathwise cyclic coordinate descent (CCD) through two R6 classes:

- **Part A** fits standalone regularization paths with `ElasticNet`,
- **Part B** selects the penalty by K-fold cross-validation with `ElasticNetCV`, including the pure-ridge
  `lambda_2_lars` conversion used internally by TRexGVS.

It mirrors the corresponding C++ example in
[cpp/ml_methods/model_selection/demo_mlm_ms_02_enet_cv_ccd/](../../../../cpp/ml_methods/model_selection/demo_mlm_ms_02_enet_cv_ccd/),
replacing Eigen matrices with ordinary R matrices.

---

## Statistical setup

The synthetic data follow

```math
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta}_{\mathrm{true}} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma^2 \boldsymbol{I}_n),
```

with $\sigma = 1$ and 10 nonzero coefficients equal to $1.0$ at evenly spaced column indices.

The elastic net solves, for each $\lambda$ on the path,

```math
\widehat{\boldsymbol{\beta}}(\lambda) =
\arg\min_{\boldsymbol{\beta}}
\frac{1}{2n} \lVert \boldsymbol{y} - \boldsymbol{X}\boldsymbol{\beta} \rVert_2^2
+ \lambda \left( \alpha \lVert \boldsymbol{\beta} \rVert_1
+ \tfrac{1 - \alpha}{2} \lVert \boldsymbol{\beta} \rVert_2^2 \right),
```

where the mixing parameter interpolates between the lasso ($\alpha = 1$, sparse fits) and ridge ($\alpha = 0$, dense
shrinkage).

---

## Part A1: Path fit, low-dimensional

```math
n = 300, \qquad p = 100,
```

with $\alpha \in \lbrace 0.0, 0.5, 1.0 \rbrace$. For each $\alpha$ the demo prints a condensed path summary — every
tenth $\lambda$ with the number of non-zero coefficients, the $L_1$-norm of the coefficient vector, and the deviance
ratio (glmnet's `dev.ratio`, rising as $\lambda$ decreases).

Two checks are added on top of the C++ demo: the deviance ratio must be non-decreasing along the path, and
`$fit_grid()` evaluated at the auto grid must reproduce the auto path exactly.

---

## Part A2: Path fit, high-dimensional

```math
n = 200, \qquad p = 500,
```

with $\alpha \in \lbrace 0.0, 1.0 \rbrace$ — the $p > n$ regime, where the lasso path stops adding variables once the
active set saturates.

---

## Part B1: K-fold CV, low-dimensional

```math
n = 300, \qquad p = 100, \qquad K = 10 \text{ folds},
```

for $\alpha \in \lbrace 0.0, 0.5, 1.0 \rbrace$. Each fit prints the CV-MSE curve with the `lambda.min` and
`lambda.1se` positions marked, and asserts the ordering `lambda.min <= lambda.1se`.

---

## Part B2: Ridge CV and the `lambda_2_lars` conversion

```math
n = 300, \qquad p = 200, \qquad K = 10 \text{ folds}, \qquad \alpha = 0,
```

i.e. pure ridge as used inside TRexGVS. The demo converts both selections to the T-Rex scale via

```math
\lambda_{2,\mathrm{lars}} = \lambda_{\lbrace \min, 1se \rbrace} \cdot \frac{p}{2},
```

mirroring R's `lm_dummy` and `cv.glmnet(alpha = 0)`. This is the value TRexGVS uses internally via `CV_1SE_CCD`.

---

## Running

```bash
Rscript R/ml_methods/model_selection/demo_mlm_ms_02_enet_cv_ccd/demo_mlm_ms_02_enet_cv_ccd.R
```

Console output only.

---

## What to look for

When reading the console output, focus on these points:

- how the non-zero count grows along the path as $\lambda$ decreases, and how it differs between $\alpha = 1$ (sparse)
  and $\alpha = 0$ (dense),
- whether the deviance ratio increases monotonically along each path,
- whether `$fit_grid()` reproduces the auto-grid path exactly,
- where `lambda.min` and `lambda.1se` land on the CV curve for each $\alpha$,
- the reported `lambda_2_lars` values, which is the quantity consumed by TRexGVS.

---

## Technical notes

- APIs used: `MLElasticNet$new()` with `$fit(X, y, alpha)` / `$fit_grid(X, y, lambda_grid, alpha)`, `$get_coef()`,
  `$get_intercepts()`, `$get_lambdas()`, `$get_dev_ratio()`, `$converged()`, `$predict(X_new)`; and
  `MLElasticNetCV$new()` with `$fit(X, y, alpha, n_folds, n_lambda, lambda_min_ratio, seed, standardize, intercept,
  max_iter, tol)`, `$cv_min()`, `$cv_1se()`, `$index_min()`, `$index_1se()`, `$get_lambdas()`, `$get_cv_errors()`,
  `$get_cv_std()`.
- `$index_min()` / `$index_1se()` return **1-based** indices into the lambda grid, following the R language
  convention — the binding converts deliberately at the Rcpp boundary (the C++ and Python accessors are 0-based).
- `$get_coef()` returns the coefficient path in the original predictor scale as a $p \times K$ matrix; the lambda grid
  is descending (glmnet convention) and may hold fewer than `n_lambda` points when glmnet's `fdev`/`devmax` early
  termination kicks in.
- `ElasticNetCV` builds the full-data lambda grid once and reuses it across folds (glmnet semantics).
- RNG seeds match the C++ demo, but R's Mersenne-Twister stream differs from `std::mt19937`, so numbers agree
  qualitatively, not bitwise.
