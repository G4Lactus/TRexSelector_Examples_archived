# Demo: `SVDSolver` Usage

## Purpose

This demo illustrates basic usage of the `SVDSolver` R6 class for computing a truncated singular value decomposition.

It studies three scenarios:

1. a direct SVD path in the regime $n \ge p$,
2. a wide-matrix Gram-path case where $p > 2n$,
3. an orthogonality check for the returned singular vectors.

The demo is intended as a numerical sanity check and interface demonstration rather than a performance benchmark.

It mirrors the corresponding C++ example in
[cpp/ml_methods/svd/demo_mlm_svd_01/](../../../../cpp/ml_methods/svd/demo_mlm_svd_01/), replacing Eigen matrices with
ordinary R matrices.

---

## Mathematical setup

For a matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$, the singular value decomposition is

```math
\mathbf{X} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\top.
```

If only the top $M$ singular values are retained, the solver returns a truncated decomposition

```math
\mathbf{X}_M = \mathbf{U}_M \boldsymbol{\Sigma}_M \mathbf{V}_M^\top,
```

where

- $`\mathbf{U}_M \in \mathbb{R}^{n \times M}`$,
- $`\boldsymbol{\Sigma}_M \in \mathbb{R}^{M \times M}`$,
- $`\mathbf{V}_M \in \mathbb{R}^{p \times M}`$.

This decomposition gives the best $\text{rank-}M$ approximation to $\mathbf{X}$ in the Frobenius norm sense.

---

## Subproblem 1: Direct path

The first part of the demo generates a random Gaussian matrix with

```math
n = 300, \qquad p = 80, \qquad M = 10,
```

and computes a truncated SVD in the regime $n \ge p$.

The code then prints the dimensions of $\mathbf{U}$, $\mathbf{S}$, and $\mathbf{V}$, reports the relative reconstruction
error

```math
\frac{\|\mathbf{X}_M - \mathbf{X}\|_F}{\|\mathbf{X}\|_F},
```

and lists the top $M$ singular values.

Two `stopifnot()` checks are added on top of the C++ demo: the singular values must be sorted decreasingly, and the
full-rank decomposition ($M = p$) must reproduce $\mathbf{X}$ up to floating-point error.

---

## Subproblem 2: Gram path

The second part studies a wide matrix with

```math
n = 100, \qquad p = 500, \qquad M = 8,
```

so that $p > 2n$ and the solver dispatches to the Gram path, working on the $n \times n$ matrix
$\mathbf{X}\mathbf{X}^\top$ instead of $\mathbf{X}$ itself.

As in the first scenario, the demo reports output dimensions and the relative Frobenius reconstruction error of the
truncated approximation.

---

## Subproblem 3: Orthogonality check

The third part tests whether the returned singular vectors are numerically orthonormal. It uses a random matrix with

```math
n = 200, \qquad p = 150, \qquad M = 20,
```

computes the truncated decomposition, and evaluates

```math
\|\mathbf{U}^\top \mathbf{U} - \mathbf{I}\|_F
\qquad \text{and} \qquad
\|\mathbf{V}^\top \mathbf{V} - \mathbf{I}\|_F
```

If both quantities are sufficiently small, the demo reports that $\mathbf{U}$ and $\mathbf{V}$ are orthonormal up to
numerical precision.

---

## Running

```bash
Rscript R/ml_methods/svd/demo_mlm_svd_01/demo_mlm_svd_01.R
```

The console output is organized into three sections matching the three subproblems above.

---

## What to look for

When reading the console output, focus on these points:

- whether the dimensions of $\mathbf{U}$, $\mathbf{S}$, and $\mathbf{V}$ match the requested truncated rank $M$,
- whether the singular values are reported in descending order,
- whether the truncated reconstruction error is numerically reasonable,
- whether the orthogonality errors $\|\mathbf{U}^\top \mathbf{U} - \mathbf{I}\|_F$ and
  $\|\mathbf{V}^\top \mathbf{V} - \mathbf{I}\|_F$ are close to zero.

---

## Technical notes

- `MLSVD$new()$compute(X, M)` returns a named list with elements `U`, `S`, and `V`; the path (direct, Gram, or
  randomized) is chosen internally from the shape of `X`.
- The demo labels the first scenario as the "Direct Path" for $n \ge p$ and the second as the "Gram Path" for the
  wide-matrix regime $p > 2n$, matching the C++ demo.
- The demo does not include a dedicated rank-deficient example; the three scenarios use full-rank random Gaussian
  matrices.
- The reconstruction reported by the demo is a **truncated** reconstruction using only the top $M$ singular components,
  so the error is not expected to be near machine precision unless $M$ is large enough to capture essentially all
  variation.
- RNG seeds match the C++ demo, but R's Mersenne-Twister stream differs from `std::mt19937`, so numbers agree
  qualitatively, not bitwise.
