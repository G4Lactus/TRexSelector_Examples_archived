# SVD: Demonstration Suite

## Overview

This folder contains R examples for the **singular value decomposition (SVD)** solver of the **TRexSelectorNeo**
package (`SVDSolver`).

SVD is used in this project for stable matrix decompositions, dimensionality reduction, and regression-related
computations. In particular, it is relevant for PCA and for numerically stable ridge-regression workflows.

Equivalent SVD demos are available in the C++ and Python packages; see
[cpp/ml_methods/svd/](../../../cpp/ml_methods/svd/README.md) and
[Python/ml_methods/svd/](../../../Python/ml_methods/svd/README.md) for the corresponding walkthroughs.

---

## Key concept

Given a matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$, the singular value decomposition is

```math
\mathbf{X} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\top,
```

where

- $\mathbf{U}$ contains left singular vectors,
- $\boldsymbol{\Sigma}$ contains singular values,
- $\mathbf{V}$ contains right singular vectors.

A truncated $\textrm{rank-}M$ approximation uses only the leading $M$ singular components:

```math
\mathbf{X}_M = \mathbf{U}_M \boldsymbol{\Sigma}_M \mathbf{V}_M^\top.
```

---

## Start here

| Folder | Purpose |
| ------ | ------- |
| [demo_mlm_svd_01/](demo_mlm_svd_01/README.md) | Demonstrates basic `SVDSolver` usage, truncated reconstruction, and orthogonality checks |

---

## Running

```bash
Rscript R/ml_methods/svd/demo_mlm_svd_01/demo_mlm_svd_01.R
```

Console output only.

---

## APIs used

- `SVDSolver$new()` with `$compute(X, M)`, returning a list with elements `U` ($n \times M$), `S` ($M$), and `V`
  ($p \times M$).

---

## Folder contents

```txt
svd/
  ├── README.md
  └── demo_mlm_svd_01/
      ├── demo_mlm_svd_01.R
      └── README.md
```
