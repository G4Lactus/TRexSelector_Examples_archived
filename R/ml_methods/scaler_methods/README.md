# Scaler Methods: Demonstration Suite

## Overview

This folder contains R examples for the **column-wise scaling and normalization** utilities of the **TRexSelectorNeo**
package (`ZScoreScaler`, `LpNormScaler`).

These preprocessing methods are used to standardize features before downstream tasks such as penalized regression, and
other numerical linear-algebra workflows. The current demo illustrates both variance-based standardization and
norm-based normalization. The C++ counterparts live in
[cpp/ml_methods/scaler_methods/](../../../cpp/ml_methods/scaler_methods/).

---

## Key concepts

Given a data matrix $\boldsymbol{X} \in \mathbb{R}^{n \times p}$, scaling is applied column by column.

Two core transformations are demonstrated:

### Z-score scaling

```math
\tilde{x}_{ij} = \frac{x_{ij} - \bar{x}_j}{s_j},
```

where $`\bar{x}_j`$ is the column mean and $s_j$ is the column scale.

### Lp-norm scaling

```math
\tilde{\boldsymbol{x}}_j =
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_p},
```

when centering is enabled.

---

## Start here

| Folder | Purpose |
| ------ | ------- |
| [demo_mlm_scaler_01/](demo_mlm_scaler_01/README.md) | Demonstrates `ZScoreScaler` and `LpNormScaler`, including inverse transforms, serialization, and memory-mapped usage |

---

## Running

```bash
Rscript R/ml_methods/scaler_methods/demo_mlm_scaler_01/demo_mlm_scaler_01.R
```

Console output only.

---

## APIs used

- `MLZScoreScaler$new(center, scale)` with `$fit()`, `$transform_inplace()`, `$inverse_transform_inplace()`,
  `$get_centers()`, `$get_scales()`, `$save()`, `$load()`.
- `MLLpNormScaler$new(norm_type, center, scale)` with the same method set.
- `mem_matrix()` for the memory-mapped subproblem; the handle is passed straight to the scaler methods, which bind
  the disk-backed buffer zero-copy (all scaler methods require `"readwrite"` mode).

---

## History note

This folder was previously named `standardization/` and held two flat scripts
(`demo_standardization_01_in_memory.R`, `demo_standardization_02_mmap.R`) that scaled an 8 x 4 toy matrix. The C++ demo
covers memory-mapped scaling as Subproblem 4 of `demo_mlm_scaler_01`, so the mmap walkthrough was folded into that demo,
the standalone mmap variant was retired, and the folder was renamed to match the C++ tree.
