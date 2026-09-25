# PCA: Demonstration Suite

## Overview

This folder contains R examples for the **Principal Component Analysis (PCA)** R6 class of the **TRexSelectorNeo**
package.

The demos illustrate how a data matrix can be projected onto a lower-dimensional orthogonal basis while retaining as
much variance as possible. In the TRexSelector project, PCA can be useful as a dimensionality-reduction and
preprocessing step when working with high-dimensional data.

The main goals of this folder are:

1. to demonstrate how PCA is fit from R,
2. to show how training data and new data are transformed into principal-component scores,
3. to illustrate how explained variance is inspected after fitting,
4. to verify that temporary in-place preprocessing can be safely reversed.

Equivalent PCA demos are available in the C++ and Python packages; see
[cpp/ml_methods/pca/](../../../cpp/ml_methods/pca/README.md) and
[Python/ml_methods/pca/](../../../Python/ml_methods/pca/README.md) for the corresponding walkthroughs.

---

## What this folder covers

Given a data matrix

```math
\boldsymbol{X} \in \mathbb{R}^{n \times p},
```

PCA seeks a low-dimensional representation using $M$ principal directions with $M \le p$.

The demo in this folder focuses on:

- **basic PCA fitting**, including score and loading dimensions,
- **explained variance inspection**, for the retained components,
- **restore and transform behavior**, for both training data and new data.

---

## Start here

If you are new to this folder, begin with:

1. **`demo_mlm_pca_01/`**
   A self-contained PCA demo covering fit, restore, and transform workflows.

---

## Demo overview

| Folder | Purpose |
| ------ | ------- |
| [demo_mlm_pca_01/](demo_mlm_pca_01/README.md) | Demonstrates PCA fitting, explained variance, restore round-trip checks, and transformation of new data |

---

## Running

```bash
Rscript R/ml_methods/pca/demo_mlm_pca_01/demo_mlm_pca_01.R
```

Console output only.

---

## APIs used

- `MLPCA$new(X, M, center = TRUE, normalize = TRUE)` with `$fit()`, `$transform(X_new)`, `$restore()`, `$get_scores()`,
  `$get_loadings()`, `$get_explained_variance()`, `$get_means()`, `$get_norms()`.

---

## Folder contents

```txt
pca/
  ├── README.md
  └── demo_mlm_pca_01/
      ├── demo_mlm_pca_01.R
      └── README.md
```

---

## Notes for new users

- This folder currently contains one introductory PCA demo.
- The demo prints its diagnostics to the console.
- `MLPCA$new(X, ...)` maps the memory of `X` directly, so `$fit()` standardizes `X` **in place**; take an explicit deep
  copy with `X + 0` where the original values are needed afterwards (a plain alias shares the same buffer).
- The demo is intended to explain PCA usage and interface behavior rather than to benchmark large-scale performance.
