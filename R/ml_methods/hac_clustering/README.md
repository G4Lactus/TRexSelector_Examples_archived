# HAC Clustering: Demonstration Suite

## Overview

This folder contains R examples for **hierarchical agglomerative clustering (HAC)** through the **TRexSelectorNeo**
function `agglomerative_cluster()` and its companion `cut_tree()`.

Within the TRexSelector project, HAC is used by the **GVS (Grouped Variable Selection)** workflow to construct variable
groups from a dendrogram cut. The examples in this folder show how clustering behaves in smaller in-memory settings and
how the same ideas can be scaled to much larger problems through memory-mapped I/O.

The main goals of this folder are:

1. to demonstrate how HAC is applied to synthetic sample and variable clustering tasks,
2. to compare clustering behavior under different distance and linkage choices,
3. to illustrate how large clustering problems can be handled with memory-mapped storage,
4. to provide a bridge between exploratory clustering demos and downstream grouped-variable workflows.

Equivalent HAC clustering demos are available in the C++ and Python packages; see
[cpp/ml_methods/hac_clustering/](../../../cpp/ml_methods/hac_clustering/README.md) and
[Python/ml_methods/hac_clustering/](../../../Python/ml_methods/hac_clustering/README.md) for the corresponding
walkthroughs.

---

## What this folder covers

Hierarchical agglomerative clustering builds a tree of successive merges and then recovers cluster assignments by
cutting that tree at a chosen level.

In this folder, the demos focus on:

- **in-memory HAC experiments**, for understanding clustering structure and linkage behavior,
- **correlation-based variable clustering**, which is relevant for grouped variable selection,
- **memory-mapped large-scale HAC**, for problems that would otherwise exceed available RAM.

---

## Start here

If you are new to this folder, begin with:

1. **`demo_mlm_hac_01/`**
   An in-memory HAC demo suite covering six clustering scenarios and linkage comparisons.

2. **`demo_mlm_hac_02_mmap/`**
   A memory-mapped demo that reproduces the LSH linkage comparison at an out-of-core scale.

---

## Demo overview

| Folder | Purpose |
| ------ | ------- |
| [demo_mlm_hac_01/](demo_mlm_hac_01/README.md) | In-memory HAC demo suite for sample clustering, variable clustering, LSH benchmarks, linkage comparisons, and AR(1) block chains |
| [demo_mlm_hac_02_mmap/](demo_mlm_hac_02_mmap/README.md) | Memory-mapped HAC demo for LSH-based linkage comparison without full in-RAM allocation |

---

## Running

```bash
Rscript R/ml_methods/hac_clustering/demo_mlm_hac_01/demo_mlm_hac_01.R
Rscript R/ml_methods/hac_clustering/demo_mlm_hac_02_mmap/demo_mlm_hac_02_mmap.R
```

Console output only.

---

## APIs used

- `agglomerative_cluster(data, method, metric, use_mmap)`, where `data` is a numeric matrix or an `mmap_matrix`, and
  the enums are `LinkageMethod$...` (Ward, Average, Complete, Single, WPGMA, Median, Centroid) and
  `DistanceMetric$...` (Euclidean, Correlation, Manhattan, Correlation_LSH_Filter, Correlation_LSH_Approx).
- `cut_tree(linkage, num_orig_objs, num_clusters)` to flatten the dendrogram.
- `convert_to_memory_mapped()` / `mmap_matrix()` for the out-of-core demo.

---

## Folder contents

```txt
hac_clustering/
  ├── README.md
  ├── demo_mlm_hac_01/
  │   ├── demo_mlm_hac_01.R
  │   └── README.md
  └── demo_mlm_hac_02_mmap/
      ├── demo_mlm_hac_02_mmap.R
      └── README.md
```

---

## Notes for new users

- `agglomerative_cluster()` clusters the **rows** of its input: pass `X` for sample clustering and `t(X)` for variable
  clustering.
- Ward, median, and centroid linkage require Euclidean geometry and are rejected for correlation or Manhattan metrics.
  The single exception is Ward + `Correlation_LSH_Approx`, which runs in a 64-dimensional SimHash-projected Euclidean
  space.
- Cluster ids are **1-based**, following the R language convention — the binding converts deliberately at the Rcpp
  boundary (the C++ and Python ports are 0-based).
- Both demos are self-contained: the helper functions (ARI, the data generators, the correlation standardizer) are
  inlined in each script, mirroring the C++ and Python demos. The former shared `hac_helpers.R` was retired.
