# Memory-Mapped Matrix Utilities: Demonstration Suite

## Overview

This folder contains R examples for the **memory-mapped matrix** utility of the **TRexSelectorNeo** package
(`mmap_matrix()` / `convert_to_memory_mapped()`).

The utility provides a disk-backed matrix interface that behaves like an ordinary R matrix, making it useful when
matrices are too large to keep entirely in RAM. In the TRexSelector project, this is especially relevant for large
design matrices and dummy matrices that exceed available RAM in solver and selector workflows — the R object passes a
pointer to the disk-backed data straight through to the C++ backend.

The main goals of this folder are:

1. to show how memory-mapped matrices are created and accessed from R,
2. to demonstrate out-of-core workflows for large matrices,
3. to illustrate safe file lifetime management with `on.exit()`.

The C++ demo suite additionally showcases C++-specific features such as RAII and OpenMP-based parallel generation; see
[cpp/memory_mapping/](../../cpp/memory_mapping/) and [Python/memory_mapping/](../../Python/memory_mapping/) for the
corresponding walkthroughs.

---

## What this folder covers

A memory-mapped matrix stores its entries in a file on disk and exposes them through a matrix-like interface:

```math
\boldsymbol{X}_{\mathrm{mmap}} \in \mathbb{R}^{n \times p}.
```

Conceptually, the matrix behaves like a regular numeric matrix, but its storage is file-backed rather than purely
in-memory. The R binding maps `double` matrices.

This is useful for:

- **large matrices**, when RAM is limited,
- **persistent intermediate data**, when matrix contents should survive beyond a single object lifetime,
- **solver interoperability**, because the mapped data is consumed by the C++ solvers and selectors through a pointer,
  without a copy.

---

## Start here

If you are new to memory mapping in TRexSelector, begin with:

1. **`demo_memory_mapping/`** — a self-contained demo covering matrix creation, read/write access, serial out-of-core
   generation, deterministic per-column seeding, and safe cleanup.

---

## Demo overview

| Folder | Purpose |
| ------ | ------- |
| [demo_memory_mapping/](demo_memory_mapping/README.md) | Demonstrates basic mmap usage, out-of-core writes, deterministic per-column seeded filling, and element-wise access checks |

---

## Running

```bash
Rscript R/memory_mapping/demo_memory_mapping/demo_memory_mapping.R
```

---

## Folder contents

```txt
memory_mapping/
  ├── README.md
  └── demo_memory_mapping/
      ├── demo_memory_mapping.R
      └── README.md
```

---

## Notes for new users

- This demo prints results to the console and does not write simulation summary files.
- Temporary backing files are created in the R session's temporary directory and removed automatically.
- All indices are 1-based; `[` / `[<-` support contiguous index ranges only.
- The demo is intended to explain the memory-mapping utility itself, not to benchmark maximum scale production
  performance.

---

## Downstream users of the same interface

The R mmap interface is exercised further by
[../ml_methods/hac_clustering/demo_mlm_hac_02_mmap/](../ml_methods/hac_clustering/demo_mlm_hac_02_mmap/),
[../ml_methods/scaler_methods/demo_mlm_scaler_01/](../ml_methods/scaler_methods/demo_mlm_scaler_01/),
[../ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/](../ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/),
and the T-Rex demos in [../trex_selector_methods/trex/](../trex_selector_methods/trex/README.md).
