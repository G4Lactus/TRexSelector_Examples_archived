# Demo: Memory-Mapped Matrix Usage

## Purpose

This demo illustrates how to use the memory-mapped matrix interface of the **TRexSelectorNeo** package from R.

It mirrors the corresponding C++ example in
[cpp/memory_mapping/demo_memory_mapping/](../../../cpp/memory_mapping/demo_memory_mapping/)
(and the Python example in [Python/memory_mapping/demo_memory_mapping/](../../../Python/memory_mapping/demo_memory_mapping/)),
replacing the C++-specific features (RAII cleanup, `Eigen::Map` interoperability, OpenMP parallel generation) with their
R equivalents: `on.exit()` cleanup, the `mmap_matrix` S3 methods (`print`, `dim`, `[`, `[<-`, `as.matrix`), and
deterministic per-column seeding.

The demo is organized into four scenarios:

1. **Basics**: create a memory-mapped matrix, write data into it, inspect metadata, extract a block, verify a full
   roundtrip, and reopen the file in read-only mode.
2. **Out-of-core streaming (serial)**: generate a matrix column by column with `writeBin()` without allocating the full
   matrix in R memory.
3. **Per-column seeded fill**: fill columns through `[<-` with deterministic per-column random seeds. The
   OpenMP-parallel version of this loop is C++-exclusive; the seeding pattern that makes it deterministic is reproduced
   here serially.
4. **Element-wise access**: test scalar reads, scalar writes, compound assignment, bounds checking, read-only
   protection, and the R block extensions (`[` / `[<-` on contiguous index ranges, `as.matrix()`).

---

## What this demo shows

The demo focuses on several practical ideas:

- how a file-backed matrix can be used like an ordinary R matrix,
- how to reduce RAM usage by writing data directly into a mapped file,
- how `on.exit()` ensures that temporary files are cleaned up safely,
- how deterministic column-wise random number generation supports reproducible parallel workflows.

Entry points used (from `TRexSelectorNeo`):

- `mem_matrix(filename, rows, cols, mode)` — open (or create) a disk-backed matrix in `"readwrite"` or `"readonly"`
  mode
- `mem_matrix_from(mat, filename)` — write an R matrix to a binary file and return an `mmap_matrix` handle
  over it
- the S3 methods `print()`, `dim()`, `[`, `[<-`, and `as.matrix()` on the returned handle

---

## Running the demo

```bash
Rscript R/memory_mapping/demo_memory_mapping/demo_memory_mapping.R
```

The demo writes its output to the console only. No `simulation_results/` folder is created for this example.

---

## What to look for

When reading the console output, check the following points:

- the created matrix has the expected dimensions,
- extracted blocks and full roundtrip checks match the original in-memory data,
- read-only reopening works as expected,
- serial out-of-core generation reproduces the expected column values,
- per-column seeds reproduce each column independently of fill order,
- invalid access triggers the intended error handling (out-of-bounds indices, writes on a read-only handle),
- temporary `.bin` backing files are removed automatically after use.

---

## Technical notes

- The backing files are created in the R session's temporary directory via `tempfile()`.
- Cleanup is managed through `on.exit()` inside each demo function; the handles go out of scope when the function
  returns, so the mapping is released before the file is removed.
- Binary files store columns as contiguous `double` blocks (column-major), matching the Eigen ColMajor storage used by
  the C++ backend.
- All indices are **1-based**, following the R language convention (the C++ and Python counterparts are 0-based).
  Per-column seeds are therefore
  `1000 + (c - 1)` and `3000 + (c - 1)` so that column 1 in R uses the same base seed as column 0 in C++.
- `[` and `[<-` accept **contiguous** index ranges only — they map onto the C++ block API; non-contiguous slices raise
  an error.
- The mapped file is handed to the C++ backend as a pointer, so solvers and selectors consume it without a copy; the
  R-level accessors (`[`, `as.matrix()`) do copy the requested block into R memory.
