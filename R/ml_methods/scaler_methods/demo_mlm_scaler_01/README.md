# Demo: `ZScoreScaler` and `LpNormScaler`

## Purpose

This demo illustrates two column-scaling utilities used for data normalization:

- `ZScoreScaler`, which standardizes columns using a mean-and-scale transformation,
- `LpNormScaler`, which centers columns and rescales them to unit $L_p$-norm.

In addition to basic in-memory usage, the demo also covers inverse-transform accuracy, serialization of fitted scaler
parameters, and scaling of data held in a memory-mapped file.

It mirrors the corresponding C++ example in
[cpp/ml_methods/scaler_methods/demo_mlm_scaler_01/](../../../../cpp/ml_methods/scaler_methods/demo_mlm_scaler_01/),
replacing Eigen matrices with R matrices and the OpenMP parallel fill of the mapped matrix with a column-wise fill.

---

## Mathematical setup

Let $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote a data matrix with columns
$`\boldsymbol{x}_1, \dots, \boldsymbol{x}_p`$.

### Z-score scaling

For each column $j$, the transformed column is

```math
\tilde{\boldsymbol{x}}_j =
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}{s_j},
```

where $`\bar{x}_j`$ is the sample mean and $s_j$ is the stored scale parameter.

### Lp-norm scaling

With centering enabled, the transformed column is

```math
\tilde{\boldsymbol{x}}_j =
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_p}.
```

The demo explicitly uses both $p = 2$ and $p = 1$ in memory, and uses the $L_2$ version again in the memory-mapped
example.

---

## Subproblem 1: Scaler comparison — center / scale switches

The first part generates a dense random matrix with

```math
n = 1000, \qquad p = 500,
```

drawn from $\mathcal{N}(5, 3^2)$ so the effect of each switch is visible, and applies

- `ZScoreScaler` under all four `center` / `scale` combinations (R `scale()` semantics),
- `LpNormScaler` with $L_2$-norm (centered and uncentered) and $L_1$-norm (centered).

Each case prints aggregate column statistics (average absolute column mean, average column SD, average column
$L_2$-norm): centering drives avg|mean| toward 0, z-score scaling drives the average SD toward 1, and $L_2$
normalization drives the average $L_2$-norm toward 1.

---

## Subproblem 2: Inverse-transform accuracy

The second part tests whether scaling can be numerically reversed after an in-place transform. It generates a random
matrix with

```math
n = 1000, \qquad p = 5000,
```

then applies `$fit()`, `$transform_inplace()`, and `$inverse_transform_inplace()` for both `ZScoreScaler` and the
$L_2$-based `LpNormScaler`.

The script reports the maximum and mean absolute reconstruction errors after the round-trip, so this subproblem checks
interface correctness and numerical reversibility rather than downstream modeling behavior.

---

## Subproblem 3: Serialization

The third part demonstrates persistence of a fitted scaler. A `ZScoreScaler` is fit on synthetic data with

```math
n = 1000, \qquad p = 500,
```

saved to disk via `$save()`, loaded back into a fresh scaler object via `$load()`, and then compared through the norms
of the differences in stored centers and scales.

This section verifies that fitted preprocessing parameters can be serialized and restored without loss, which is useful
when training and deployment happen in separate stages.

---

## Subproblem 4: Memory-mapped matrix

The fourth part applies scaling to data held in a memory-mapped file with

```math
n = 1000, \qquad p = 5000
```

(the C++ demo uses $n = 20000$, $p = 100000$; the sizes are scaled down here to keep the R demo fast). The script
creates a writable mapped matrix with `mem_matrix()`, fills it column by column with the deterministic pattern
$X_{ij} = 2i + j/100$ (0-based $i$, $j$, matching the C++ demo), fits an $L_2$-based `LpNormScaler`, transforms and
inverse-transforms the data in place on disk, and finally removes the backing file.

The `mmap_matrix` handle is passed straight to the scaler methods: the binding binds the disk-backed buffer zero-copy
through the same `Eigen::Map` interface the C++ demo uses, so no block is materialized in R memory. The scaler methods
(including `$fit()`) require the handle to be open in `"readwrite"` mode. The round-trip is verified against the
analytic fill pattern of the first column.

---

## Running

```bash
Rscript R/ml_methods/scaler_methods/demo_mlm_scaler_01/demo_mlm_scaler_01.R
```

The demo prints status messages for fitting, transforming, inverse transforming, serialization checks, and
mapped-matrix cleanup.

---

## What to look for

When reading the console output, focus on these points:

- whether the aggregate column statistics react to the `center` / `scale` switches as described above,
- whether inverse transforms return the data to the original values up to very small numerical error,
- whether serialized and reloaded `ZScoreScaler` parameters match essentially exactly,
- whether the memory-mapped example completes fit, transform, inverse transform, round-trip check, and cleanup without
  error — all directly on the mapped file, without materializing the matrix in R memory.

---

## Technical notes

- The scalers operate through Eigen on the C++ side and mutate a double-storage R matrix **in place** (zero-copy),
  returning the same buffer. A plain alias (`X2 <- X`) therefore shares storage; the demo takes explicit deep copies
  with `X + 0` before transforming so the original data survives for comparison.
- `mmap_matrix` handles are accepted by the same methods and bound zero-copy as well; because the scaler core binds a
  mutable map, all scaler methods (including `$fit()`) require the handle in `"readwrite"` mode.
- The serialization and mmap backing files are created with `tempfile()` in the session temporary directory and removed
  through `on.exit()` (the C++ demo uses the fallback filenames `scaler_test.bin` and `mmap_matrix.bin`).
- The mapped-matrix example uses `LpNormScaler` with $L_2$-norm rather than `ZScoreScaler`, matching the C++ demo.
- The C++ demo pins the OpenMP thread count to 6; the R binding does not expose a thread-count setter, so the backend
  default is used.
- RNG seeds match the C++ demo, but R's Mersenne-Twister stream differs from `std::mt19937`, so numbers agree
  qualitatively, not bitwise.
