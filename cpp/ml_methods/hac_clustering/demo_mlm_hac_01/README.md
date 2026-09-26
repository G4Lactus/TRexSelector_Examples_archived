# Demo: In-Memory HAC Clustering Suite

## Purpose

This demo illustrates several in-memory uses of hierarchical agglomerative clustering (HAC) that are relevant for the
TRexSelector project and, in particular, for grouped variable selection workflows.

Six subproblems run in `main()`:

1. clustering samples with Euclidean distance,
2. clustering variables with correlation-based distance,
3. accelerating correlation-based clustering with Locality-Sensitive Hashing (LSH),
4. comparing linkage strategies on structured synthetic data (exact correlation),
5. comparing linkage strategies under the pure $O(1)$ LSH correlation approximation,
6. clustering AR(1) Toeplitz block structures, where correlation decays with distance.

Throughout the demo, clustering quality is assessed by the **Adjusted Rand Index (ARI)**, which compares the recovered
cluster labels against the known ground-truth labels from the synthetic data generator. Each subproblem section below
ends with the console output of the recorded run (library 8ad3376, 8 OpenMP threads, 12-core Apple Silicon, 64 GB).

---

## Statistical model and notation

Let

- $n$ be the number of samples,
- $p$ be the number of variables,
- $K$ be the number of true clusters or blocks,
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote the data matrix,
- $z_i \in \lbrace1,\dots,K\rbrace$ the true cluster label of sample $i$,
- $b(j) \in \lbrace1,\dots,K\rbrace$ the true block label of variable $j$.

The recovered clustering is evaluated with the Adjusted Rand Index

```math
\mathrm{ARI} :=
\frac{\mathrm{Index} - \mathbb{E}\lbrack\mathrm{Index}\rbrack}
{\mathrm{MaxIndex} - \mathbb{E}\lbrack\mathrm{Index}\rbrack},
```

where $\mathrm{ARI} = 1$ indicates perfect agreement between true and predicted labels.

---

## Subproblem 1: Sample clustering

### SP1 Data model

In the sample-clustering scenario, each row belongs to one of $K$ clusters and is generated according to

```math
X_{ij} = \mu_{z_i} + \sigma_{z_i} \varepsilon_{ij},
\qquad
\varepsilon_{ij} \sim \mathcal{N}(0,1).
```

Here, $\mu_{z_i}$ is the cluster-specific mean and $\sigma_{z_i}$ is the cluster-specific scale.

### SP1 Mathematical approach

The goal is to cluster the rows $`\boldsymbol{x}_1, \dots, \boldsymbol{x}_n`$ using Euclidean distance,

```math
d(i,i') = \lVert \boldsymbol{x}_i - \boldsymbol{x}_{i'} \rVert_2 .
```

A **single-linkage** dendrogram is built with the SLINK algorithm. The dendrogram is then cut to recover $K$ clusters,
and the predicted labels are compared with the true sample labels using ARI.

### SP1 Console output

1,500 samples in 50,000 dimensions, three clusters of 500. The matrix is generated row-major and clustered through a
transposed view, so the same column-wise SLINK serves both sample and variable clustering. Under four seconds.

```text
============================================================
Demo 1: SLINK Sample Clustering (Euclidean)
============================================================
Generating Row-Major data (N=1500, P=50000)...
Creating transposed view for column-wise clustering...
Set distance policy for sample clustering...
Building Single Linkage tree (SLINK)...
✓ Tree built successfully in 3.703249 seconds.

Cutting tree to form 3 clusters...
Resulting Variable Clustering:
  Predicted Cluster 0: 500 variables
  Predicted Cluster 1: 500 variables
  Predicted Cluster 2: 500 variables

--- CLUSTERING PERFORMANCE ---
Adjusted Rand Index (ARI): 1.0000 / 1.0000

[Demo 1 Complete]
```

---

## Subproblem 2: Variable clustering with correlation distance

### SP2 Data model

In the variable-clustering scenario, variables are generated in correlated blocks.
Each variable $j$ belongs to a block $b(j)$, and the entries follow a latent-factor model of the form

```math
X_{ij} =
\mu_{b(j)} +
\sigma_{b(j)}
\left(
\alpha L_{i,b(j)} + \gamma \varepsilon_{ij}
\right)
```

where

- $L_{i,b(j)} \sim \mathcal{N}(0,1)$ is the latent block signal,
- $\varepsilon_{ij} \sim \mathcal{N}(0,1)$ is independent noise,
- $\alpha$ controls the shared within-block signal,
- $\gamma$ controls the noise contribution.

This creates blockwise collinearity among variables.

### SP2 Standardization

Before correlation-based clustering, each column is standardized by

```math
\widetilde{\boldsymbol{x}}_j =
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_2}
```

so that variables are mean-centered and scaled to unit $\ell_2$-norm.

### SP2 Mathematical approach

After standardization, variables are clustered using the correlation distance

```math
d(j,\ell) = 1 - \left| \mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell) \right|.
```

Thus, variables with either strong positive or strong negative correlation are treated as close. The demo applies SLINK
to this distance structure, builds a single-linkage dendrogram on the variables, and cuts the tree into $K$ groups. The
resulting block assignments are compared against the known true block labels.

### SP2 Console output

50,000 variables in five imbalanced blocks (20,000 down to 5,000) over 10,000 samples, exact Pearson correlation. SLINK
evaluates every pair once, so the 13.4 minutes here are the exact-distance baseline that Subproblem 3 speeds up.

```text
============================================================
Demo 2: SLINK Variable Clustering (Correlation)
============================================================
Generating Correlated Variables (N=10000, P=50000)...
Standardizing variables (Mean=0, L2=1) for correlation...
Set distance policy for variable clustering (Correlation)...
Building Single Linkage tree for Variables...
✓ Tree built successfully in 805.6827 seconds.

Cutting tree to form 5 clusters...
Resulting Variable Cluster Sizes:
  Predicted Cluster 0: 20000 variables
  Predicted Cluster 1: 10000 variables
  Predicted Cluster 2: 8000 variables
  Predicted Cluster 3: 7000 variables
  Predicted Cluster 4: 5000 variables

--- CLUSTERING PERFORMANCE ---
Adjusted Rand Index (ARI): 1.0000 / 1.0000

[Demo 2 Complete]
```

---

## Subproblem 3: Variable clustering with LSH approximation

### SP3 Data model

The data model is the same block-correlated latent-factor model as in Subproblem 2, but the focus now shifts from exact
clustering to **approximate similarity search** for large variable sets.

### SP3 Mathematical approach

This subproblem compares three correlation-based strategies:

1. **Exact correlation**, using the full correlation geometry.
2. **LSH filter**, where hashing acts as a gatekeeper before more exact comparisons.
3. **LSH approximation**, where correlation is approximated through a SimHash-style angular sketch.

Conceptually, the LSH approximation replaces expensive high-dimensional similarity comparisons by fast hash-based
surrogates:

```math
\mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell)
\approx
\widehat{\mathrm{corr}}_{\mathrm{LSH}}(\boldsymbol{x}_j, \boldsymbol{x}_\ell).
```

The purpose is to study whether the block structure is still recovered accurately when correlation distances are
approximated rather than computed exactly.

### SP3 Console output

Same block model, now 40,000 samples, which makes every exact dot product four times as expensive as in Subproblem 2.
The three strategies recover the blocks perfectly and differ only in cost: the exact baseline needs 53 minutes, the LSH
filter (Hamming gate, exact recheck of the survivors) 20 minutes, the pure $O(1)$ SimHash estimate 16 seconds.

```text
============================================================
Demo 3: SLINK Variable Clustering (LSH Benchmark)
============================================================
Generating Imbalanced Correlated Variables (N=40000, P=50000)...
Standardizing variables (Mean=0, L2=1) for correlation...

------------------------------------------------------------
[Test A] Exact Pearson Correlation (The Baseline)
✓ Tree built successfully in 3161.6318 seconds.
ARI: 1.0000 / 1.0000

------------------------------------------------------------
[Test B] LSH Filter (Gatekeeper to Exact Dot Product)
✓ Tree built successfully in 1210.5001 seconds.
ARI: 1.0000 / 1.0000

------------------------------------------------------------
[Test C] LSH Approx (Pure O(1) SimHash Angular Estimation)
✓ Tree built successfully in 15.9493 seconds.

Resulting Variable Cluster Sizes (Approx):
  Predicted Cluster 0: 20000 variables
  Predicted Cluster 1: 10000 variables
  Predicted Cluster 2: 8000 variables
  Predicted Cluster 3: 7000 variables
  Predicted Cluster 4: 5000 variables
ARI: 1.0000 / 1.0000

[Demo 3 Complete]
```

---

## Subproblem 4: Linkage-method comparison

### SP4 Data model

This comparison is again carried out on structured synthetic data with known cluster labels, so that different linkage
rules can be judged by their ARI and runtime.

### SP4 Mathematical approach

The main question is how clustering changes when the dendrogram update rule changes. The demo compares several linkage
principles, including:

- **Single linkage**, based on nearest-neighbor inter-cluster distance,
- **Complete linkage**, based on farthest-neighbor inter-cluster distance,
- **Average linkage (UPGMA)**, based on average inter-cluster distance,
- **Weighted average linkage (WPGMA)**, based on weighted average inter-cluster distance,
- **Ward linkage**, based on minimum increase in within-cluster variance,
- and additional centroid/median-style variants, which don't have the reducibility property.

In abstract form, HAC iteratively merges two clusters $A$ and $B$ minimizing a linkage criterion

```math
d(A,B),
```

where the definition of $d(A,B)$ depends on the chosen linkage rule. This subproblem studies how that choice affects
both accuracy and computational behavior on the same synthetic structure. The comparison here is run under the **exact**
Pearson-correlation distance, and it also exercises centroid- and median-style variants that lack the reducibility
property.

### SP4 Console output

Clean blocks, 1,500 samples, 50,000 variables. All eight rules recover the five blocks. The times split into three
groups: SLINK and the dispatcher's single-linkage route at 125 s; the six matrix-based rules, which first fill the
50,000 x 50,000 distance matrix (the 19,208 MB buffer) and then run the NN-chain, at 128 to 135 s; and the two rules
without the reducibility property, centroid and median, which take the generic linkage path at 259 s and 200 s. Ward
switches to squared Euclidean distance, as its update rule requires.

```text
============================================================
Demo 4: Linkage Comparison
============================================================
1. Generating Clean Correlated Variables (N=1500, P=50000)...
2. Standardizing variables (Mean=0, L2=1) for correlation...
3. Setting distance policy (Correlation)...

[TEST A] SLINK:
Time: 125.4618 seconds
ARI:  1.0000 / 1.0000

[Test B] Single Linkage (via Dispatcher):
Time: 125.4572 seconds
ARI:  1.0000 / 1.0000

[Test C] Complete Linkage (via Dispatcher -> Disk/Matrix Bound):
  -> [Backend] RAM allocated 19208 MB.

Time: 129.1597 seconds
ARI:  1.0000 / 1.0000

[Test D] UPGMA (via Dispatcher -> RAM/Geometric Bound):
  -> [Backend] RAM allocated 19208 MB.

Time: 129.4517 seconds
ARI:  1.0000 / 1.0000

[Test E] WPGMA (via Dispatcher -> Disk/Matrix Bound):
  -> [Backend] RAM allocated 19208 MB.

Time: 128.4328 seconds
ARI:  1.0000 / 1.0000

[Test F] WARD (via Dispatcher -> RAM/Geometric Bound):
 -> Note: Switching to Squared Euclidean Distance!
  -> [Backend] RAM allocated 19208 MB.

Time: 135.4874 seconds
ARI:  1.0000 / 1.0000

[Test G] Centroid Linkage (via Dispatcher -> Generic Linkage):
  -> [Backend] RAM allocated 19208 MB.

Time: 259.0649 seconds
ARI:  1.0000 / 1.0000

[Test H] Median Linkage (via Dispatcher -> Generic Linkage):
  -> [Backend] RAM allocated 19208 MB.

Time: 200.3145 seconds
ARI:  1.0000 / 1.0000

[Demo 4 Complete]
```

---

## Subproblem 5: LSH-approximate linkage comparison

### SP5 Data model

This subproblem uses the same block-correlated latent-factor model as Subproblems 2–4, but on a larger, imbalanced
configuration. Because the latent factors create dense, well-separated blocks in correlation space, the geometry is
favorable for approximate hashing.

### SP5 Mathematical approach

The comparison repeats the linkage study of Subproblem 4, but replaces the exact correlation distance with the pure
$O(1)$ SimHash-style LSH approximation

```math
d_{\mathrm{LSH}}(j,\ell)
\approx
1 - \left| \mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell) \right|.
```

Single (SLINK), complete, average (UPGMA), and weighted average (WPGMA) linkage are each run under this approximation.
The purpose is to show that, when the block geometry is well separated, the approximate distance still recovers the
correct structure at a substantially reduced per-distance cost.

### SP5 Console output

The largest case: 50,000 samples, 100,000 variables. With the $O(1)$ distance, single linkage over all pairs takes 57 s;
the three matrix-based rules fill a 100,000 x 100,000 buffer of 76,440 MB before the NN-chain runs. The buffer is
allocated lazily, so the operating system backs only the pages that are touched; in this run that still exceeded the
64 GB of the machine and put about 37 GB into swap for the duration of tests B to D, which is the price of the exact
$O(p^2)$ matrix machinery at this size. All four rules recover the blocks.

```text
============================================================
Demo 5: LSH Approx Linkage Comparison
============================================================
1. Generating Imbalanced Correlated Variables (N=50000, P=100000)...
2. Standardizing variables (Mean=0, L2=1) for correlation...
3. Setting distance policy (Pure O(1) Correlation LSH Approx)...

[TEST A] Single Linkage (via SLINK with LSH):
Time: 56.6953 seconds
ARI:  1.0000 / 1.0000

[Test B] Complete Linkage (via Dispatcher -> Disk/Matrix Bound):
  -> [Backend] RAM allocated 76440 MB.

Time: 235.1925 seconds
ARI:  1.0000 / 1.0000

[Test C] UPGMA (via Dispatcher -> RAM/Geometric Bound):
  -> [Backend] RAM allocated 76440 MB.

Time: 381.8243 seconds
ARI:  1.0000 / 1.0000

[Test D] WPGMA (via Dispatcher -> Disk/Matrix Bound):
  -> [Backend] RAM allocated 76440 MB.

Time: 384.7498 seconds
ARI:  1.0000 / 1.0000

[Demo 5 Complete]
```

---

## Subproblem 6: AR(1) Toeplitz block clustering

### SP6 Data model

Here the variables are generated with a **block-diagonal AR(1) Toeplitz** correlation structure. Each block $k$ is built
by the stationary recursion
$`\boldsymbol{x}_{j} = \rho_{k}\,\boldsymbol{x}_{j-1} + \sqrt{1 - \rho_{k}^{2}}\,\boldsymbol{\varepsilon}_{j}`$ with
$`\boldsymbol{\varepsilon}_{j}`$ standard normal, which leaves every column at unit variance. The correlation therefore
decays *exactly* with the index separation inside a block, and vanishes between blocks:

```math
\mathrm{corr}(\boldsymbol{x}_{j}, \boldsymbol{x}_{j + m}) =
\begin{cases}
\rho_{k}^{m}, & b(j) = b(j+m) = k, \\
0, & \text{otherwise},
\end{cases}
```

so a variable is strongly correlated with its immediate neighbors but nearly uncorrelated with distant variables in the
same block. Each block carries its own decay rate, $\rho_{k} \in \lbrace 0.95,\, 0.90,\, 0.85,\, 0.80,\, 0.75 \rbrace$,
so even the slowest-decaying block is down to $0.95^{100} \approx 0.006$ after a hundred steps.

### SP6 Mathematical approach

This chaining geometry is used to contrast linkage rules under the **exact** Pearson-correlation distance. Single
linkage (SLINK) follows the correlation gradient and chains each block together, whereas complete and average linkage
tend to fragment the gradually decaying blocks. LSH is deliberately not used here, because the coarse angular
quantization destroys the continuous correlation gradient. In the block-factor models of Subproblems 3 and 5 the
variables of a block share one latent direction, so a 64-bit SimHash sketch separates the blocks as well as the exact
distance does. In this temporal AR(1) structure the correlation lives in the ordering of small angles between
neighbors, which the sketch quantizes away, and LSH does not reach the results of the exact-correlation run.

### SP6 Console output

This is the one subproblem where the linkage rule decides the outcome. Single linkage follows the chain of neighbors and
recovers all five blocks. Complete linkage judges two groups by their most distant members, which in a decaying chain
are nearly uncorrelated, so it never joins a block end to end and scores ARI 0.005. UPGMA averages over the whole chain
and does a little better at 0.109, still far from the structure. The three runs cost the same 4 minutes.

```text
============================================================
Demo 6: AR(1) Block Linkage Comparison
============================================================
1. Generating AR(1) Toeplitz Variables (N=3000, P=50000)...
2. Standardizing variables (Mean=0, L2=1) for correlation...
3. Setting distance policy (Exact Pearson Correlation)...

[TEST A] Single Linkage (via SLINK):
Time: 250.4720 seconds
ARI:  1.0000 / 1.0000

[Test B] Complete Linkage (via Dispatcher -> Disk/Matrix Bound):
  -> [Backend] RAM allocated 19208 MB.

Time: 242.9061 seconds
ARI:  0.0047 / 1.0000

[Test C] UPGMA (via Dispatcher -> RAM/Geometric Bound):
  -> [Backend] RAM allocated 19208 MB.

Time: 243.8951 seconds
ARI:  0.1089 / 1.0000

[Demo 6 Complete]
```

---

## Results at a glance

| Subproblem | Size (n x p) | Distance | Rule | Time | ARI |
| --- | --- | --- | --- | --- | --- |
| 1 Sample clustering | 1,500 x 50,000 | Euclidean | SLINK | 3.7 s | 1.000 |
| 2 Variable clustering | 10,000 x 50,000 | exact correlation | SLINK | 806 s | 1.000 |
| 3 LSH benchmark | 40,000 x 50,000 | exact / LSH filter / LSH approx | SLINK | 3162 s / 1211 s / 16 s | 1.000 each |
| 4 Linkage comparison | 1,500 x 50,000 | exact correlation | 8 rules | 125 to 259 s | 1.000 each |
| 5 LSH linkage comparison | 50,000 x 100,000 | LSH approx | SLINK, complete, UPGMA, WPGMA | 57 / 235 / 382 / 385 s | 1.000 each |
| 6 AR(1) Toeplitz blocks | 3,000 x 50,000 | exact correlation | SLINK / complete / UPGMA | 250 / 243 / 244 s | 1.000 / 0.005 / 0.109 |

The whole suite took 2 h 20 min on 8 threads. The same scenarios are also implemented in
[Python](../../../../Python/ml_methods/hac_clustering/demo_mlm_hac_01/) and
[R](../../../../R/ml_methods/hac_clustering/demo_mlm_hac_01/) at smaller sizes.

---

## Running the demo

Through the folder's launch script, which rebuilds the target, keeps the console output under `logs/` and reports the
wall-clock time:

```bash
cd cpp/ml_methods/hac_clustering
THREADS=8 ./run_hac_demos.sh 01
```

Or directly, the OpenMP thread count as the optional first argument (default 6):

```bash
cpp/build/release/bin/ml_methods/hac_clustering/demo_mlm_hac_01/demo_mlm_hac_01 8
```

The demo prints to the console only. Plan for hours, not minutes: the exact-correlation baseline of Subproblem 3 alone
runs 53 minutes, and Subproblem 5 needs more than 64 GB for its distance matrix.

---

## Technical notes

- All six subproblems are enabled in `main()`; each is guarded by its own `if(true)` toggle so an individual part can be
  switched off while iterating on another.
- Subproblems 4 and 6 use the exact Pearson-correlation distance, while Subproblems 3 and 5 exercise the LSH filter
  and/or pure $O(1)$ LSH approximation.
- The matrix-based linkage rules (complete, UPGMA, WPGMA, Ward, centroid, median) fill the full $p \times p$ distance
  matrix as a `double` buffer: 19,208 MB at $p = 50{,}000$, 76,440 MB at $p = 100{,}000$. The dispatcher's `use_mmap`
  flag would place that buffer in a memory-mapped file instead; the demo leaves it at its default, RAM, in every test,
  including those whose labels read "Disk/Matrix Bound".
- The parallel work is the data generation, the standardization, the per-column pre-computation and the distance rows;
  the merge loop of every linkage rule is sequential. More threads shorten the former, not the latter: on 6 threads the
  tree builds of Subproblems 2 and 3 took 932 s, 3595 s and 1392 s, on 8 threads 806 s, 3162 s and 1211 s.
