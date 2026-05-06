# Johnson–Lindenstrauss idea + sparse (Achlioptas-style) projections

This is a short theory note for the final write-up (not a formal proof).

## What we want

After preprocessing, each cell is a vector \(x \in \mathbb{R}^G\) (only HVG columns, \(G\) in the thousands). We pick a random matrix \(R \in \mathbb{R}^{G \times d}\) with \(d \ll G\) and form the sketch \(y = R^\top x\) (same as multiplying \(x\) by \(R^\top\) depending on row/column layout). For many constructions, pairwise Euclidean distances are preserved in expectation up to a small distortion \(\varepsilon\) if \(d\) is large enough on the order of \(\varepsilon^{-2} \log n\) for \(n\) cells (JL-type lemmas).

## Achlioptas / sparse coins

Achlioptas (2003) showed you can use a **sparse** random \(R\) with entries in \(\{-1,0,+1\}\) (with careful probabilities and scaling) and still get JL-style guarantees while touching fewer entries per projection—cheaper for database-style matrix multiply. sklearn’s `SparseRandomProjection` with `density="auto"` implements a standard sparse RP recipe in that spirit.

## What we actually run

We use `SparseRandomProjection` for **SRP** and `GaussianRandomProjection` for the **dense RP** baseline. Both are fit with a fixed `random_state` for reproducibility. Block-wise I/O is an **engineering** detail (out-of-core); it does not change the linear map being applied to each cell block.

## Caveats for biology

JL statements are about **Euclidean geometry** on the vectors you feed in. After log1p + HVG subset, distances are a **proxy** for biology, not a biological axiom. ARI / kNN overlap in `evaluate.py` are **empirical** checks that neighbors and clusters mostly survive compression.
