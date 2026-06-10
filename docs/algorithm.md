# Algorithm details

## Overview

`euk-bin-refiner` implements an iterative greedy algorithm that selects a set of non-redundant eukaryotic bins from overlapping bins produced by multiple metagenomic binners. The algorithm is conceptually identical to [DAS_Tool](https://github.com/cmks/DAS_Tool) (Sieber et al., 2018), but uses eukaryotic BUSCO markers and an *exact* re-evaluation of completeness at every step.

## Scoring

Each bin receives a score based on its BUSCO statistics:

```
score = (C − α·D) / 100
```

where:

* `C` = % of BUSCO eukaryotic markers found as Complete (single-copy or duplicated)
* `D` = % of BUSCO markers found as Duplicated only
* `α` = duplication penalty weight (default 5.0)

The intuition: a bin with high `C` and low `D` is *complete and clean*. A bin with high `C` but also high `D` is likely a **chimera** — two or more organisms mixed together, producing duplicated markers. The `α·D` penalty therefore reject chimeric bins.

### Example

A bin with `C = 80%, D = 20%`:

* with α = 1: score = (80 − 20)/100 = 0.60 (good)
* with α = 5: score = (80 − 100)/100 = −0.20 (rejected as chimera)

## Iterative selection

```
Initialize:
    For each bin i, compute C_i, D_i, score_i from its initial contigs.
    Define candidate_pool = {i | score_i ≥ min_score}.
    Define tracking_pool  = all bins (used for monitoring losses).

Iterate:
    1. Pick winner w = argmax(score) in candidate_pool.
    2. The contigs currently assigned to w are frozen as its final contig set.
    3. For every other bin j (in tracking_pool):
        - Remove contigs that are now claimed by w.
        - Exactly recompute C_j, D_j, score_j.
        - Drop j from candidate_pool if score_j < min_score or contigs < min_contigs.
    4. Remove w from both pools.

Stop when candidate_pool is empty.
```

This guarantees that:

* No two output bins share any contig.
* The winning bin at every iteration is locally optimal for the current state.
* The final solution is a strong heuristic for the (NP-hard) problem of partitioning contigs to maximize per-bin BUSCO quality.

## Exact BUSCO recomputation

The most important technical detail is how `C_j` and `D_j` are recomputed when bin `j` loses some contigs.

### The naive approach (used by some tools)

Scale the score by the fraction of length retained:

```
C_new = C_old × (length_kept / length_initial)
```

This is **fast but wrong**. If a bin loses 10% of its length but those 10% happened to contain 90% of the BUSCO markers, the bin is effectively gutted yet its score barely drops. This produces *phantom winners*: bins that appear to retain high quality after losing most of their actual signal.

### Our exact approach

At every iteration, the BUSCO `full_table.tsv` is re-evaluated against the bin's *current* contigs:

```
For each BUSCO marker:
    Look at the contig(s) listed in full_table.tsv.
    Keep only the rows whose contig is still in the bin.
    Re-classify the marker as:
        Complete (single-copy) — 1 Complete hit remaining
        Complete (duplicated)  — ≥2 Complete hits remaining
        Fragmented             — only Fragmented hits remaining
        Missing                — no non-missing entries remaining
```

This costs more compute time (we re-iterate over all markers every step), but is the only way to get the right answer. Empirically the entire pipeline still completes in 5–15 minutes on real metagenomic datasets with hundreds of bins.

## Why tracking excluded bins matters

Bins with `score < min_score` (typically chimeric ones, e.g. `C = 75%, D = 35%, score = −1.00`) **cannot win** but they may still share contigs with winners. Without tracking these losses, you lose all visibility on chimeric bins.

In `euk-bin-refiner` such bins remain in `tracking_pool`, so every contig redistribution is recorded in `bin_relationships.tsv`. This lets you answer questions like:

* "Bin X was a chimera with `D = 35%`. Where did its contigs end up?"
* "Did winners Y and Z absorb the actual organisms hidden inside the chimera?"

A typical pattern: a chimeric `MetaBAT` bin with `C = 75%, D = 35%` is excluded, but its contigs get split across 2-3 winners (e.g. one CONCOCT bin and two REMAG bins), each representing the underlying real organisms more cleanly.

## Reproducibility

The algorithm is fully deterministic given:

* The set of input bins
* The BUSCO full_tables
* `α`, `min_score`, `min_contigs`

Tiebreaks (when two bins have identical score) are resolved by alphabetical bin name, so the result is reproducible across runs and machines.
