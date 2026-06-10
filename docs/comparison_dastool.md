# Comparison with DAS_Tool

[DAS_Tool](https://github.com/cmks/DAS_Tool) (Sieber et al., 2018, *Nature Microbiology*) is the de-facto standard for bin refinement in **bacterial** metagenomics. `euk-bin-refiner` shares the same core idea (iteratively pick the best bin, redistribute contigs) but is tailored for **eukaryotic** datasets.

## Side-by-side

| Feature | DAS_Tool | euk-bin-refiner |
|---------|----------|-----------------|
| **Target taxa** | Prokaryotes (bacteria, archaea) | Eukaryotes (protists, microalgae, fungi…) |
| **Marker set** | 107 archaeal/bacterial single-copy genes | BUSCO `eukaryota_odb10` (255 markers) |
| **Score formula** | `score = (C − α·D)` with α dataset-dependent | `score = (C − α·D)/100`, default α = 5.0 |
| **C/D recomputation** | Length-weighted approximation | Exact: re-read BUSCO `full_table.tsv` per iteration |
| **Loser tracking** | Not exported | Yes (`bin_relationships.tsv`, `bin_summary_loss.tsv`) |
| **Dependencies** | DIAMOND, USEARCH, prodigal, R | Python ≥ 3.8, pandas |
| **External tool calls** | Yes (runs marker search internally) | No (BUSCO must be run upstream) |
| **License** | BSD-3-clause | MIT |

## Why not just run DAS_Tool on eukaryotes?

DAS_Tool's marker set (107 single-copy prokaryotic genes) is biologically meaningless for eukaryotic bins:

* Most of these genes are **bacterial-specific**.
* Eukaryotic single-copy markers are different (BUSCO, EukCC use ≥ 255 markers).
* Running DAS_Tool on a eukaryotic bin returns near-zero completeness for *every* bin, so the algorithm cannot rank them.

A practical fix would be to plug eukaryotic markers into DAS_Tool, but DAS_Tool's internal scoring code expects bacterial conventions and the score recomputation uses length-weighted scaling, which produces artifacts in our setting (see below).

## Why exact recomputation matters

In bacterial datasets, length-weighted scaling is a reasonable approximation: a 5 Mb bacterial genome with 107 markers has them roughly uniformly distributed, so losing 10% of length usually means losing ≈10% of markers.

In eukaryotic datasets this fails badly:

* Eukaryotic genomes are **larger** (10–500+ Mb) and contain large regions devoid of single-copy markers (introns, repeats, plastid/mitochondrial DNA, transposons).
* Markers tend to be **clustered** in a few protein-coding-rich contigs.
* A bin that loses 10% of length may lose 60% of its markers — or none. Length-weighted scaling can't tell.

The **phantom winner** artifact we observed in development: `metabat2_bin.379` (C = 75%, D = 35%, score = −1.0) was correctly excluded, but other chimeric bins with intermediate scores (e.g. `concoct_18` with C = 90%, D = 45%) were incorrectly retained in early iterations because length-weighted scaling underestimated their post-redistribution duplication.

`euk-bin-refiner` solves this by **directly re-reading the BUSCO full_table** every time a bin's contig set changes.

## Performance

DAS_Tool runs USEARCH and DIAMOND internally, which is fast but disk-heavy. `euk-bin-refiner` separates marker computation (BUSCO, done once upstream) from refinement (lightweight Python).

Typical timings on a 286-bin dataset:

* BUSCO on 286 bins (parallel): ~8 hours of cluster time
* DAS_Tool refinement (if applicable): ~30 minutes
* euk-bin-refiner refinement: ~10 minutes

The total wall-clock time is dominated by BUSCO, not by the refinement itself.

## When to use DAS_Tool vs euk-bin-refiner

| Use case | Recommended tool |
|----------|------------------|
| Bacterial / archaeal MAGs | DAS_Tool |
| Eukaryotic MAGs | euk-bin-refiner |
| Mixed prokaryote+eukaryote bins | Run both: DAS_Tool on prokaryotic candidates, euk-bin-refiner on eukaryotic ones |
| You also need full tracking of which winner stole from which loser | euk-bin-refiner (this is a unique feature) |

## Reference

Sieber, C.M.K., Probst, A.J., Sharrar, A. et al. **Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy.** *Nat Microbiol* 3, 836–843 (2018). https://doi.org/10.1038/s41564-018-0171-1
