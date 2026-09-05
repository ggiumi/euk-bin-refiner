**A bin refinement algorithm for eukaryotic Metagenome-Assembled Genomes (MAGs), inspired by [DAS_Tool](https://github.com/cmks/DAS_Tool) but adapted for BUSCO eukaryotic markers.**

## Why this tool?

Recovering eukaryotic MAGs from environmental metagenomes typically requires running multiple binners in parallel (e.g. CONCOCT, MaxBin, MetaBAT, REMAG). These tools produce **overlapping bins** that share many contigs. DAS_Tool is the standard solution, but its bacterial single-copy markers do not work for eukaryotes.

`euk-bin-refiner` solves the same problem for eukaryotic data by:

1. **Using BUSCO eukaryotic markers** (`eukaryota_odb10`) as the quality signal.
2. **Exactly recomputing completeness (C) and duplication (D)** at every iteration of contig redistribution using the full_table.tsv given by BUSCO.
3. **Tracking the full history** of which losing bins gave up contigs to which winners, so chimeric / overlapping bins can be diagnosed.

The result is a **non-redundant set of high quality eukaryotic bins**, each owning a unique set of contigs.

## Quick start

### Install

```bash
git clone https://github.com/ggiumi/euk-bin-refiner.git
cd euk-bin-refiner
pip install -r requirements.txt
```

Dependencies: `python>=3.8`, `pandas>=1.3`.

### Run

```bash
python euk_bin_refiner.py \
    --bins-dir   data/all_bins/ \
    --busco-dirs data/busco_classical/ data/busco_remag/ \
    --out-dir    results/ \
    --alpha-d    5.0 \
    --min-score  0.05 \
    --min-contigs 10 \
    --write-fasta
```

A typical run on ~300 input bins with BUSCO already computed takes **5–15 minutes** on a single core (memory-bound, not CPU-bound).

## Input requirements

### Bin FASTA files

All input bin FASTA files must be in **one directory** (`--bins-dir`) with **distinguishing filename prefixes** so the tool can identify the originating binner:

```
data/all_bins/
├── concoct_1.fa
├── concoct_2.fa
├── maxbin2_bin.30.fa
├── metabat2_bin.99.fa
├── remag_bin_5.fa
└── ...
```

Recognized prefixes: `concoct`, `maxbin`, `metabat`, `remag`, `merged`. Anything else falls under `Unknown`.

### BUSCO results

For every bin, you must have a corresponding **BUSCO `eukaryota_odb10` result directory**. Pass one or more parent directories with `--busco-dirs`. Each parent directory should contain per-bin subfolders matching the bin filename:

```
data/busco_classical/
├── concoct_1/
│   └── run_eukaryota_odb10/
│       └── full_table.tsv
├── concoct_2/
│   └── run_eukaryota_odb10/
│       └── full_table.tsv
├── maxbin2_bin.30/
│   └── run_eukaryota_odb10/
│       └── full_table.tsv
└── ...
```

`full_table.tsv` is the standard BUSCO output containing per-marker contig assignments.

## Output files

All outputs land in `--out-dir/`:

| File | Description |
|------|-------------|
| `winners.tsv` | One row per refined (non-redundant) bin. Columns: bin name, originating tool, initial vs final BUSCO C and D, score before/after, contigs kept. |
| `contig_assignments.tsv` | Per-contig assignment to the winning bin. |
| `refinement_log.tsv` | Per-iteration log: which bin won, with which score, what BUSCO. |
| `bin_relationships.tsv` | Per (loser, winner) pair: how many contigs were taken, with the loser's BUSCO before/after. |
| `bin_summary_loss.tsv` | Per losing bin: how many winners took contigs from it, dominant winner, total loss. |
| `dereplicated_bins/*.fa` | Cleaned FASTA files of the refined bins (only if `--write-fasta`). |

## Algorithm

At each iteration the bin with the highest current `score = (C − α·D)/100` is selected as a winner. All its current contigs are assigned to it. Every other bin then **loses** those overlapping contigs, and its BUSCO is **recomputed from scratch** by re-evaluating the per-marker `full_table.tsv` against the contigs that remain. The loop continues until no remaining bin has `score ≥ min_score`.

The detail is the **exact BUSCO recomputation**: when a bin loses contigs, its new C and D are computed by looking at which markers' contigs are still present, not by naive proportional scaling. This avoids the *phantom winner* artifact where bins appear to retain quality they have actually lost.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--alpha-d` | 5.0 | Duplication penalty. A bin with 50% completeness and 5% duplication has score `(50 − 5·5)/100 = 0.25`. Higher α punishes chimeric bins more strongly. |
| `--min-score` | 0.05 | Bins with score below this are not eligible to win (but are still tracked, in case they donate contigs). |
| `--min-contigs` | 10 | Bins with fewer contigs are dropped. |
| `--write-fasta` | off | If set, write cleaned FASTA files of the winning bins. |

## Acknowledgements

This tool was developed during a research traineeship at the [Institut de Ciències del Mar (ICM-CSIC), Barcelona](https://www.icm.csic.es/), in the context of the PACMAN project. It was inspired by the original [DAS_Tool](https://github.com/cmks/DAS_Tool) (Sieber et al., 2018, *Nature Microbiology*).
