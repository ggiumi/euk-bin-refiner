# Examples

This folder contains a minimal usage example.

## `run_example.sh`

A bash wrapper showing how to call `euk_bin_refiner.py` with typical parameters.

To use:

1. Edit the paths at the top of `run_example.sh` to point to your data.
2. Activate an environment with Python ≥ 3.8 and pandas.
3. Run:

```bash
bash run_example.sh
```

## Expected input layout

```
my_project/
├── all_bins/                    ← --bins-dir
│   ├── concoct_1.fa
│   ├── concoct_2.fa
│   ├── metabat2_bin.99.fa
│   ├── remag_bin_5.fa
│   └── ...
└── busco_results/               ← --busco-dirs
    ├── concoct_1/run_eukaryota_odb10/full_table.tsv
    ├── concoct_2/run_eukaryota_odb10/full_table.tsv
    ├── metabat2_bin.99/run_eukaryota_odb10/full_table.tsv
    ├── remag_bin_5/run_eukaryota_odb10/full_table.tsv
    └── ...
```

## Expected output layout

```
results/refined/
├── winners.tsv
├── contig_assignments.tsv
├── refinement_log.tsv
├── bin_relationships.tsv
├── bin_summary_loss.tsv
└── dereplicated_bins/
    ├── concoct_2.fa
    ├── remag_bin_5.fa
    └── ...
```
