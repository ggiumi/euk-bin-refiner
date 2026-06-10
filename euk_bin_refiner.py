#!/usr/bin/env python3
"""
euk-bin-refiner: A DAS_Tool-equivalent bin refinement algorithm for eukaryotic MAGs.

Given multiple metagenomic bin sets produced by different binners (e.g. CONCOCT,
MaxBin, MetaBAT, REMAG), this tool selects a non-redundant subset of bins by
maximizing eukaryotic BUSCO completeness while penalizing duplication. Contig
assignments are resolved iteratively so that no contig is assigned to more than
one final bin.

Unlike DAS_Tool (which uses single-copy bacterial markers), this tool uses BUSCO
eukaryotic markers and *exactly* recomputes completeness/duplication at every
iteration of contig redistribution — avoiding the "phantom winner" artifact of
naive score scaling.

USAGE EXAMPLE:
    python euk_bin_refiner.py \\
        --bins-dir   data/all_bins/ \\
        --busco-dirs data/busco_classical/ data/busco_remag/ \\
        --out-dir    results/refined_bins/ \\
        --alpha-d    5.0 \\
        --min-score  0.05 \\
        --min-contigs 10 \\
        --write-fasta

OUTPUTS:
    winners.tsv               One row per winning (non-redundant) bin.
    contig_assignments.tsv    Per-contig assignment to its winning bin.
    refinement_log.tsv        Iteration log of the greedy selection.
    bin_relationships.tsv     Per (loser, winner) pair: how many contigs stolen.
    bin_summary_loss.tsv      Per loser bin: total contigs lost and to whom.
    dereplicated_bins/        FASTA files of refined bins (if --write-fasta).
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd


# ============================================================================
# Constants
# ============================================================================

# Default BUSCO eukaryota_odb10 marker count (used if not parsed from full_table)
DEFAULT_N_MARKERS = 255

# Tool detection patterns: maps regex pattern -> tool name
TOOL_PATTERNS = [
    (re.compile(r"^remag_"), "REMAG"),
    (re.compile(r"^concoct"), "CONCOCT"),
    (re.compile(r"^maxbin"), "MaxBin"),
    (re.compile(r"^metabat"), "MetaBAT"),
    (re.compile(r"^merged"), "Merged"),
]


# ============================================================================
# I/O helpers
# ============================================================================

def read_contigs_from_fa(fa_path: Path) -> Set[str]:
    """Extract contig IDs from a FASTA file (everything after '>' up to whitespace)."""
    contigs: Set[str] = set()
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs.add(line[1:].split()[0])
    return contigs


def load_busco_full_table(table_path: Path) -> Tuple[List[Tuple[str, str, Optional[str]]], Optional[int]]:
    """Parse a BUSCO full_table.tsv. Returns (marker rows, total markers).

    Each marker row is a tuple (marker_id, status, contig_or_None) where status
    is one of: Complete, Duplicated, Fragmented, Missing.
    """
    rows: List[Tuple[str, str, Optional[str]]] = []
    n_total: Optional[int] = None
    with open(table_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                if "number of BUSCOs" in line:
                    match = re.search(r"number of BUSCOs:\s*(\d+)", line)
                    if match:
                        n_total = int(match.group(1))
                continue
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            marker_id = parts[0]
            status = parts[1]
            contig = parts[2] if len(parts) > 2 and status != "Missing" else None
            rows.append((marker_id, status, contig))
    return rows, n_total


def detect_tool(genome_name: str) -> str:
    """Identify the binner from the bin filename based on naming conventions."""
    name_lower = genome_name.lower()
    for pattern, tool in TOOL_PATTERNS:
        if pattern.match(name_lower):
            return tool
    return "Unknown"


# ============================================================================
# BUSCO score computation
# ============================================================================

def compute_busco_C_D(
    markers: List[Tuple[str, str, Optional[str]]],
    current_contigs: Set[str],
    n_total_markers: int,
) -> Tuple[float, float, int]:
    """Compute BUSCO Complete% and Duplicated% for a bin given its current contigs.

    Crucially, this re-evaluates BUSCO from scratch given which contigs are
    currently assigned to the bin (the "phantom winner" fix: naive scaling
    would give wrong numbers when a bin loses contigs).

    A marker is:
      - Complete (single-copy) if exactly one Complete (or Duplicated) hit is on a kept contig
      - Duplicated if two or more Complete/Duplicated hits are on kept contigs
      - Fragmented if only Fragmented hits remain on kept contigs
      - Missing otherwise

    Args:
        markers: Output of load_busco_full_table()
        current_contigs: Set of contig IDs currently owned by the bin
        n_total_markers: Total marker count (denominator for C, D)

    Returns:
        (C_percent, D_percent, n_complete_markers_absolute)
    """
    by_marker: Dict[str, List[Tuple[str, Optional[str]]]] = defaultdict(list)
    for marker_id, status, contig in markers:
        by_marker[marker_id].append((status, contig))

    n_complete_single = 0
    n_complete_dup = 0

    for marker_id, entries in by_marker.items():
        # Keep only entries whose contig is still in this bin
        # (Missing entries are kept too, as they always apply)
        kept = [(s, c) for (s, c) in entries
                if (c is None) or (c in current_contigs)]
        kept_non_missing = [(s, c) for (s, c) in kept if s != "Missing"]

        if not kept_non_missing:
            continue  # Marker is Missing now

        statuses = [s for (s, _) in kept_non_missing]
        if "Complete" in statuses:
            n_comp = statuses.count("Complete")
            if n_comp >= 2:
                n_complete_dup += 1
            else:
                n_complete_single += 1
        elif "Duplicated" in statuses:
            n_dup = statuses.count("Duplicated")
            if n_dup >= 2:
                n_complete_dup += 1
            else:
                n_complete_single += 1
        # else: only Fragmented or only Missing -> not counted in C

    n_complete_total = n_complete_single + n_complete_dup
    C = 100 * n_complete_total / n_total_markers
    D = 100 * n_complete_dup / n_total_markers
    return C, D, n_complete_total


def compute_score(C: float, D: float, alpha_d: float) -> float:
    """Compute the refinement score: (C - alpha_d * D) / 100.

    Higher = better. alpha_d weighs the penalty for duplicated markers
    (chimeric bins). Default alpha_d=5 gives a strong penalty.
    """
    return (C - alpha_d * D) / 100


# ============================================================================
# Bin loading
# ============================================================================

def load_bins_with_busco(
    bins_dir: Path,
    busco_dirs: List[Path],
    min_contigs: int,
    logger: logging.Logger,
) -> Tuple[Dict[str, Set[str]],
           Dict[str, List[Tuple[str, str, Optional[str]]]],
           Dict[str, int]]:
    """Load FASTA bins and pair them with their BUSCO results.

    Returns:
        bin_contigs:  bin_filename -> set of contig IDs
        bin_markers:  bin_filename -> BUSCO marker rows
        bin_n_total:  bin_filename -> number of BUSCO markers (denominator)
    """
    # ---- Load FASTA contigs ----
    logger.info("Loading FASTA contigs from %s", bins_dir)
    bin_contigs: Dict[str, Set[str]] = {}
    for fa_path in sorted(bins_dir.glob("*.fa")):
        contigs = read_contigs_from_fa(fa_path)
        if len(contigs) >= min_contigs:
            bin_contigs[fa_path.name] = contigs
    logger.info("  Loaded %d bins with >= %d contigs", len(bin_contigs), min_contigs)

    # ---- Load BUSCO ----
    logger.info("Loading BUSCO full_tables from %d directories", len(busco_dirs))
    bin_markers: Dict[str, List[Tuple[str, str, Optional[str]]]] = {}
    bin_n_total: Dict[str, int] = {}

    for busco_dir in busco_dirs:
        if not busco_dir.is_dir():
            logger.warning("  Skipping non-existent BUSCO dir: %s", busco_dir)
            continue
        for sub in sorted(busco_dir.iterdir()):
            if not sub.is_dir():
                continue
            bin_name_noext = sub.name
            if bin_name_noext.endswith("_busco"):
                bin_name_noext = bin_name_noext[:-len("_busco")]
            candidates = [f"{bin_name_noext}.fa", f"remag_{bin_name_noext}.fa"]

            full_table = sub / "run_eukaryota_odb10" / "full_table.tsv"
            if not full_table.is_file():
                continue
            markers, n_total = load_busco_full_table(full_table)
            if n_total is None:
                n_total = DEFAULT_N_MARKERS

            for fa_name in candidates:
                if fa_name in bin_contigs:
                    bin_markers[fa_name] = markers
                    bin_n_total[fa_name] = n_total
                    break

    logger.info("  BUSCO loaded for %d / %d bins", len(bin_markers), len(bin_contigs))

    # Drop bins with no BUSCO data
    missing = set(bin_contigs) - set(bin_markers)
    if missing:
        logger.warning("  %d bins without BUSCO data will be excluded", len(missing))
        for m in missing:
            bin_contigs.pop(m, None)

    return bin_contigs, bin_markers, bin_n_total


# ============================================================================
# Greedy refinement
# ============================================================================

def greedy_refinement(
    bin_contigs_initial: Dict[str, Set[str]],
    bin_markers: Dict[str, List[Tuple[str, str, Optional[str]]]],
    bin_n_total: Dict[str, int],
    alpha_d: float,
    min_score: float,
    min_contigs: int,
    logger: logging.Logger,
) -> Tuple[List[Tuple[str, Set[str]]],
           Dict[str, str],
           List[dict],
           List[dict],
           Dict[str, int]]:
    """Run the greedy iterative refinement.

    At each iteration:
      1. Pick the bin with the highest current score (>= min_score).
      2. Assign all its current contigs to it (it's a 'winner').
      3. Remove those contigs from every other bin.
      4. For each affected bin, exactly recompute C, D, and score.
      5. Drop bins that fall below min_score or min_contigs.
      6. CRUCIAL: even bins excluded from candidacy (initial score < min_score)
         are tracked, to record which winners 'stole' their contigs.

    Returns:
        final_winners:    list of (bin, contigs)
        contig_owner:     contig_id -> winning bin
        log_rows:         per-iteration log
        relationships:    per (loser, winner, iteration) event
        winner_iteration: bin -> iteration when it won
    """
    # Compute initial scores
    logger.info("Computing initial BUSCO and scores for %d bins",
                len(bin_contigs_initial))
    initial_busco: Dict[str, Tuple[float, float]] = {}

    # current_score/current_contigs: only candidates that can still win
    current_score: Dict[str, float] = {}
    current_contigs: Dict[str, Set[str]] = {}

    # tracking_contigs: ALL bins, including those with initial score < min_score.
    # Used so we don't miss tracking when a low-quality bin loses contigs to a winner.
    tracking_contigs: Dict[str, Set[str]] = {}

    for b, ctgs in bin_contigs_initial.items():
        C, D, _ = compute_busco_C_D(bin_markers[b], ctgs, bin_n_total[b])
        initial_busco[b] = (C, D)
        score = compute_score(C, D, alpha_d)
        tracking_contigs[b] = set(ctgs)
        # Only enter the candidate pool if score is good enough
        if score >= min_score:
            current_score[b] = score
            current_contigs[b] = set(ctgs)

    initial_n_contigs = {b: len(ctgs) for b, ctgs in bin_contigs_initial.items()}

    # Iterative selection
    logger.info("Starting greedy selection (min_score=%.3f, min_contigs=%d)",
                min_score, min_contigs)
    final_winners: List[Tuple[str, Set[str]]] = []
    contig_owner: Dict[str, str] = {}
    log_rows: List[dict] = []
    relationships: List[dict] = []
    winner_iteration: Dict[str, int] = {}

    iteration = 0
    while True:
        candidates = {b: s for b, s in current_score.items()
                      if s >= min_score and len(current_contigs[b]) >= min_contigs}
        if not candidates:
            break

        # Deterministic tiebreak: highest score, then alphabetical bin name
        best = max(candidates.keys(), key=lambda b: (candidates[b], b))
        best_score = candidates[best]
        best_contigs = current_contigs[best].copy()
        best_C_init, best_D_init = initial_busco[best]
        best_C_now, best_D_now, _ = compute_busco_C_D(
            bin_markers[best], best_contigs, bin_n_total[best]
        )

        # Record the win
        final_winners.append((best, best_contigs))
        winner_iteration[best] = iteration
        for c in best_contigs:
            contig_owner[c] = best

        log_rows.append({
            "iteration": iteration,
            "bin": best,
            "tool": detect_tool(best),
            "score_at_selection": round(best_score, 4),
            "C_initial": round(best_C_init, 1),
            "D_initial": round(best_D_init, 1),
            "C_at_selection": round(best_C_now, 1),
            "D_at_selection": round(best_D_now, 1),
            "n_contigs_owned": len(best_contigs),
            "n_contigs_initial": initial_n_contigs[best],
            "pct_kept": round(100 * len(best_contigs) / initial_n_contigs[best], 1),
        })

        # Steal contigs from every other bin (CANDIDATE or NON-CANDIDATE)
        # We iterate over tracking_contigs to also capture low-quality bins
        for other in list(tracking_contigs.keys()):
            if other == best:
                continue

            other_ctgs = tracking_contigs[other]
            lost = other_ctgs & best_contigs
            if not lost:
                continue

            other_ctgs -= lost
            tracking_contigs[other] = other_ctgs

            # Recompute BUSCO after the loss
            new_C, new_D, _ = compute_busco_C_D(
                bin_markers[other], other_ctgs, bin_n_total[other]
            )
            new_score = compute_score(new_C, new_D, alpha_d)

            other_C_init, other_D_init = initial_busco[other]

            # Record the relationship
            was_candidate = other in current_score
            relationships.append({
                "loser_bin": other,
                "loser_tool": detect_tool(other),
                "loser_was_candidate": was_candidate,
                "winner_bin": best,
                "winner_tool": detect_tool(best),
                "winner_iteration": iteration,
                "n_contigs_taken": len(lost),
                "n_contigs_loser_initial": initial_n_contigs[other],
                "pct_of_loser_contigs": round(100 * len(lost) / initial_n_contigs[other], 2),
                "loser_C_initial": round(other_C_init, 1),
                "loser_D_initial": round(other_D_init, 1),
                "loser_C_after_this_steal": round(new_C, 1),
                "loser_D_after_this_steal": round(new_D, 1),
                "loser_score_after_this_steal": round(new_score, 4),
                "loser_contigs_remaining_after": len(other_ctgs),
            })

            # Update candidacy
            if other in current_score:
                current_contigs[other] = other_ctgs
                if len(other_ctgs) < min_contigs or new_score < min_score:
                    del current_score[other]
                    del current_contigs[other]
                else:
                    current_score[other] = new_score

        # Remove the winner from any further consideration
        current_score.pop(best, None)
        current_contigs.pop(best, None)
        tracking_contigs.pop(best, None)

        iteration += 1
        if iteration % 20 == 0:
            logger.info("  Iteration %d: %d winners, %d candidates left",
                        iteration, len(final_winners), len(current_score))

    logger.info("Greedy finished: %d iterations, %d final winners",
                iteration, len(final_winners))

    return final_winners, contig_owner, log_rows, relationships, winner_iteration


# ============================================================================
# Output writers
# ============================================================================

def write_outputs(
    out_dir: Path,
    bin_contigs_initial: Dict[str, Set[str]],
    bin_markers: Dict[str, List[Tuple[str, str, Optional[str]]]],
    bin_n_total: Dict[str, int],
    initial_busco: Dict[str, Tuple[float, float]],
    initial_n_contigs: Dict[str, int],
    alpha_d: float,
    final_winners: List[Tuple[str, Set[str]]],
    contig_owner: Dict[str, str],
    log_rows: List[dict],
    relationships: List[dict],
    winner_iteration: Dict[str, int],
    bins_dir: Path,
    write_fasta: bool,
    logger: logging.Logger,
) -> None:
    """Write all output TSVs and (optionally) cleaned FASTA bins."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- winners.tsv ----
    winners_rows = []
    for b, ctgs in final_winners:
        C_init, D_init = initial_busco[b]
        C_final, D_final, _ = compute_busco_C_D(
            bin_markers[b], ctgs, bin_n_total[b]
        )
        winners_rows.append({
            "bin": b,
            "tool": detect_tool(b),
            "BUSCO_C_initial": round(C_init, 1),
            "BUSCO_D_initial": round(D_init, 1),
            "BUSCO_C_final": round(C_final, 1),
            "BUSCO_D_final": round(D_final, 1),
            "delta_C": round(C_final - C_init, 1),
            "delta_D": round(D_final - D_init, 1),
            "initial_score": round(compute_score(C_init, D_init, alpha_d), 4),
            "final_score": round(compute_score(C_final, D_final, alpha_d), 4),
            "n_contigs_initial": initial_n_contigs[b],
            "n_contigs_final": len(ctgs),
            "pct_kept": round(100 * len(ctgs) / initial_n_contigs[b], 1),
        })
    winners_df = pd.DataFrame(winners_rows).sort_values(
        "BUSCO_C_final", ascending=False
    )
    winners_df.to_csv(out_dir / "winners.tsv", sep="\t", index=False)
    logger.info("  Wrote: winners.tsv (%d rows)", len(winners_df))

    # ---- contig_assignments.tsv ----
    pd.DataFrame(
        [{"contig": c, "winner_bin": b} for c, b in contig_owner.items()]
    ).to_csv(out_dir / "contig_assignments.tsv", sep="\t", index=False)
    logger.info("  Wrote: contig_assignments.tsv (%d rows)", len(contig_owner))

    # ---- refinement_log.tsv ----
    pd.DataFrame(log_rows).to_csv(
        out_dir / "refinement_log.tsv", sep="\t", index=False
    )

    # ---- bin_relationships.tsv ----
    rel_df = pd.DataFrame(relationships).sort_values(
        ["loser_bin", "winner_iteration"]
    )
    rel_df.to_csv(out_dir / "bin_relationships.tsv", sep="\t", index=False)
    logger.info("  Wrote: bin_relationships.tsv (%d rows)", len(rel_df))

    # ---- bin_summary_loss.tsv ----
    winners_set = {b for b, _ in final_winners}
    summary_rows = []
    for b in bin_contigs_initial:
        bin_rels = rel_df[rel_df["loser_bin"] == b]
        if len(bin_rels) == 0:
            continue

        total_taken = int(bin_rels["n_contigs_taken"].sum())
        n_winners_took = len(bin_rels)
        n_init = initial_n_contigs[b]
        c_init, d_init = initial_busco[b]
        was_candidate = bool(bin_rels["loser_was_candidate"].iloc[0])

        # Status
        if b in winners_set:
            status = f"won_at_iter_{winner_iteration[b]}"
        elif was_candidate:
            status = "excluded_during_greedy"
        else:
            status = "excluded_low_initial_score"

        # Final C/D
        if b in winners_set:
            owned = next(ctgs for (bin_, ctgs) in final_winners if bin_ == b)
            c_final, d_final, _ = compute_busco_C_D(
                bin_markers[b], owned, bin_n_total[b]
            )
            remaining = len(owned)
        else:
            still_unassigned = bin_contigs_initial[b] - set(contig_owner.keys())
            c_final, d_final, _ = compute_busco_C_D(
                bin_markers[b], still_unassigned, bin_n_total[b]
            )
            remaining = len(still_unassigned)

        # Dominant winner
        dom = bin_rels.loc[bin_rels["n_contigs_taken"].idxmax()]
        summary_rows.append({
            "loser_bin": b,
            "loser_tool": detect_tool(b),
            "status": status,
            "loser_was_candidate": was_candidate,
            "loser_initial_score": round(compute_score(c_init, d_init, alpha_d), 4),
            "n_winners_that_took_from": n_winners_took,
            "total_contigs_taken": total_taken,
            "n_contigs_initial": n_init,
            "pct_lost": round(100 * total_taken / n_init, 1),
            "contigs_remaining": remaining,
            "C_initial": round(c_init, 1),
            "D_initial": round(d_init, 1),
            "C_final": round(c_final, 1),
            "D_final": round(d_final, 1),
            "delta_C": round(c_final - c_init, 1),
            "delta_D": round(d_final - d_init, 1),
            "dominant_winner": dom["winner_bin"],
            "dominant_winner_tool": detect_tool(dom["winner_bin"]),
            "dominant_winner_pct": dom["pct_of_loser_contigs"],
        })

    summary_df = pd.DataFrame(summary_rows).sort_values(
        ["status", "pct_lost"], ascending=[True, False]
    )
    summary_df.to_csv(out_dir / "bin_summary_loss.tsv", sep="\t", index=False)
    logger.info("  Wrote: bin_summary_loss.tsv (%d rows)", len(summary_df))

    # ---- dereplicated_bins/ FASTA ----
    if write_fasta:
        fasta_dir = out_dir / "dereplicated_bins"
        fasta_dir.mkdir(parents=True, exist_ok=True)
        for b, ctgs in final_winners:
            src = bins_dir / b
            dst = fasta_dir / b
            with open(src) as fin, open(dst, "w") as fout:
                write_this = False
                for line in fin:
                    if line.startswith(">"):
                        name = line[1:].split()[0]
                        write_this = name in ctgs
                    if write_this:
                        fout.write(line)
        logger.info("  Wrote: %d cleaned FASTA bins to %s/", len(final_winners), fasta_dir)


# ============================================================================
# Summary printing
# ============================================================================

def print_summary(out_dir: Path, logger: logging.Logger) -> None:
    """Print a human-readable summary at the end of the run."""
    winners_df = pd.read_csv(out_dir / "winners.tsv", sep="\t")
    summary_df = pd.read_csv(out_dir / "bin_summary_loss.tsv", sep="\t")
    rel_df = pd.read_csv(out_dir / "bin_relationships.tsv", sep="\t")

    print()
    print("=" * 70)
    print("  REFINEMENT SUMMARY")
    print("=" * 70)
    print(f"\nFinal non-redundant bins: {len(winners_df)}")
    print("\nWinners by binning tool:")
    print(winners_df["tool"].value_counts().to_string())

    print(f"\nTotal loser-winner relationships tracked: {len(rel_df)}")

    print("\n--- TOP 10 LOSERS by % contigs lost ---")
    cols = ["loser_bin", "loser_tool", "status", "C_initial", "C_final",
            "delta_C", "n_winners_that_took_from", "pct_lost",
            "dominant_winner", "dominant_winner_pct"]
    print(summary_df.head(10)[cols].to_string(index=False))

    hq_lost = summary_df[
        (summary_df["C_initial"] >= 50)
        & (summary_df["status"].str.startswith("excluded"))
    ].sort_values("C_initial", ascending=False)
    if len(hq_lost):
        print("\n--- HIGH-QUALITY LOSERS (C_init >= 50%, excluded) ---")
        print(hq_lost.head(10)[cols].to_string(index=False))


# ============================================================================
# CLI
# ============================================================================

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bins-dir", required=True, type=Path,
        help="Directory containing input bin FASTA files (*.fa). "
             "All bins from all binners should be in the same folder, "
             "with distinguishing filename prefixes (e.g. concoct_, metabat_, remag_).",
    )
    parser.add_argument(
        "--busco-dirs", required=True, nargs="+", type=Path,
        help="One or more directories containing per-bin BUSCO output subfolders. "
             "Each subfolder must contain run_eukaryota_odb10/full_table.tsv. "
             "The subfolder name must match the bin filename (without .fa).",
    )
    parser.add_argument(
        "--out-dir", required=True, type=Path,
        help="Output directory. Created if it does not exist.",
    )
    parser.add_argument(
        "--alpha-d", type=float, default=5.0,
        help="Weight of the duplication penalty. score = (C - alpha_d * D) / 100. "
             "Higher values reject chimeric bins more aggressively. Default: 5.0.",
    )
    parser.add_argument(
        "--min-score", type=float, default=0.05,
        help="Minimum score for a bin to be considered as a candidate winner. "
             "Default: 0.05 (i.e. a bin must have at least 5%% completeness equivalent).",
    )
    parser.add_argument(
        "--min-contigs", type=int, default=10,
        help="Minimum number of contigs for a bin to remain a candidate. Default: 10.",
    )
    parser.add_argument(
        "--write-fasta", action="store_true",
        help="If set, write cleaned FASTA files of the winning bins to "
             "<out-dir>/dereplicated_bins/.",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable debug-level logging.",
    )
    return parser.parse_args()


def setup_logger(verbose: bool) -> logging.Logger:
    """Configure a basic stderr logger."""
    logger = logging.getLogger("euk_bin_refiner")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.addHandler(handler)
    return logger


def main() -> int:
    args = parse_args()
    logger = setup_logger(args.verbose)

    # Load input
    bin_contigs, bin_markers, bin_n_total = load_bins_with_busco(
        args.bins_dir, args.busco_dirs, args.min_contigs, logger
    )

    if not bin_contigs:
        logger.error("No bins to process. Exiting.")
        return 1

    # Compute initial BUSCO once (for later reporting)
    initial_busco: Dict[str, Tuple[float, float]] = {}
    for b, ctgs in bin_contigs.items():
        C, D, _ = compute_busco_C_D(bin_markers[b], ctgs, bin_n_total[b])
        initial_busco[b] = (C, D)
    initial_n_contigs = {b: len(ctgs) for b, ctgs in bin_contigs.items()}

    # Run refinement
    final_winners, contig_owner, log_rows, relationships, winner_iteration = (
        greedy_refinement(
            bin_contigs, bin_markers, bin_n_total,
            args.alpha_d, args.min_score, args.min_contigs, logger,
        )
    )

    # Write outputs
    write_outputs(
        args.out_dir, bin_contigs, bin_markers, bin_n_total,
        initial_busco, initial_n_contigs, args.alpha_d,
        final_winners, contig_owner, log_rows, relationships, winner_iteration,
        args.bins_dir, args.write_fasta, logger,
    )

    print_summary(args.out_dir, logger)
    logger.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
