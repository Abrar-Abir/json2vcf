#!/usr/bin/env python3
"""Compare variant positions between two VCF files.

Extracts (CHROM, POS, REF, ALT) tuples from each VCF, expanding
multi-allelic rows into per-allele tuples, and reports set overlap.

Usage:
    python3 scripts/compare_positions.py \
        --json2vcf /tmp/json2vcf_chrY.vcf \
        --reference /tmp/ref_extracts/chrY.vcf \
        [--max-examples 20]
"""

import argparse
import gzip
import sys
import time
from collections import defaultdict


# ---------------------------------------------------------------------------
# VCF reading helpers
# ---------------------------------------------------------------------------

def open_vcf(path):
    """Return a line iterator for a VCF file (plain or gzipped)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def extract_variant_tuples(path):
    """Stream a VCF and return a set of (CHROM, POS, REF, ALT) tuples.

    Multi-allelic ALTs (comma-separated) are expanded into individual tuples.
    ALT values of '.' and '*' are skipped.

    Returns (tuples_set, row_count).
    """
    tuples_set = set()
    row_count = 0

    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            row_count += 1

            # Only need first 5 columns: CHROM POS ID REF ALT
            fields = line.rstrip("\n").split("\t", 5)
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt_field = fields[4]

            if alt_field == ".":
                continue

            for alt in alt_field.split(","):
                alt = alt.strip()
                if alt in (".", "*"):
                    continue
                tuples_set.add((chrom, pos, ref, alt))

    return tuples_set, row_count


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------

def compute_pos_range(tuples_set):
    """Return (min_pos, max_pos) from a set of variant tuples."""
    if not tuples_set:
        return (0, 0)
    positions = [t[1] for t in tuples_set]
    return min(positions), max(positions)


def filter_by_range(tuples_set, min_pos, max_pos):
    """Keep only tuples where min_pos <= POS <= max_pos."""
    return {t for t in tuples_set if min_pos <= t[1] <= max_pos}


def analyze_discrepancies(only_j2v, only_ref):
    """Categorize discrepancies into position-matched vs truly-unique.

    Position-matched: same (CHROM, POS) appears in both exclusive sets,
    meaning the variant exists at that position in both files but with
    different REF/ALT — likely allele normalization differences.
    """
    j2v_positions = defaultdict(list)
    for t in only_j2v:
        j2v_positions[(t[0], t[1])].append(t)

    ref_positions = defaultdict(list)
    for t in only_ref:
        ref_positions[(t[0], t[1])].append(t)

    shared_positions = set(j2v_positions.keys()) & set(ref_positions.keys())

    pos_matched_j2v = []
    pos_matched_ref = []
    for pos_key in sorted(shared_positions):
        pos_matched_j2v.extend(j2v_positions[pos_key])
        pos_matched_ref.extend(ref_positions[pos_key])

    truly_unique_j2v = [t for t in only_j2v
                        if (t[0], t[1]) not in shared_positions]
    truly_unique_ref = [t for t in only_ref
                        if (t[0], t[1]) not in shared_positions]

    return {
        "shared_positions": len(shared_positions),
        "pos_matched_j2v": sorted(pos_matched_j2v),
        "pos_matched_ref": sorted(pos_matched_ref),
        "truly_unique_j2v": sorted(truly_unique_j2v),
        "truly_unique_ref": sorted(truly_unique_ref),
    }


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

def fmt_tuple(t):
    """Format a variant tuple as a readable string."""
    return f"  {t[0]}\t{t[1]}\t{t[2]}\t{t[3]}"


def pct(num, denom):
    """Format a percentage string, handling zero denominator."""
    if denom == 0:
        return "N/A"
    return f"{100.0 * num / denom:.2f}%"


def print_report(j2v_path, ref_path, j2v_tuples, ref_tuples,
                 j2v_rows, ref_rows, j2v_range, ref_range,
                 ref_range_filtered, shared, only_j2v, only_ref,
                 analysis, max_examples):
    """Print the full concordance report to stdout."""
    union_size = len(shared) + len(only_j2v) + len(only_ref)

    print("=" * 60)
    print("POSITION-LEVEL CONCORDANCE REPORT")
    print("=" * 60)
    print(f"json2vcf file:  {j2v_path}")
    print(f"Reference file: {ref_path}")
    print()

    # --- Variant counts ---
    print("VARIANT COUNTS")
    print(f"  json2vcf VCF rows:        {j2v_rows:>12,}")
    print(f"  Reference VCF rows:       {ref_rows:>12,}")
    print(f"  json2vcf variant tuples:  {len(j2v_tuples):>12,}  "
          f"(after multi-allelic expansion)")
    print(f"  Reference variant tuples: {len(ref_tuples):>12,}  "
          f"(after multi-allelic expansion)")
    if ref_range_filtered is not None:
        print(f"  Reference tuples (range-filtered): "
              f"{ref_range_filtered:>8,}")
    print()

    # --- Position ranges ---
    print("POSITION RANGES")
    print(f"  json2vcf:  {j2v_range[0]:,} - {j2v_range[1]:,}")
    print(f"  Reference: {ref_range[0]:,} - {ref_range[1]:,}")
    if ref_range[1] > j2v_range[1] or ref_range[0] < j2v_range[0]:
        print("  NOTE: Reference covers a wider range; auto-filtered "
              "to json2vcf range for comparison.")
    print()

    # --- Set comparison ---
    print("SET COMPARISON")
    print(f"  Shared (in both):    {len(shared):>12,}  ({pct(len(shared), union_size)})")
    print(f"  Only in json2vcf:    {len(only_j2v):>12,}  ({pct(len(only_j2v), len(j2v_tuples))})")
    print(f"  Only in reference:   {len(only_ref):>12,}  ({pct(len(only_ref), ref_range_filtered or len(ref_tuples))})")
    if union_size > 0:
        jaccard = len(shared) / union_size
        print(f"  Jaccard similarity:  {jaccard:>12.6f}")
    print()

    # --- Discrepancy analysis ---
    print("DISCREPANCY ANALYSIS")
    print(f"  Position-matched but allele-different: "
          f"{analysis['shared_positions']:,}")
    print(f"  Truly unique to json2vcf:  {len(analysis['truly_unique_j2v']):,}")
    print(f"  Truly unique to reference: {len(analysis['truly_unique_ref']):,}")
    print()

    # --- Examples ---
    if analysis["pos_matched_j2v"]:
        n = min(max_examples, len(analysis["pos_matched_j2v"]))
        print(f"EXAMPLES: Position-matched, allele-different (first {n})")
        # Group by (CHROM, POS)
        seen = set()
        count = 0
        for t in analysis["pos_matched_j2v"]:
            key = (t[0], t[1])
            if key in seen:
                continue
            seen.add(key)
            j2v_at_pos = [x for x in analysis["pos_matched_j2v"]
                          if (x[0], x[1]) == key]
            ref_at_pos = [x for x in analysis["pos_matched_ref"]
                          if (x[0], x[1]) == key]
            print(f"  POS {key[0]}:{key[1]}")
            for x in j2v_at_pos:
                print(f"    json2vcf:  REF={x[2]}  ALT={x[3]}")
            for x in ref_at_pos:
                print(f"    reference: REF={x[2]}  ALT={x[3]}")
            count += 1
            if count >= max_examples:
                break
        print()

    if analysis["truly_unique_j2v"]:
        n = min(max_examples, len(analysis["truly_unique_j2v"]))
        print(f"EXAMPLES: Truly unique to json2vcf (first {n})")
        for t in analysis["truly_unique_j2v"][:n]:
            print(fmt_tuple(t))
        print()

    if analysis["truly_unique_ref"]:
        n = min(max_examples, len(analysis["truly_unique_ref"]))
        print(f"EXAMPLES: Truly unique to reference (first {n})")
        for t in analysis["truly_unique_ref"][:n]:
            print(fmt_tuple(t))
        print()

    # --- Pass/fail summary ---
    overlap_pct = 100.0 * len(shared) / len(j2v_tuples) if j2v_tuples else 0
    print("-" * 60)
    if overlap_pct >= 99.0:
        print(f"RESULT: PASS  (overlap = {overlap_pct:.2f}% of json2vcf tuples)")
    else:
        print(f"RESULT: REVIEW NEEDED  (overlap = {overlap_pct:.2f}% of json2vcf tuples, target >= 99%)")
    print("-" * 60)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare variant positions between two VCF files.",
    )
    parser.add_argument(
        "--json2vcf", required=True,
        help="Path to json2vcf output VCF (plain or .gz)",
    )
    parser.add_argument(
        "--reference", required=True,
        help="Path to reference VCF extract (plain or .gz)",
    )
    parser.add_argument(
        "--max-examples", type=int, default=20,
        help="Max examples to print per discrepancy category (default: 20)",
    )
    args = parser.parse_args()

    t0 = time.time()

    # --- Load json2vcf tuples ---
    print("Loading json2vcf variants...", file=sys.stderr)
    j2v_tuples, j2v_rows = extract_variant_tuples(args.json2vcf)
    j2v_range = compute_pos_range(j2v_tuples)
    print(f"  {len(j2v_tuples):,} tuples from {j2v_rows:,} rows "
          f"(range: {j2v_range[0]:,}-{j2v_range[1]:,})",
          file=sys.stderr)

    # --- Load reference tuples ---
    print("Loading reference variants...", file=sys.stderr)
    ref_tuples_all, ref_rows = extract_variant_tuples(args.reference)
    ref_range = compute_pos_range(ref_tuples_all)
    print(f"  {len(ref_tuples_all):,} tuples from {ref_rows:,} rows "
          f"(range: {ref_range[0]:,}-{ref_range[1]:,})",
          file=sys.stderr)

    # --- Auto-range filter reference to json2vcf range ---
    ref_tuples = filter_by_range(ref_tuples_all, j2v_range[0], j2v_range[1])
    ref_filtered_count = len(ref_tuples)
    if ref_filtered_count < len(ref_tuples_all):
        print(f"  Filtered reference to json2vcf range: "
              f"{ref_filtered_count:,} tuples kept",
              file=sys.stderr)
    else:
        ref_filtered_count = None  # no filtering needed

    # --- Compare ---
    print("Comparing...", file=sys.stderr)
    shared = j2v_tuples & ref_tuples
    only_j2v = j2v_tuples - ref_tuples
    only_ref = ref_tuples - j2v_tuples

    # --- Analyze discrepancies ---
    analysis = analyze_discrepancies(only_j2v, only_ref)

    elapsed = time.time() - t0
    print(f"Done in {elapsed:.1f}s\n", file=sys.stderr)

    # --- Print report ---
    print_report(
        j2v_path=args.json2vcf,
        ref_path=args.reference,
        j2v_tuples=j2v_tuples,
        ref_tuples=ref_tuples_all,
        j2v_rows=j2v_rows,
        ref_rows=ref_rows,
        j2v_range=j2v_range,
        ref_range=ref_range,
        ref_range_filtered=ref_filtered_count,
        shared=shared,
        only_j2v=only_j2v,
        only_ref=only_ref,
        analysis=analysis,
        max_examples=args.max_examples,
    )


if __name__ == "__main__":
    main()
