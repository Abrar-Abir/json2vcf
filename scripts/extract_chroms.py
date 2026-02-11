#!/usr/bin/env python3
"""Extract specific chromosomes from a large VCF into separate files.

Single-pass streaming: reads a gzipped VCF once, writes matching lines
to per-chromosome output files.  Much faster than running zcat|awk N times.

Usage:
    python3 scripts/extract_chroms.py \
        --input /path/to/ref.vcf.gz \
        --chromosomes chr1,chr21,chrY \
        --output-dir /tmp/ref_extracts/
"""

import argparse
import gzip
import os
import sys
import time


def main():
    parser = argparse.ArgumentParser(
        description="Extract chromosomes from a (gzipped) VCF into separate files.",
    )
    parser.add_argument(
        "--input", required=True,
        help="Input VCF file (.vcf or .vcf.gz)",
    )
    parser.add_argument(
        "--chromosomes", required=True,
        help="Comma-separated list of chromosomes to extract (e.g. chr1,chr21,chrY)",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory for per-chromosome VCF files",
    )
    args = parser.parse_args()

    chroms = set(args.chromosomes.split(","))
    os.makedirs(args.output_dir, exist_ok=True)

    # Open per-chromosome output files
    out_files = {}
    for chrom in chroms:
        path = os.path.join(args.output_dir, f"{chrom}.vcf")
        out_files[chrom] = open(path, "w")

    # Stream input (gzipped or plain)
    if args.input.endswith(".gz"):
        fh = gzip.open(args.input, "rt", encoding="utf-8")
    else:
        fh = open(args.input, "r", encoding="utf-8")

    total_lines = 0
    written = {c: 0 for c in chroms}
    t0 = time.time()

    try:
        for line in fh:
            if line.startswith("#"):
                continue
            total_lines += 1

            # Fast CHROM extraction — only find the first tab
            tab_idx = line.index("\t")
            chrom = line[:tab_idx]

            if chrom in out_files:
                out_files[chrom].write(line)
                written[chrom] += 1

            if total_lines % 5_000_000 == 0:
                elapsed = time.time() - t0
                print(
                    f"  {total_lines:,} lines processed in {elapsed:.0f}s",
                    file=sys.stderr,
                )
    finally:
        fh.close()
        for f in out_files.values():
            f.close()

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s. {total_lines:,} total lines processed.",
          file=sys.stderr)
    for chrom in sorted(written.keys()):
        print(f"  {chrom}: {written[chrom]:,} lines written -> "
              f"{os.path.join(args.output_dir, chrom + '.vcf')}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
