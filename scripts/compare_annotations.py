#!/usr/bin/env python3
"""Compare annotation values between json2vcf output and a reference VCF.

Two-pass approach:
  Pass 1: Load json2vcf VCF into a dict keyed by (CHROM, POS, REF, ALT)
  Pass 2: Stream reference VCF, look up each variant, compare annotation fields

Usage:
    python3 scripts/compare_annotations.py \
        --json2vcf /tmp/json2vcf_chrY.vcf \
        --reference /tmp/ref_extracts/chrY.vcf \
        [--output report.txt] \
        [--tsv discordant.tsv] \
        [--max-examples 20] \
        [--verbose]
"""

import argparse
import gzip
import math
import sys
import time
from collections import Counter


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# CSQ subfield indices (json2vcf format, 19 pipe-separated fields)
CSQ_ALLELE_IDX = 0
CSQ_SYMBOL_IDX = 2
CSQ_CANONICAL_IDX = 16
CSQ_POLYPHEN_IDX = 17
CSQ_SIFT_IDX = 18

# ANN subfield indices (snpEff format)
ANN_ALLELE_IDX = 0
ANN_GENE_IDX = 3

# Map json2vcf SIFT/PolyPhen text predictions to single-letter codes
SIFT_MAP = {
    "deleterious": "D",
    "tolerated": "T",
    "deleterious_low_confidence": "D",
    "tolerated_low_confidence": "T",
}

POLYPHEN_MAP = {
    "probably_damaging": "D",
    "possibly_damaging": "P",
    "benign": "B",
}


class FieldSpec:
    """Describes one annotation field comparison."""

    __slots__ = ("name", "j2v_key", "ref_key", "method", "abs_tol", "rel_tol")

    def __init__(self, name, j2v_key, ref_key, method,
                 abs_tol=0.01, rel_tol=0.05):
        self.name = name
        self.j2v_key = j2v_key
        self.ref_key = ref_key
        self.method = method
        self.abs_tol = abs_tol
        self.rel_tol = rel_tol


FIELD_SPECS = [
    FieldSpec("gnomAD AF (all)", "gnomAD_AF", "gnomadWGS_AF", "numeric"),
    FieldSpec("gnomAD AF (AFR)", "gnomAD_AFR_AF", "gnomadWGS_AF_AFR", "numeric"),
    FieldSpec("gnomAD AF (AMR)", "gnomAD_AMR_AF", "gnomadWGS_AF_AMR", "numeric"),
    FieldSpec("gnomAD AF (EAS)", "gnomAD_EAS_AF", "gnomadWGS_AF_EAS", "numeric"),
    FieldSpec("gnomAD AF (EUR/NFE)", "gnomAD_EUR_AF", "gnomadWGS_AF_NFE", "numeric"),
    FieldSpec("gnomAD AF (SAS)", "gnomAD_SAS_AF", "gnomadWGS_AF_SAS", "numeric"),
    FieldSpec("ClinVar significance", "CLINVAR_SIG", "CLNSIG", "categorical"),
    FieldSpec("ClinVar review status", "CLINVAR_REVSTAT", "CLNREVSTAT", "categorical"),
    FieldSpec("dbSNP rsID", "dbSNP_IDs", "rs_ids", "set"),
    FieldSpec("REVEL score", "REVEL", "dbNSFP_REVEL_score", "numeric",
              abs_tol=0.01, rel_tol=0.0),
    FieldSpec("TOPMed AF", "TOPMed_AF", "TOPMED", "numeric"),
    FieldSpec("Gene symbol", "GENE_SYMBOLS", "ANN_GENES", "set"),
    FieldSpec("SIFT", "SIFT_pred", "dbNSFP_SIFT_pred", "categorical"),
    FieldSpec("PolyPhen", "PolyPhen_pred", "dbNSFP_Polyphen2_HDIV_pred", "categorical"),
]


class FieldStats:
    """Accumulates per-field comparison statistics."""

    def __init__(self, name, max_examples=20):
        self.name = name
        self.max_examples = max_examples
        self.both_present = 0
        self.concordant = 0
        self.discordant = 0
        self.only_in_j2v = 0
        self.only_in_ref = 0
        # Numeric running sums (for Pearson correlation)
        self.sum_x = 0.0
        self.sum_y = 0.0
        self.sum_x2 = 0.0
        self.sum_y2 = 0.0
        self.sum_xy = 0.0
        self.sum_abs_diff = 0.0
        self.max_abs_diff = 0.0
        # Categorical mismatch categories
        self.mismatch_counts = Counter()
        # Discordant examples: list of (variant_key, j2v_val, ref_val)
        self.examples = []

    def add_example(self, key, j2v_val, ref_val):
        if len(self.examples) < self.max_examples:
            self.examples.append((key, j2v_val, ref_val))

    def correlation(self):
        n = self.both_present
        if n < 2:
            return None
        num = n * self.sum_xy - self.sum_x * self.sum_y
        dx = n * self.sum_x2 - self.sum_x ** 2
        dy = n * self.sum_y2 - self.sum_y ** 2
        if dx <= 0 or dy <= 0:
            return None
        return num / math.sqrt(dx * dy)

    def mean_abs_diff(self):
        if self.both_present == 0:
            return None
        return self.sum_abs_diff / self.both_present


# ---------------------------------------------------------------------------
# VCF I/O helpers
# ---------------------------------------------------------------------------

def open_vcf(path):
    """Return a line iterator for a VCF file (plain or gzipped)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_info(info_str):
    """Parse VCF INFO column into {key: value_string} dict."""
    if info_str == ".":
        return {}
    result = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            result[k] = v
        else:
            result[part] = True
    return result


def _unescape_info_value(value):
    """Reverse percent-encoding from json2vcf mapper (decode in reverse order)."""
    return (
        value.replace("%2C", ",")
        .replace("%3D", "=")
        .replace("%3B", ";")
        .replace("%20", " ")
        .replace("%25", "%")
    )


def _get_per_allele_value(info_dict, key, allele_idx):
    """Extract a string value for a specific allele from a Number=A field."""
    raw = info_dict.get(key)
    if raw is None:
        return None
    parts = raw.split(",")
    if allele_idx >= len(parts):
        return None
    val = parts[allele_idx]
    return None if val == "." else val


def _per_allele_float(info_dict, key, allele_idx):
    """Extract a float for a specific allele from a Number=A field."""
    val = _get_per_allele_value(info_dict, key, allele_idx)
    if val is None:
        return None
    try:
        return float(val)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# json2vcf VCF loader
# ---------------------------------------------------------------------------

def _extract_prediction(pred_str):
    """Extract prediction name from 'prediction(score)' format."""
    if not pred_str:
        return None
    paren = pred_str.find("(")
    if paren >= 0:
        return pred_str[:paren].strip().lower()
    return pred_str.strip().lower()


def _parse_csq_for_allele(csq_str, alt_allele):
    """Parse CSQ string and extract annotations for a specific alt allele.

    Returns dict with keys: GENE_SYMBOLS (set), SIFT_pred (str), PolyPhen_pred (str).
    """
    result = {"GENE_SYMBOLS": set(), "SIFT_pred": None, "PolyPhen_pred": None}
    if not csq_str:
        return result

    transcripts = csq_str.split(",")
    canonical_sift = None
    canonical_pp = None
    any_sift = None
    any_pp = None

    for entry in transcripts:
        fields = entry.split("|")
        if len(fields) < 19:
            continue

        # Match allele: for SNVs the CSQ allele matches alt_allele directly;
        # for deletions it may be "-"
        csq_allele = fields[CSQ_ALLELE_IDX]
        if csq_allele != alt_allele and csq_allele != "-":
            continue

        symbol = fields[CSQ_SYMBOL_IDX]
        if symbol:
            result["GENE_SYMBOLS"].add(symbol)

        is_canonical = fields[CSQ_CANONICAL_IDX] == "YES"

        pp_raw = fields[CSQ_POLYPHEN_IDX]
        sift_raw = fields[CSQ_SIFT_IDX]

        if pp_raw:
            pred = _extract_prediction(pp_raw)
            mapped = POLYPHEN_MAP.get(pred)
            if mapped is not None:
                if is_canonical:
                    canonical_pp = mapped
                elif any_pp is None:
                    any_pp = mapped

        if sift_raw:
            pred = _extract_prediction(sift_raw)
            mapped = SIFT_MAP.get(pred)
            if mapped is not None:
                if is_canonical:
                    canonical_sift = mapped
                elif any_sift is None:
                    any_sift = mapped

    result["SIFT_pred"] = canonical_sift if canonical_sift else any_sift
    result["PolyPhen_pred"] = canonical_pp if canonical_pp else any_pp
    return result


def _normalize_clinvar(raw, separator="&"):
    """Normalize a ClinVar value: unescape, split, lowercase, deduplicate, sort, rejoin."""
    if raw is None:
        return None
    decoded = _unescape_info_value(raw)
    parts = decoded.replace("/", separator).split(separator)
    parts = sorted(set(p.strip().lower().replace(" ", "_") for p in parts if p.strip()))
    return separator.join(parts) if parts else None


def load_json2vcf(path, verbose=False):
    """Load json2vcf VCF into {(chrom, pos, ref, alt): annotation_dict}.

    Multi-allelic rows are expanded to per-allele entries.
    """
    variants = {}
    row_count = 0

    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            row_count += 1

            cols = line.rstrip("\n").split("\t", 8)
            chrom = cols[0]
            pos = int(cols[1])
            id_col = cols[2]
            ref = cols[3]
            alt_field = cols[4]
            info_str = cols[7]

            if alt_field == ".":
                continue

            alt_list = alt_field.split(",")
            info = parse_info(info_str)

            # dbSNP IDs from ID column (shared across all alleles at this pos)
            dbsnp_ids = set()
            if id_col != ".":
                dbsnp_ids = {x for x in id_col.split(";") if x.startswith("rs")}

            # ClinVar (position-level, not per-allele)
            clinvar_sig = _normalize_clinvar(info.get("CLINVAR_SIG"))
            clinvar_revstat = _normalize_clinvar(info.get("CLINVAR_REVSTAT"))

            # CSQ string (will be parsed per-allele)
            csq_str = info.get("CSQ")

            for i, alt in enumerate(alt_list):
                alt = alt.strip()
                if alt in (".", "*"):
                    continue

                key = (chrom, pos, ref, alt)

                # Per-allele numeric fields
                annot = {
                    "gnomAD_AF": _per_allele_float(info, "gnomAD_AF", i),
                    "gnomAD_AFR_AF": _per_allele_float(info, "gnomAD_AFR_AF", i),
                    "gnomAD_AMR_AF": _per_allele_float(info, "gnomAD_AMR_AF", i),
                    "gnomAD_EAS_AF": _per_allele_float(info, "gnomAD_EAS_AF", i),
                    "gnomAD_EUR_AF": _per_allele_float(info, "gnomAD_EUR_AF", i),
                    "gnomAD_SAS_AF": _per_allele_float(info, "gnomAD_SAS_AF", i),
                    "REVEL": _per_allele_float(info, "REVEL", i),
                    "TOPMed_AF": _per_allele_float(info, "TOPMed_AF", i),
                    "CLINVAR_SIG": clinvar_sig,
                    "CLINVAR_REVSTAT": clinvar_revstat,
                    "dbSNP_IDs": dbsnp_ids,
                }

                # CSQ-derived fields
                csq_annot = _parse_csq_for_allele(csq_str, alt) if csq_str else {
                    "GENE_SYMBOLS": set(), "SIFT_pred": None, "PolyPhen_pred": None
                }
                annot["GENE_SYMBOLS"] = csq_annot["GENE_SYMBOLS"]
                annot["SIFT_pred"] = csq_annot["SIFT_pred"]
                annot["PolyPhen_pred"] = csq_annot["PolyPhen_pred"]

                variants[key] = annot

            if verbose and row_count % 100_000 == 0:
                print(f"  json2vcf: {row_count:,} rows loaded...",
                      file=sys.stderr)

    if verbose:
        print(f"  json2vcf: {row_count:,} rows -> {len(variants):,} variant entries",
              file=sys.stderr)
    return variants


# ---------------------------------------------------------------------------
# Reference VCF parser
# ---------------------------------------------------------------------------

def _safe_float(val):
    """Parse a float, returning None on failure or missing."""
    if val is None or val == "." or val == "":
        return None
    try:
        return float(val)
    except ValueError:
        return None


def _parse_ref_float(info_dict, key, allele_idx=0):
    """Extract float from reference INFO. Handle possible per-allele values."""
    raw = info_dict.get(key)
    if raw is None:
        return None
    parts = raw.split(",")
    if allele_idx < len(parts):
        return _safe_float(parts[allele_idx])
    return _safe_float(parts[0]) if parts else None


def _parse_ref_dbsnp(info_dict):
    """Extract rsID set from reference INFO rs_ids field."""
    raw = info_dict.get("rs_ids")
    if raw is None or raw == ".":
        return set()
    return {x.strip() for x in raw.split(",") if x.strip().startswith("rs")}


def _parse_ann_genes(ann_str, alt_allele):
    """Extract gene symbols from snpEff ANN field for a specific allele."""
    if not ann_str:
        return set()
    genes = set()
    for entry in ann_str.split(","):
        fields = entry.split("|")
        if len(fields) <= ANN_GENE_IDX:
            continue
        # snpEff allele field may differ from VCF ALT in complex variants;
        # collect all genes as fallback
        gene = fields[ANN_GENE_IDX]
        if gene:
            genes.add(gene)
    return genes


def _normalize_ref_clinvar(raw, separator="&"):
    """Normalize reference ClinVar value: lowercase, deduplicate, split on / and _, sort."""
    if raw is None or raw == "." or raw == "":
        return None
    parts = raw.replace("/", separator).split(separator)
    parts = sorted(set(p.strip().lower().replace(" ", "_") for p in parts if p.strip()))
    return separator.join(parts) if parts else None


def _parse_ref_prediction(raw):
    """Parse dbNSFP prediction field that may have multiple values.

    dbNSFP uses ',' or ';' or '.' to separate predictions from multiple
    transcripts, e.g. "T,", ",D,D,D,D", "B,B,".
    Returns the most common non-empty prediction (uppercase single letter).
    """
    if raw is None or raw == "." or raw == "":
        return None
    # Split on comma first (most common dbNSFP separator)
    parts = raw.replace(";", ",").split(",")
    # Collect valid single-letter predictions
    preds = []
    for p in parts:
        p = p.strip()
        if p and p != ".":
            preds.append(p.upper())
    if not preds:
        return None
    # Return the most common prediction
    counts = Counter(preds)
    return counts.most_common(1)[0][0]


def parse_ref_line(line):
    """Parse one reference VCF data line.

    Returns list of (key, annotation_dict) tuples, one per alt allele.
    """
    cols = line.rstrip("\n").split("\t", 8)
    chrom = cols[0]
    pos = int(cols[1])
    ref = cols[3]
    alt_field = cols[4]
    info_str = cols[7]

    if alt_field == ".":
        return []

    alt_list = alt_field.split(",")
    info = parse_info(info_str)

    # Position-level fields
    clinvar_sig = _normalize_ref_clinvar(info.get("CLNSIG"))
    clinvar_revstat = _normalize_ref_clinvar(info.get("CLNREVSTAT"))
    ref_dbsnp = _parse_ref_dbsnp(info)
    ann_str = info.get("ANN")
    sift_pred = _parse_ref_prediction(info.get("dbNSFP_SIFT_pred"))
    pp_pred = _parse_ref_prediction(info.get("dbNSFP_Polyphen2_HDIV_pred"))

    results = []
    for i, alt in enumerate(alt_list):
        alt = alt.strip()
        if alt in (".", "*"):
            continue

        key = (chrom, pos, ref, alt)
        annot = {
            "gnomadWGS_AF": _parse_ref_float(info, "gnomadWGS_AF", i),
            "gnomadWGS_AF_AFR": _parse_ref_float(info, "gnomadWGS_AF_AFR", i),
            "gnomadWGS_AF_AMR": _parse_ref_float(info, "gnomadWGS_AF_AMR", i),
            "gnomadWGS_AF_EAS": _parse_ref_float(info, "gnomadWGS_AF_EAS", i),
            "gnomadWGS_AF_NFE": _parse_ref_float(info, "gnomadWGS_AF_NFE", i),
            "gnomadWGS_AF_SAS": _parse_ref_float(info, "gnomadWGS_AF_SAS", i),
            "CLNSIG": clinvar_sig,
            "CLNREVSTAT": clinvar_revstat,
            "rs_ids": ref_dbsnp,
            "dbNSFP_REVEL_score": _parse_ref_float(info, "dbNSFP_REVEL_score", i),
            "TOPMED": _parse_ref_float(info, "TOPMED", i),
            "ANN_GENES": _parse_ann_genes(ann_str, alt),
            "dbNSFP_SIFT_pred": sift_pred,
            "dbNSFP_Polyphen2_HDIV_pred": pp_pred,
        }
        results.append((key, annot))

    return results


# ---------------------------------------------------------------------------
# Comparison engine
# ---------------------------------------------------------------------------

def _compare_numeric(j2v_val, ref_val, abs_tol, rel_tol):
    """Compare two numeric values. Returns result string or None if both absent."""
    if j2v_val is None and ref_val is None:
        return None
    if j2v_val is None:
        return "only_ref"
    if ref_val is None:
        return "only_j2v"
    diff = abs(j2v_val - ref_val)
    if diff < abs_tol:
        return "concordant"
    denom = max(abs(j2v_val), abs(ref_val))
    if rel_tol > 0 and denom > 0 and diff / denom < rel_tol:
        return "concordant"
    return "discordant"


def _compare_categorical(j2v_val, ref_val):
    """Compare two categorical values. Returns result string or None."""
    if j2v_val is None and ref_val is None:
        return None
    if j2v_val is None:
        return "only_ref"
    if ref_val is None:
        return "only_j2v"
    if j2v_val == ref_val:
        return "concordant"
    return "discordant"


def _compare_sets(j2v_set, ref_set):
    """Compare two sets. Concordant if any overlap. Returns result or None."""
    j2v_empty = not j2v_set
    ref_empty = not ref_set
    if j2v_empty and ref_empty:
        return None
    if j2v_empty:
        return "only_ref"
    if ref_empty:
        return "only_j2v"
    if j2v_set & ref_set:
        return "concordant"
    return "discordant"


def compare_field(spec, j2v_annot, ref_annot, stats):
    """Compare one field between two annotation dicts, updating stats."""
    j2v_val = j2v_annot.get(spec.j2v_key)
    ref_val = ref_annot.get(spec.ref_key)

    if spec.method == "numeric":
        result = _compare_numeric(j2v_val, ref_val, spec.abs_tol, spec.rel_tol)
    elif spec.method == "categorical":
        result = _compare_categorical(j2v_val, ref_val)
    elif spec.method == "set":
        result = _compare_sets(j2v_val, ref_val)
    else:
        return

    if result is None:
        return

    if result == "concordant":
        stats.both_present += 1
        stats.concordant += 1
    elif result == "discordant":
        stats.both_present += 1
        stats.discordant += 1
        stats.mismatch_counts[(str(j2v_val), str(ref_val))] += 1
    elif result == "only_j2v":
        stats.only_in_j2v += 1
    elif result == "only_ref":
        stats.only_in_ref += 1

    # Numeric running sums for correlation
    if spec.method == "numeric" and result in ("concordant", "discordant"):
        x, y = float(j2v_val), float(ref_val)
        stats.sum_x += x
        stats.sum_y += y
        stats.sum_x2 += x * x
        stats.sum_y2 += y * y
        stats.sum_xy += x * y
        diff = abs(x - y)
        stats.sum_abs_diff += diff
        if diff > stats.max_abs_diff:
            stats.max_abs_diff = diff

    # Collect discordant examples
    if result == "discordant":
        stats.add_example(None, j2v_val, ref_val)


def compare_with_reference(j2v_dict, ref_path, field_specs, max_examples,
                           verbose=False):
    """Stream reference VCF, compare against json2vcf dict.

    Returns (stats_dict, total_shared, total_only_j2v, total_only_ref).
    """
    stats = {s.name: FieldStats(s.name, max_examples) for s in field_specs}
    seen_keys = set()
    total_shared = 0
    total_only_ref = 0
    ref_row_count = 0

    with open_vcf(ref_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            ref_row_count += 1

            for key, ref_annot in parse_ref_line(line):
                if key in j2v_dict:
                    seen_keys.add(key)
                    total_shared += 1
                    j2v_annot = j2v_dict[key]

                    # Store variant key in examples
                    for spec in field_specs:
                        st = stats[spec.name]
                        old_len = len(st.examples)
                        compare_field(spec, j2v_annot, ref_annot, st)
                        # Backfill the variant key into the last example
                        if len(st.examples) > old_len:
                            idx = len(st.examples) - 1
                            _, jv, rv = st.examples[idx]
                            st.examples[idx] = (key, jv, rv)
                else:
                    total_only_ref += 1

            if verbose and ref_row_count % 100_000 == 0:
                print(f"  reference: {ref_row_count:,} rows processed...",
                      file=sys.stderr)

    total_only_j2v = len(j2v_dict) - len(seen_keys)

    if verbose:
        print(f"  reference: {ref_row_count:,} rows total",
              file=sys.stderr)
        print(f"  shared: {total_shared:,}  only_j2v: {total_only_j2v:,}  "
              f"only_ref: {total_only_ref:,}", file=sys.stderr)

    return stats, total_shared, total_only_j2v, total_only_ref


# ---------------------------------------------------------------------------
# Report output
# ---------------------------------------------------------------------------

def _pct(num, denom):
    if denom == 0:
        return "N/A"
    return f"{100.0 * num / denom:.1f}%"


def _fmt_key(key):
    return f"{key[0]}:{key[1]} {key[2]}>{key[3]}"


def print_report(stats, total_j2v, total_shared, total_only_j2v,
                 total_only_ref, j2v_path, ref_path, field_specs,
                 max_examples, out=None):
    """Print the annotation concordance report."""
    if out is None:
        out = sys.stdout

    p = lambda *a, **kw: print(*a, file=out, **kw)

    p("=" * 65)
    p("ANNOTATION CONCORDANCE REPORT")
    p("=" * 65)
    p(f"json2vcf file:  {j2v_path}")
    p(f"Reference file: {ref_path}")
    p()

    # Position match summary
    total_ref = total_shared + total_only_ref
    total_union = total_shared + total_only_j2v + total_only_ref
    p("POSITION MATCH")
    p(f"  Total json2vcf variants:  {total_j2v:>12,}")
    p(f"  Total reference variants: {total_ref:>12,}")
    p(f"  Shared:                   {total_shared:>12,}  "
      f"({_pct(total_shared, total_union)})")
    p(f"  Only in json2vcf:         {total_only_j2v:>12,}  "
      f"({_pct(total_only_j2v, total_j2v)})")
    p(f"  Only in reference:        {total_only_ref:>12,}  "
      f"({_pct(total_only_ref, total_ref)})")
    p()

    # Per-field blocks
    for spec in field_specs:
        st = stats[spec.name]
        total_compared = st.both_present + st.only_in_j2v + st.only_in_ref
        if total_compared == 0:
            p(f"{spec.name} ({spec.j2v_key} vs {spec.ref_key})")
            p("  No data in either file for this field.")
            p()
            continue

        p(f"{spec.name} ({spec.j2v_key} vs {spec.ref_key})")
        p(f"  Both present:     {st.both_present:>10,}")
        p(f"  Concordant:       {st.concordant:>10,}  "
          f"({_pct(st.concordant, st.both_present)})")
        p(f"  Discordant:       {st.discordant:>10,}  "
          f"({_pct(st.discordant, st.both_present)})")
        p(f"  Only in json2vcf: {st.only_in_j2v:>10,}")
        p(f"  Only in reference:{st.only_in_ref:>10,}")

        if spec.method == "numeric" and st.both_present > 0:
            mean_diff = st.mean_abs_diff()
            corr = st.correlation()
            if mean_diff is not None:
                p(f"  Mean abs diff:    {mean_diff:.6g}")
            p(f"  Max abs diff:     {st.max_abs_diff:.6g}")
            if corr is not None:
                p(f"  Correlation:      {corr:.6f}")

        if spec.method == "categorical" and st.mismatch_counts:
            p("  Top mismatches:")
            for (jv, rv), cnt in st.mismatch_counts.most_common(5):
                p(f"    \"{jv}\" vs \"{rv}\": {cnt}")

        if spec.method == "set" and st.mismatch_counts:
            p("  Top mismatches:")
            for (jv, rv), cnt in st.mismatch_counts.most_common(5):
                p(f"    {jv} vs {rv}: {cnt}")

        p()

    # Discordant examples
    has_examples = any(stats[s.name].examples for s in field_specs)
    if has_examples:
        p("=" * 65)
        p(f"DISCORDANT EXAMPLES (first {max_examples} per field)")
        p("=" * 65)
        for spec in field_specs:
            st = stats[spec.name]
            if not st.examples:
                continue
            p(f"\n{spec.name}:")
            for key, jv, rv in st.examples:
                loc = _fmt_key(key) if key else "?"
                if spec.method == "numeric" and jv is not None and rv is not None:
                    diff = abs(float(jv) - float(rv))
                    p(f"  {loc}  json2vcf={jv}  ref={rv}  (diff={diff:.6g})")
                else:
                    p(f"  {loc}  json2vcf={jv}  ref={rv}")
        p()

    # Summary pass/fail per field
    p("-" * 65)
    p("FIELD SUMMARY")
    p("-" * 65)
    for spec in field_specs:
        st = stats[spec.name]
        if st.both_present == 0:
            p(f"  {spec.name:<30s}  NO DATA")
            continue
        conc_pct = 100.0 * st.concordant / st.both_present
        p(f"  {spec.name:<30s}  {conc_pct:6.1f}% concordant  "
          f"({st.concordant:,}/{st.both_present:,})")
    p("-" * 65)


def write_tsv(stats, field_specs, path):
    """Write discordant examples to a TSV file."""
    with open(path, "w") as f:
        f.write("CHROM\tPOS\tREF\tALT\tFIELD\tJSON2VCF_VALUE\t"
                "REFERENCE_VALUE\tDIFF\n")
        for spec in field_specs:
            st = stats[spec.name]
            for key, jv, rv in st.examples:
                chrom, pos, ref, alt = key if key else (".", 0, ".", ".")
                diff = ""
                if spec.method == "numeric" and jv is not None and rv is not None:
                    diff = f"{abs(float(jv) - float(rv)):.6g}"
                f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{spec.name}\t"
                        f"{jv}\t{rv}\t{diff}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare annotation values between json2vcf output "
                    "and a reference VCF.",
    )
    parser.add_argument(
        "--json2vcf", required=True,
        help="Path to json2vcf output VCF (plain or .gz)",
    )
    parser.add_argument(
        "--reference", required=True,
        help="Path to reference VCF (plain or .gz)",
    )
    parser.add_argument(
        "--output",
        help="Write text report to this file (default: stdout)",
    )
    parser.add_argument(
        "--tsv",
        help="Write discordant examples to this TSV file",
    )
    parser.add_argument(
        "--max-examples", type=int, default=20,
        help="Max discordant examples per field (default: 20)",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print progress to stderr",
    )
    args = parser.parse_args()

    t0 = time.time()

    # Pass 1: load json2vcf
    if args.verbose:
        print("Loading json2vcf variants...", file=sys.stderr)
    j2v_dict = load_json2vcf(args.json2vcf, verbose=args.verbose)

    # Pass 2: stream reference, compare
    if args.verbose:
        print("Comparing with reference...", file=sys.stderr)
    stats, total_shared, total_only_j2v, total_only_ref = compare_with_reference(
        j2v_dict, args.reference, FIELD_SPECS, args.max_examples,
        verbose=args.verbose,
    )

    elapsed = time.time() - t0
    if args.verbose:
        print(f"Done in {elapsed:.1f}s", file=sys.stderr)

    # Report
    out = open(args.output, "w") if args.output else None
    try:
        print_report(
            stats, len(j2v_dict), total_shared, total_only_j2v,
            total_only_ref, args.json2vcf, args.reference, FIELD_SPECS,
            args.max_examples, out=out,
        )
    finally:
        if out:
            out.close()

    # TSV
    if args.tsv:
        write_tsv(stats, FIELD_SPECS, args.tsv)
        if args.verbose:
            print(f"TSV written to {args.tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
