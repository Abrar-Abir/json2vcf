"""Maps Nirvana Position objects to VCF record dicts."""

from typing import Any, Dict, List, Optional, Tuple

from .constants import CSQ_FIELDS
from .models import (
    NirvanaHeader,
    PopulationFrequency,
    Position,
    Sample,
    Variant,
)


def _escape_info_value(value: str) -> str:
    """Percent-encode special characters in VCF INFO values."""
    return (
        value.replace("%", "%25")
        .replace(" ", "%20")
        .replace(";", "%3B")
        .replace("=", "%3D")
        .replace(",", "%2C")
    )


def _fmt_float(value: float) -> str:
    """Format a float with 6 significant figures."""
    return f"{value:.6g}"


def _get_variant_for_allele(
    variants: Optional[List[Variant]], alt_allele: str
) -> Optional[Variant]:
    """Find the variant matching a specific alt allele."""
    if not variants:
        return None
    for v in variants:
        if v.alt_allele == alt_allele:
            return v
    return None


def _per_allele_values(
    position: Position, extractor, formatter=str
) -> Optional[str]:
    """Build a comma-separated per-allele value string for Number=A INFO fields.

    extractor: callable that takes a Variant and returns a value or None.
    Returns None if all values are missing.
    """
    values = []
    all_missing = True
    for alt in position.alt_alleles:
        variant = _get_variant_for_allele(position.variants, alt)
        val = extractor(variant) if variant else None
        if val is not None:
            values.append(formatter(val))
            all_missing = False
        else:
            values.append(".")
    if all_missing:
        return None
    return ",".join(values)


def build_csq_string(variant: Variant) -> Optional[str]:
    """Build VEP-style CSQ string for one variant (all transcripts).

    Returns pipe-delimited transcript annotations separated by commas
    for multiple transcripts.
    """
    if not variant.transcripts:
        return None

    csq_parts = []
    for t in variant.transcripts:
        fields = [
            variant.alt_allele,                                   # Allele
            "&".join(t.consequence) if t.consequence else "",     # Consequence
            t.hgnc or "",                                         # SYMBOL
            t.gene_id or "",                                      # Gene
            "Transcript",                                         # Feature_type
            t.transcript or "",                                   # Feature
            t.bio_type or "",                                     # BIOTYPE
            t.exons or "",                                        # EXON
            t.introns or "",                                      # INTRON
            t.hgvsc or "",                                        # HGVSc
            t.hgvsp or "",                                        # HGVSp
            t.cdna_pos or "",                                     # cDNA_position
            t.cds_pos or "",                                      # CDS_position
            t.protein_pos or "",                                  # Protein_position
            t.amino_acids or "",                                  # Amino_acids
            t.codons or "",                                       # Codons
            "YES" if t.is_canonical else "",                      # CANONICAL
            _format_prediction_score(t.poly_phen_prediction, t.poly_phen_score),  # PolyPhen
            _format_prediction_score(t.sift_prediction, t.sift_score),        # SIFT
        ]
        csq_parts.append("|".join(fields))

    return ",".join(csq_parts)


def _format_prediction_score(prediction: Optional[str], score: Optional[float]) -> str:
    """Format a prediction tool result (PolyPhen, SIFT) as 'prediction(score)'."""
    if prediction is None:
        return ""
    result = prediction
    if score is not None:
        result += f"({score})"
    return result


def _popfreq_extractor(source_attr: str, freq_attr: str):
    """Return an extractor for a population frequency field on a Variant."""
    def extract(v):
        source = getattr(v, source_attr, None)
        return getattr(source, freq_attr, None) if source else None
    return extract


_POPFREQ_FIELDS = [
    ("gnomAD_AF",     "gnomad", "all_af", _fmt_float),
    ("gnomAD_AC",     "gnomad", "all_ac", str),
    ("gnomAD_AN",     "gnomad", "all_an", str),
    ("gnomAD_AFR_AF", "gnomad", "afr_af", _fmt_float),
    ("gnomAD_AMR_AF", "gnomad", "amr_af", _fmt_float),
    ("gnomAD_EUR_AF", "gnomad", "eur_af", _fmt_float),
    ("gnomAD_EAS_AF", "gnomad", "eas_af", _fmt_float),
    ("gnomAD_SAS_AF", "gnomad", "sas_af", _fmt_float),
    ("oneKG_AF",      "one_kg", "all_af", _fmt_float),
    ("oneKG_AFR_AF",  "one_kg", "afr_af", _fmt_float),
    ("oneKG_AMR_AF",  "one_kg", "amr_af", _fmt_float),
    ("oneKG_EUR_AF",  "one_kg", "eur_af", _fmt_float),
    ("oneKG_EAS_AF",  "one_kg", "eas_af", _fmt_float),
    ("oneKG_SAS_AF",  "one_kg", "sas_af", _fmt_float),
    ("TOPMed_AF",     "topmed", "all_af", _fmt_float),
]


def build_info_field(position: Position, csq_only: bool = False) -> str:
    """Build the INFO column string from position + variant annotations."""
    parts = []

    if not csq_only:
        # Position-level fields
        if position.cytogenetic_band:
            parts.append(f"CytoBand={_escape_info_value(position.cytogenetic_band)}")

        # SV fields (position-level)
        if position.sv_end is not None:
            parts.append(f"SVEND={position.sv_end}")
        if position.sv_length is not None:
            parts.append(f"SVLEN={position.sv_length}")
        if position.ci_pos is not None:
            parts.append(f"CIPOS={','.join(str(x) for x in position.ci_pos)}")
        if position.ci_end is not None:
            parts.append(f"CIEND={','.join(str(x) for x in position.ci_end)}")

        # Per-allele fields from variants
        _add_per_allele(parts, position, "phyloP", lambda v: v.phylop_score, _fmt_float)
        _add_per_allele(parts, position, "DANN", lambda v: v.dann_score, _fmt_float)
        _add_per_allele(parts, position, "GERP", lambda v: v.gerp_score, _fmt_float)
        _add_per_allele(parts, position, "REVEL", lambda v: v.revel_score, _fmt_float)

        # SVTYPE per-allele
        _add_per_allele(
            parts, position, "SVTYPE",
            lambda v: v.variant_type if v.is_structural_variant else None,
            str,
        )

        # Population frequencies (gnomAD, 1000 Genomes, TOPMed)
        for info_key, source_attr, freq_attr, formatter in _POPFREQ_FIELDS:
            _add_per_allele(parts, position, info_key,
                            _popfreq_extractor(source_attr, freq_attr), formatter)

        # ClinVar (from first variant that has it — variable number)
        _add_clinvar_info(parts, position)

        # SpliceAI (from first variant that has it)
        _add_splice_ai_info(parts, position)

    # CSQ (always included)
    csq_strings = []
    if position.variants:
        for variant in position.variants:
            csq = build_csq_string(variant)
            if csq:
                csq_strings.append(csq)
    if csq_strings:
        parts.append(f"CSQ={','.join(csq_strings)}")

    return ";".join(parts) if parts else "."


def _add_per_allele(
    parts: list, position: Position, key: str, extractor, formatter
) -> None:
    """Add a per-allele INFO field if any allele has a value."""
    val = _per_allele_values(position, extractor, formatter)
    if val is not None:
        parts.append(f"{key}={val}")


def _add_clinvar_info(parts: list, position: Position) -> None:
    """Add ClinVar INFO fields, collecting from all variants."""
    if not position.variants:
        return

    sigs = []
    ids = []
    revstats = []
    for variant in position.variants:
        if not variant.clinvar:
            continue
        for cv in variant.clinvar:
            if cv.significance:
                sigs.extend(cv.significance)
            if cv.id:
                ids.append(cv.id)
            if cv.review_status:
                revstats.append(cv.review_status)

    if ids:
        parts.append(f"CLINVAR_ID={_escape_info_value('&'.join(ids))}")
    if sigs:
        parts.append(f"CLINVAR_SIG={_escape_info_value('&'.join(sigs))}")
    if revstats:
        parts.append(f"CLINVAR_REVSTAT={_escape_info_value('&'.join(revstats))}")


_SPLICE_AI_FIELDS = [
    ("SpliceAI_AG_SCORE", "acceptor_gain_score", _fmt_float),
    ("SpliceAI_AG_DIST", "acceptor_gain_distance", str),
    ("SpliceAI_AL_SCORE", "acceptor_loss_score", _fmt_float),
    ("SpliceAI_AL_DIST", "acceptor_loss_distance", str),
    ("SpliceAI_DG_SCORE", "donor_gain_score", _fmt_float),
    ("SpliceAI_DG_DIST", "donor_gain_distance", str),
    ("SpliceAI_DL_SCORE", "donor_loss_score", _fmt_float),
    ("SpliceAI_DL_DIST", "donor_loss_distance", str),
]


def _add_splice_ai_info(parts: list, position: Position) -> None:
    """Add SpliceAI INFO fields from variants."""
    if not position.variants:
        return

    for variant in position.variants:
        if not variant.splice_ai:
            continue
        for sai in variant.splice_ai:
            for info_key, attr_name, formatter in _SPLICE_AI_FIELDS:
                val = getattr(sai, attr_name)
                if val is not None:
                    parts.append(f"{info_key}={formatter(val)}")


def build_sample_columns(
    samples: Optional[List[Sample]],
) -> Tuple[str, List[str]]:
    """Build FORMAT string and per-sample value strings.

    Returns (format_string, [sample1_values, sample2_values, ...]).
    FORMAT fields are included dynamically based on which fields have data.
    GT is always first.
    """
    if not samples:
        return "", []

    # Determine which FORMAT fields are present across all samples
    format_keys = ["GT"]  # GT always first
    field_extractors = _get_format_extractors()

    for key, extractor in field_extractors:
        for sample in samples:
            val = extractor(sample)
            if val is not None:
                format_keys.append(key)
                break

    # Build per-sample value strings
    sample_strings = []
    for sample in samples:
        values = []
        for key in format_keys:
            if key == "GT":
                values.append(sample.genotype)
            else:
                extractor = dict(field_extractors)[key]
                val = extractor(sample)
                if val is None:
                    values.append(".")
                else:
                    values.append(val)
        sample_strings.append(":".join(values))

    return ":".join(format_keys), sample_strings


def _scalar_fmt(attr: str, formatter=str):
    """Create extractor for a scalar sample attribute."""
    def extract(s):
        val = getattr(s, attr)
        return formatter(val) if val is not None else None
    return extract


def _list_fmt(attr: str, formatter=str):
    """Create extractor for a list sample attribute (comma-joined)."""
    def extract(s):
        val = getattr(s, attr)
        return ",".join(formatter(x) for x in val) if val else None
    return extract


def _get_format_extractors() -> List[Tuple[str, Any]]:
    """Return list of (FORMAT_key, extractor_function) pairs."""
    return [
        ("DP", _scalar_fmt("total_depth")),
        ("GQ", _scalar_fmt("genotype_quality")),
        ("AD", _list_fmt("allele_depths")),
        ("VF", _list_fmt("variant_frequencies", _fmt_float)),
        ("CN", _scalar_fmt("copy_number")),
        ("FT", lambda s: ("FAIL" if s.failed_filter else "PASS") if s.failed_filter is not None else None),
        ("DN", lambda s: ("true" if s.is_de_novo else "false") if s.is_de_novo is not None else None),
        ("DQ", _scalar_fmt("de_novo_quality", _fmt_float)),
        ("SR", _list_fmt("split_read_counts")),
        ("PR", _list_fmt("paired_end_read_counts")),
        ("SQ", _scalar_fmt("somatic_quality", _fmt_float)),
    ]


def map_position_to_vcf_record(
    position: Position,
    header: NirvanaHeader,
    csq_only: bool = False,
    include_samples: bool = True,
) -> Dict[str, Any]:
    """Convert a Position into a VCF record dict.

    Returns dict with keys: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO,
    and optionally FORMAT, samples.
    """
    # ID: collect dbsnp from all variants, deduplicate preserving order
    ids = list(dict.fromkeys(
        rsid
        for v in (position.variants or [])
        if v.dbsnp
        for rsid in v.dbsnp
    ))
    id_str = ";".join(ids) if ids else "."

    # ALT
    alt_str = ",".join(position.alt_alleles) if position.alt_alleles else "."

    # QUAL
    if position.quality is not None:
        qual_str = _fmt_float(position.quality) if position.quality != int(position.quality) else str(int(position.quality))
    else:
        qual_str = "."

    # FILTER
    if position.filters is None:
        filter_str = "."
    elif position.filters == [] or position.filters == ["."]:
        filter_str = "."
    elif position.filters == ["PASS"]:
        filter_str = "PASS"
    else:
        filter_str = ";".join(position.filters)

    # INFO
    info_str = build_info_field(position, csq_only=csq_only)

    record = {
        "CHROM": position.chromosome,
        "POS": position.position,
        "ID": id_str,
        "REF": position.ref_allele,
        "ALT": alt_str,
        "QUAL": qual_str,
        "FILTER": filter_str,
        "INFO": info_str,
    }

    # FORMAT + sample columns
    if include_samples and position.samples:
        format_str, sample_strs = build_sample_columns(position.samples)
        record["FORMAT"] = format_str
        record["samples"] = sample_strs
    else:
        record["FORMAT"] = ""
        record["samples"] = []

    return record
