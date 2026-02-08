"""Tests for the JSON-to-VCF field mapping engine."""

import json

import pytest

from json2vcf.mapper import (
    build_csq_string,
    build_info_field,
    build_sample_columns,
    map_position_to_vcf_record,
    _escape_info_value,
)
from json2vcf.parser import parse_position_line, parse_variant, parse_sample
from json2vcf.models import NirvanaHeader
from tests.conftest import (
    MINIMAL_HEADER,
    MINIMAL_POSITION_SNV,
    MULTI_ALLELIC_POSITION,
    REFERENCE_ONLY_POSITION,
    SV_POSITION,
    MISSING_QUAL_POSITION,
    TWO_SAMPLE_POSITION,
    EMPTY_SAMPLE_POSITION,
    FAILED_FILTER_POSITION,
    SPLICE_AI_POSITION,
)


def _make_header():
    return NirvanaHeader(
        annotator=MINIMAL_HEADER["annotator"],
        creation_time=MINIMAL_HEADER["creationTime"],
        genome_assembly=MINIMAL_HEADER["genomeAssembly"],
        schema_version=MINIMAL_HEADER["schemaVersion"],
        data_version=MINIMAL_HEADER["dataVersion"],
        data_sources=MINIMAL_HEADER["dataSources"],
        samples=MINIMAL_HEADER["samples"],
    )


def _parse_pos(position_dict):
    return parse_position_line(json.dumps(position_dict))


class TestBasicMapping:
    def test_snv_chrom_pos_ref_alt(self):
        """Minimal SNV maps to correct CHROM/POS/REF/ALT."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["CHROM"] == "chr1"
        assert record["POS"] == 12345
        assert record["REF"] == "A"
        assert record["ALT"] == "T"

    def test_snv_id_from_dbsnp(self):
        """dbsnp array maps to ID column."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["ID"] == "rs12345"

    def test_missing_id_produces_dot(self):
        """Position with no dbsnp outputs ID='.'"""
        pos = _parse_pos(MISSING_QUAL_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["ID"] == "."

    def test_multi_allelic_alt(self):
        """Two alt alleles produce comma-separated ALT."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["ALT"] == "A,GT"

    def test_multi_allelic_id_dedup(self):
        """dbsnp from multiple variants are deduplicated."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        # Only first variant has rs9999, second has none
        assert record["ID"] == "rs9999"


class TestQualFilter:
    def test_quality_present(self):
        """quality=200.0 outputs QUAL=200."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["QUAL"] == "200"

    def test_quality_missing(self):
        """quality=None outputs QUAL='.'"""
        pos = _parse_pos(MISSING_QUAL_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["QUAL"] == "."

    def test_filter_pass(self):
        """filters=["PASS"] outputs FILTER=PASS."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["FILTER"] == "PASS"

    def test_filter_multiple(self):
        """filters=["LowQual","LowDP"] outputs semicolon-joined."""
        pos = _parse_pos(FAILED_FILTER_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["FILTER"] == "LowQual;LowDP"

    def test_filter_empty(self):
        """filters=[] outputs FILTER='.'"""
        pos = _parse_pos(REFERENCE_ONLY_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        # Empty filter list should produce "."
        assert record["FILTER"] == "."


class TestInfoField:
    def test_info_gnomad_af(self):
        """gnomAD allele frequency appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "gnomAD_AF=0.00012" in info

    def test_info_gnomad_ac_an(self):
        """gnomAD AC and AN appear in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "gnomAD_AC=15" in info
        assert "gnomAD_AN=125000" in info

    def test_info_phylop(self):
        """phyloP score appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "phyloP=3.5" in info

    def test_info_onekg(self):
        """1000 Genomes AF appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "oneKG_AF=0.0002" in info

    def test_info_multi_allelic_per_alt(self):
        """Per-allele fields output comma-separated values in alt order."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        info = build_info_field(pos)

        assert "gnomAD_AF=0.05,0.001" in info

    def test_info_multi_allelic_phylop(self):
        """Per-allele phyloP with values for both alleles."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        info = build_info_field(pos)

        assert "phyloP=1.2,-0.5" in info

    def test_info_clinvar_significance(self):
        """ClinVar significance appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "CLINVAR_SIG=pathogenic" in info
        assert "CLINVAR_ID=RCV000012345" in info

    def test_info_clinvar_review_status_escaped(self):
        """ClinVar review status with special chars is escaped."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        # "criteria provided, single submitter" has a comma and space
        assert "CLINVAR_REVSTAT=criteria%20provided%2C%20single%20submitter" in info

    def test_info_sv_fields(self):
        """SV position outputs SVEND, SVLEN, CIPOS, CIEND."""
        pos = _parse_pos(SV_POSITION)
        info = build_info_field(pos)

        assert "SVEND=55050000" in info
        assert "SVLEN=-50000" in info
        assert "CIPOS=-50,50" in info
        assert "CIEND=-100,100" in info

    def test_info_cytoband(self):
        """CytoBand appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "CytoBand=1p36.33" in info

    def test_info_no_annotations(self):
        """Position with no annotations outputs '.' INFO."""
        pos = _parse_pos(REFERENCE_ONLY_POSITION)
        info = build_info_field(pos)

        assert info == "."

    def test_info_splice_ai(self):
        """SpliceAI scores appear in INFO."""
        pos = _parse_pos(SPLICE_AI_POSITION)
        info = build_info_field(pos)

        assert "SpliceAI_DG_SCORE=0.8" in info
        assert "SpliceAI_DG_DIST=-2" in info
        assert "SpliceAI_AG_SCORE=0.1" in info

    def test_csq_only_mode(self):
        """csq_only=True outputs only CSQ, no flat INFO fields."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos, csq_only=True)

        assert info.startswith("CSQ=")
        assert "gnomAD_AF" not in info
        assert "phyloP=" not in info
        assert "CytoBand" not in info


class TestCSQString:
    def test_csq_format(self):
        """CSQ field contains pipe-delimited transcript annotation."""
        variant_data = MINIMAL_POSITION_SNV["variants"][0]
        variant = parse_variant(variant_data)
        csq = build_csq_string(variant)

        assert csq is not None
        parts = csq.split("|")
        assert parts[0] == "T"               # Allele
        assert parts[1] == "missense_variant" # Consequence
        assert parts[2] == "GENE1"           # SYMBOL
        assert parts[3] == "1234"            # Gene
        assert parts[4] == "Transcript"      # Feature_type
        assert parts[5] == "NM_001234.5"     # Feature
        assert parts[6] == "protein_coding"  # BIOTYPE

    def test_csq_canonical(self):
        """Canonical transcript is marked YES."""
        variant_data = MINIMAL_POSITION_SNV["variants"][0]
        variant = parse_variant(variant_data)
        csq = build_csq_string(variant)

        parts = csq.split("|")
        assert parts[16] == "YES"  # CANONICAL

    def test_csq_no_transcripts(self):
        """Variant without transcripts returns None."""
        variant_data = MISSING_QUAL_POSITION["variants"][0]
        variant = parse_variant(variant_data)
        csq = build_csq_string(variant)

        assert csq is None


class TestSampleColumns:
    def test_format_and_sample(self):
        """Sample data maps to FORMAT string and per-sample values."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        format_str, samples = build_sample_columns(pos.samples)

        assert "GT" in format_str
        assert "DP" in format_str
        assert "GQ" in format_str
        assert "AD" in format_str
        assert "VF" in format_str
        assert len(samples) == 1
        # Sample values: 0/1:40:99:20,20:0.5
        parts = samples[0].split(":")
        assert parts[0] == "0/1"   # GT
        assert "40" in parts       # DP
        assert "99" in parts       # GQ

    def test_two_samples(self):
        """Two samples produce two sample columns."""
        pos = _parse_pos(TWO_SAMPLE_POSITION)
        format_str, samples = build_sample_columns(pos.samples)

        assert len(samples) == 2
        # First sample: 0/1, second: 0/0
        assert samples[0].startswith("0/1")
        assert samples[1].startswith("0/0")

    def test_empty_sample(self):
        """isEmpty sample outputs ./. for GT."""
        pos = _parse_pos(EMPTY_SAMPLE_POSITION)
        format_str, samples = build_sample_columns(pos.samples)

        assert len(samples) == 1
        assert samples[0].startswith("./.")

    def test_sv_sample_format(self):
        """SV sample includes SR and PR in FORMAT."""
        pos = _parse_pos(SV_POSITION)
        format_str, samples = build_sample_columns(pos.samples)

        assert "SR" in format_str
        assert "PR" in format_str

    def test_no_samples(self):
        """No samples returns empty strings."""
        format_str, samples = build_sample_columns(None)

        assert format_str == ""
        assert samples == []


class TestEscaping:
    def test_escape_semicolon(self):
        assert _escape_info_value("a;b") == "a%3Bb"

    def test_escape_equals(self):
        assert _escape_info_value("a=b") == "a%3Db"

    def test_escape_space(self):
        assert _escape_info_value("a b") == "a%20b"

    def test_escape_comma(self):
        assert _escape_info_value("a,b") == "a%2Cb"
