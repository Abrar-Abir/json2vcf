"""End-to-end integration tests: JSON in → VCF out."""

import gzip
import io
import json

import pytest

from json2vcf.cli import main
from json2vcf.mapper import map_position_to_vcf_record
from json2vcf.parser import parse_position_line, stream_positions
from json2vcf.vcf_writer import write_vcf_header, write_vcf_record
from tests.conftest import (
    make_test_header as _make_header,
    MINIMAL_HEADER,
    TWO_SAMPLE_HEADER,
    MINIMAL_POSITION_SNV,
    MULTI_ALLELIC_POSITION,
    SV_POSITION,
    MISSING_QUAL_POSITION,
    TWO_SAMPLE_POSITION,
    build_nirvana_json_string,
)


def _full_pipeline(positions, header_dict=None, **kwargs):
    """Run the full parse→map→write pipeline and return the VCF string."""
    header_dict = header_dict or MINIMAL_HEADER
    header = _make_header(header_dict)

    out = io.StringIO()
    write_vcf_header(out, header, header.genome_assembly, header.samples, **kwargs)

    for pos_dict in positions:
        pos = parse_position_line(json.dumps(pos_dict))
        record = map_position_to_vcf_record(pos, header, **kwargs)
        write_vcf_record(out, record)

    return out.getvalue()


class TestEndToEndMinimalSNV:
    def test_header_present(self):
        """Complete conversion has proper VCF header."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])

        assert vcf.startswith("##fileformat=VCFv4.2\n")
        assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE001\n" in vcf

    def test_single_data_line(self):
        """Single position produces exactly 1 data line."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]

        assert len(data_lines) == 1

    def test_data_line_fields(self):
        """Data line has correct CHROM, POS, ID, REF, ALT, QUAL, FILTER."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert fields[0] == "chr1"          # CHROM
        assert fields[1] == "12345"         # POS
        assert fields[2] == "rs12345"       # ID
        assert fields[3] == "A"             # REF
        assert fields[4] == "T"             # ALT
        assert fields[5] == "200"           # QUAL
        assert fields[6] == "PASS"          # FILTER

    def test_info_contains_annotations(self):
        """INFO field contains expected annotation keys."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        info = data_lines[0].split("\t")[7]

        assert "phyloP=3.5" in info
        assert "gnomAD_AF=0.00012" in info
        assert "oneKG_AF=0.0002" in info
        assert "CLINVAR_SIG=pathogenic" in info
        assert "CSQ=" in info

    def test_sample_column_present(self):
        """Sample column has genotype data."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert len(fields) == 10  # 8 fixed + FORMAT + 1 sample
        assert fields[8].startswith("GT")
        assert fields[9].startswith("0/1")


class TestEndToEndMultiAllelic:
    def test_comma_separated_alt(self):
        """Multi-allelic produces comma-separated ALT."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert fields[4] == "A,GT"

    def test_per_allele_info(self):
        """Per-allele INFO values are comma-separated."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        info = data_lines[0].split("\t")[7]

        assert "gnomAD_AF=0.05,0.001" in info
        assert "phyloP=1.2,-0.5" in info


class TestEndToEndMultiplePositions:
    def test_three_positions(self):
        """Three positions produce three data lines in order."""
        vcf = _full_pipeline(
            [MINIMAL_POSITION_SNV, MULTI_ALLELIC_POSITION, SV_POSITION]
        )
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]

        assert len(data_lines) == 3
        assert data_lines[0].startswith("chr1\t12345")
        assert data_lines[1].startswith("chr2\t48010488")
        assert data_lines[2].startswith("chr7\t55000000")


class TestEndToEndSV:
    def test_sv_info_fields(self):
        """SV position has SVEND, SVLEN, CIPOS, CIEND in INFO."""
        vcf = _full_pipeline([SV_POSITION])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        info = data_lines[0].split("\t")[7]

        assert "SVEND=55050000" in info
        assert "SVLEN=-50000" in info
        assert "CIPOS=-50,50" in info
        assert "CIEND=-100,100" in info


class TestEndToEndMissingFields:
    def test_missing_qual_filter_id(self):
        """Position with missing quality/filter/dbsnp outputs dots."""
        vcf = _full_pipeline([MISSING_QUAL_POSITION])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert fields[2] == "."     # ID
        assert fields[5] == "."     # QUAL
        assert fields[6] == "."     # FILTER


class TestEndToEndTwoSamples:
    def test_two_sample_columns(self):
        """Two samples produce two sample columns."""
        vcf = _full_pipeline(
            [TWO_SAMPLE_POSITION],
            header_dict=TWO_SAMPLE_HEADER,
        )
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        # 8 fixed + FORMAT + 2 samples = 11
        assert len(fields) == 11


class TestEndToEndOutputParseable:
    def test_all_lines_valid(self):
        """Every output line is either a header or has correct tab count."""
        vcf = _full_pipeline(
            [MINIMAL_POSITION_SNV, MULTI_ALLELIC_POSITION, SV_POSITION, MISSING_QUAL_POSITION]
        )
        lines = vcf.strip().split("\n")
        num_samples = len(MINIMAL_HEADER["samples"])

        for line in lines:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                cols = line.split("\t")
                assert cols[0] == "#CHROM"
            else:
                cols = line.split("\t")
                # Should be 8 fixed + FORMAT + num_samples
                expected = 8 + 1 + num_samples
                assert len(cols) == expected, f"Bad column count: {len(cols)} in line: {line[:80]}"


class TestEndToEndViaFile:
    def test_file_round_trip(self, minimal_snv_json_gz, tmp_path):
        """Full file-based round trip: .json.gz → .vcf"""
        out_path = str(tmp_path / "output.vcf")
        main(["-i", minimal_snv_json_gz, "-o", out_path])

        with open(out_path) as f:
            content = f.read()

        assert content.startswith("##fileformat=VCFv4.2")
        data_lines = [l for l in content.strip().split("\n") if not l.startswith("#")]
        assert len(data_lines) == 1
        assert data_lines[0].startswith("chr1\t12345\trs12345\tA\tT")
