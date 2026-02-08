# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Converts Nirvana/Illumina Connected Annotations JSON output to VCF (Variant Call Format) files. Pure Python, zero external dependencies.

## Build & Test Commands

Use `python3` (not `python`) and `python3 -m pytest` (not bare `pytest`).

```bash
pip install -e ".[dev]"        # install in dev mode with pytest
python3 -m pytest -v           # run all 109 tests
python3 -m pytest tests/test_parser.py -v      # run parser tests only
python3 -m pytest tests/test_mapper.py -v      # run mapper tests only
python3 -m pytest tests/test_integration.py -v # run end-to-end tests
python3 -m pytest -k "test_csq"                # run tests matching keyword
```

## CLI Usage

```bash
json2vcf -i input.json.gz -o output.vcf
json2vcf -i input.json -o output.vcf --csq-only      # VEP-style CSQ only
json2vcf -i input.json.gz --no-samples                # omit genotype columns
json2vcf -i input.json.gz --assembly GRCh37           # override assembly
```

## Architecture

Streaming pipeline: **parse** → **map** → **write** (one position at a time, no full-file load).

- `json2vcf/parser.py` — Streams Nirvana's line-based JSON format (`.json` or `.json.gz`). Each position line is parsed independently. Key entry point: `stream_positions(path)` yields `(NirvanaHeader, Position)` tuples.
- `json2vcf/mapper.py` — Transforms `Position` dataclasses into VCF record dicts. Handles per-allele (Number=A) field ordering for multi-allelic sites, VEP-style CSQ construction, ClinVar/gnomAD/SpliceAI mapping, and VCF INFO value escaping.
- `json2vcf/vcf_writer.py` — Writes VCF 4.2 plain text (no PyVCF dependency). Outputs `##` header lines, `#CHROM` column line, and tab-separated data lines.
- `json2vcf/models.py` — Dataclasses (`NirvanaHeader`, `Position`, `Variant`, `Sample`, `TranscriptAnnotation`, `ClinVarEntry`, `PopulationFrequency`, `SpliceAIEntry`) that serve as the contract between parser and mapper.
- `json2vcf/constants.py` — All `##INFO`/`##FORMAT` VCF header definitions, GRCh37/GRCh38 contig maps, CSQ field name list.
- `json2vcf/cli.py` — argparse CLI wiring the pipeline together.

## Nirvana JSON Format

Line-based streaming structure (not standard JSON — one object per line):
- Line 1: `{"header":{...},"positions":[`
- Lines 2–N: one position JSON object per line (comma-separated)
- End: `],"genes":[...]}`

Each position maps to one VCF row. Variants are nested per-alt-allele within each position.

## Key Conventions & Gotchas

- **Multi-allelic sites:** `_get_variant_for_allele()` in mapper.py matches variants by `altAllele` string against `position.altAlleles` order. Number=A INFO fields output `.` for any allele missing a variant.
- **FILTER mapping:** Empty list `[]` → `.` (missing), `["PASS"]` → `PASS`, others → semicolon-joined.
- **INFO value escaping:** `_escape_info_value()` percent-encodes `%`, space, `;`, `=`, `,` (in that order) for ClinVar IDs, significances, etc.
- **FORMAT field inclusion:** Dynamic — a FORMAT field (GT, DP, GQ, AD, etc.) appears only if at least one sample has data for it.
- **Assembly contig prefixes:** GRCh38 uses `chr1`…`chrM`; GRCh37 uses `1`…`MT` (no `chr` prefix).
- **REVEL score parsing:** Can arrive as `{"score": X}` dict or raw float — parser normalizes both.
- **Float formatting:** `.6g` format (6 significant figures, no trailing zeros).

## Test Strategy

All 109 tests are self-contained with embedded test data in `tests/conftest.py` — no external files or network access needed. Test data covers: SNVs, multi-allelic sites, structural variants, missing fields, empty samples, multiple samples, phased genotypes, de novo variants, ClinVar, SpliceAI, REVEL, and failed filters.
