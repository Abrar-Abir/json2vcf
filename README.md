# json2vcf

Convert [Nirvana/Illumina Connected Annotations](https://illumina.github.io/NirvanaDocumentation/) JSON output to VCF 4.2 format.

## Features

- Pure Python, zero external dependencies
- Streaming pipeline — processes one position at a time, no full-file load
- Reads `.json` and `.json.gz` input
- Supports GRCh37 and GRCh38 assemblies (auto-detected from header)
- VEP-style CSQ field with per-transcript annotations
- Annotations: gnomAD, ClinVar, SpliceAI, REVEL, DANN, GERP, phyloP, 1000 Genomes, TOPMed

## Installation

```bash
pip install -e .
```

## Usage

```bash
# Basic conversion
json2vcf -i input.json.gz -o output.vcf

# VEP-style CSQ only (no flat INFO fields)
json2vcf -i input.json -o output.vcf --csq-only

# Omit sample/genotype columns
json2vcf -i input.json.gz -o output.vcf --no-samples

# Override genome assembly
json2vcf -i input.json.gz -o output.vcf --assembly GRCh37

# Output to stdout
json2vcf -i input.json.gz
```

## Development

```bash
pip install -e ".[dev]"
python3 -m pytest -v
```

## Architecture

Streaming pipeline: **parse** → **map** → **write**

- `json2vcf/parser.py` — Streams Nirvana's line-based JSON format, yielding `(NirvanaHeader, Position)` tuples
- `json2vcf/mapper.py` — Transforms positions into VCF record dicts (per-allele fields, CSQ, INFO escaping)
- `json2vcf/vcf_writer.py` — Writes VCF 4.2 plain text
- `json2vcf/models.py` — Dataclass contracts between parser and mapper
- `json2vcf/constants.py` — VCF header definitions, contig maps, CSQ field names
