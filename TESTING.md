# TESTING.md — Real-Data Validation Plan for json2vcf

## Goal

Validate `json2vcf` against real 1000 Genomes (1kG) data: convert Nirvana JSON shards to VCF, then compare against an independently-annotated VCF produced by an in-house pipeline on the same variant call set.

---

## Remote VM Access

```
SSH: ssh advapp@172.32.79.51
```

No `bcftools` or `module` system is available on this VM. Use `zcat`, `gzip`, `awk`, `head`, `tail`, `wc`, `python3` for all file operations.

---

## Data Inventory

### Nirvana JSON Shards (Input to json2vcf)

```
Location: /gpfs/scratch/daily-process/externs/Muiz/joint-calling/1kg-dragen-nirvana/
Files:    dragen.anno_001.json.gz  through  dragen.anno_099.json.gz  (99 shards)
          Plus 99 .jsi index files and a few ancillary files (204 total).
Format:   Nirvana line-based JSON (.json.gz)
Assembly: GRCh38
Annotator: Nirvana 3.18.1
Samples:  NONE — annotation-only (no genotype/sample data in JSON)
```

**Shard-to-chromosome mapping** (determined by reading position lines):

| Shards       | Chromosome |
|-------------|-----------|
| 001–008     | chr1      |
| 009–015     | chr2      |
| 016         | EMPTY/ERROR (skip this shard) |
| 017–022     | chr3      |
| 023–028     | chr4      |
| 029–034     | chr5      |
| 035–039     | chr6      |
| 040–044     | chr7      |
| 045–049     | chr8      |
| 050–053     | chr9      |
| 054–057     | chr10     |
| 058–061     | chr11     |
| 062–065     | chr12     |
| 066–069     | chr13     |
| 070–073     | chr14     |
| 074–076     | chr15     |
| 077–079     | chr16     |
| 080–082     | chr17     |
| 083–085     | chr18     |
| 086–087     | chr19     |
| 088–089     | chr20     |
| 090–091     | chr21     |
| 092–093     | chr22     |
| 094–098     | chrX      |
| 099         | chrY      |

**JSON structure** (each position line):
```json
{
  "chromosome": "chr1",
  "position": 10178,
  "refAllele": "CCTAA",
  "altAlleles": ["C"],
  "cytogeneticBand": "1p36.33",
  "variants": [{
    "vid": "1-10178-CCTAA-C",
    "chromosome": "chr1",
    "begin": 10178,
    "end": 10183,
    "refAllele": "CCTAA",
    "altAllele": "C",
    "variantType": "deletion",
    "hgvsg": "NC_000001.11:g.10179_10182del",
    "inLowComplexityRegion": true,
    "dbsnp": ["rs775809821"],
    "gnomad": {"allAf": 0.423, ...},
    "transcripts": {"refSeq": [...], "ensembl": [...]},
    ...
  }]
}
```

### Reference VCF (Comparison Target)

```
Location: /gpfs/scratch/daily-process/externs/Muiz/joint-calling/3202_samples_cohort_annotation_only.vcf.gz
Size:     ~14 GB compressed
Format:   VCFv4.2
Assembly: GRCh38 (chr-prefixed contigs)
Samples:  NONE — annotation-only (8 columns: CHROM through INFO, no FORMAT/sample columns)
```

**In-house pipeline annotations** (INFO fields in the reference VCF):

| INFO Field | Source | Description |
|---|---|---|
| `ANN` | snpEff | Functional annotation (pipe-separated subfields) |
| `AC`, `AN`, `NS`, etc. | DRAGEN joint-caller | Allele counts, sample counts |
| `gnomAD_Exome_AF` | gnomAD v4.0.0 exome | Exome allele frequency |
| `gnomadWGS_AF` | gnomAD v4.0.0 WGS | Whole-genome allele frequency |
| `gnomadWGS_AF_afr`, `_amr`, `_eas`, `_nfe`, `_sas` | gnomAD v4.0.0 WGS | Population-specific WGS AF |
| `CLNSIG` | ClinVar 20240528 | Clinical significance |
| `CLNDN` | ClinVar 20240528 | Disease name |
| `CLNREVSTAT` | ClinVar 20240528 | Review status |
| `CLNVI` | ClinVar 20240528 | Variant ID |
| `CADD` | CADD | Combined annotation-dependent depletion score |
| `rs_ids` | dbSNP 151 | dbSNP rsIDs |
| `dbNSFP_SIFT_pred` | dbNSFP | SIFT prediction |
| `dbNSFP_Polyphen2_HDIV_pred` | dbNSFP | PolyPhen2 prediction |
| `dbNSFP_MutationTaster_pred` | dbNSFP | MutationTaster prediction |
| `dbNSFP_REVEL_score` | dbNSFP | REVEL score |
| `TOPMED` | TOPMed | TOPMed allele frequency |
| `REPEATMASKER` | RepeatMasker | Repeat region flags |
| `LOF`, `NMD` | snpEff | Loss-of-function / NMD predictions |

**Verified position match** (first 4 positions are identical in both files):
```
chr1:10178  CCTAA → C
chr1:10230  AC    → A
chr1:10327  T     → C
chr1:10352  T     → TA
```

### Other Files on VM (NOT Useful for Comparison)

| Path | Why Not |
|---|---|
| `1kg-dragen-anno/` subdirectories (10-annotate through 13-annotate) | All VCFs are EMPTY (0 data lines, headers only — failed pipeline runs) |
| `3202_samples_cohort_mergedsnpEff.vcf.gz` (1.3 TB) | Full multi-sample VCF — too large, has genotype columns |
| `3202_samples_cohort_mergeddbNFSPAnno.vcf.gz` (1.3 TB) | Full multi-sample VCF — too large |
| `60samples_joint_calling/` | Different dataset (60 samples, not same cohort) |
| `1kg-3202-splitted-by-sample/` | Per-sample VCFs, only chr13/18/19/Y |
| `1kg-3202-unanno-splitted-by-sample/` | Unannotated — no annotations to compare |

---

## Annotation Field Mapping (Nirvana → In-House)

These annotations come from **different pipelines** but draw from overlapping databases. Not all fields have a direct counterpart.

### Fields With Direct Counterparts (Comparable)

| Concept | json2vcf Output (Nirvana) | In-House VCF | Match Type |
|---|---|---|---|
| Variant position | `CHROM`, `POS`, `REF`, `ALT` | `CHROM`, `POS`, `REF`, `ALT` | Exact |
| dbSNP ID | `ID` column | `rs_ids` INFO field | Exact (same rsID) |
| gnomAD WGS AF (all) | `gnomAD_AF` | `gnomadWGS_AF` | Numeric (may differ by gnomAD version) |
| gnomAD WGS AF (AFR) | `gnomAD_AFR_AF` | `gnomadWGS_AF_afr` | Numeric |
| gnomAD WGS AF (AMR) | `gnomAD_AMR_AF` | `gnomadWGS_AF_amr` | Numeric |
| gnomAD WGS AF (EAS) | `gnomAD_EAS_AF` | `gnomadWGS_AF_eas` | Numeric |
| gnomAD WGS AF (EUR/NFE) | `gnomAD_EUR_AF` | `gnomadWGS_AF_nfe` | Numeric |
| gnomAD WGS AF (SAS) | `gnomAD_SAS_AF` | `gnomadWGS_AF_sas` | Numeric |
| ClinVar significance | `CLINVAR_SIG` | `CLNSIG` | Categorical (wording may differ) |
| ClinVar review status | `CLINVAR_REVSTAT` | `CLNREVSTAT` | Categorical |
| REVEL score | `REVEL` | `dbNSFP_REVEL_score` | Numeric (float) |
| SIFT prediction | CSQ `SIFT` subfield | `dbNSFP_SIFT_pred` | Categorical |
| PolyPhen prediction | CSQ `PolyPhen` subfield | `dbNSFP_Polyphen2_HDIV_pred` | Categorical |
| Gene symbol | CSQ `SYMBOL` subfield | ANN gene subfield | Exact |
| Consequence | CSQ `Consequence` subfield | ANN consequence subfield | Categorical (different ontology terms possible) |
| TOPMed AF | `TOPMed_AF` | `TOPMED` | Numeric |

### Fields Only in json2vcf Output (No In-House Counterpart)

- `phyloP`, `DANN`, `GERP` (conservation/pathogenicity scores)
- `gnomAD_AC`, `gnomAD_AN` (allele counts — in-house has `AC`/`AN` from DRAGEN, not gnomAD)
- `oneKG_*` (1000 Genomes population frequencies)
- `SpliceAI_*` (splice site predictions)
- `CytoBand`
- `CLINVAR_ID` (accession ID)
- Full CSQ transcript-level detail (exon/intron numbers, HGVS, codons, amino acids)

### Fields Only in In-House VCF (No json2vcf Counterpart)

- `ANN` (snpEff full annotation string — json2vcf uses CSQ instead)
- `CADD` (json2vcf has DANN, not CADD)
- `gnomAD_Exome_AF` (json2vcf outputs gnomAD WGS, not exome)
- `AC`, `AN`, `NS`, `HWE`, `ExcHet`, etc. (cohort-level statistics from DRAGEN)
- `LOF`, `NMD` (snpEff-specific loss-of-function predictions)
- `REPEATMASKER`
- `dbNSFP_MutationTaster_pred` and many other dbNSFP sub-scores

---

## Phase 1 — Smoke Test (Single Small Shard)

### Objective
Verify json2vcf can process a real Nirvana JSON shard without crashing, and produces syntactically valid VCF output.

### Steps

1. **Copy a small shard to the local machine:**
   ```bash
   scp advapp@172.32.79.51:/gpfs/scratch/daily-process/externs/Muiz/joint-calling/1kg-dragen-nirvana/dragen.anno_099.json.gz /tmp/dragen.anno_099.json.gz
   ```
   Shard 099 = chrY (smallest shard). If too large, alternatively use `head` on the remote to extract first ~500 position lines from shard 001 into a mini JSON file.

2. **Run json2vcf:**
   ```bash
   json2vcf -i /tmp/dragen.anno_099.json.gz -o /tmp/test_099.vcf --no-samples
   ```
   Use `--no-samples` because the JSON has no sample data.

3. **Validate VCF output:**
   - Check exit code is 0
   - Count output lines: `wc -l /tmp/test_099.vcf`
   - Check header is present: `grep '^##' /tmp/test_099.vcf | wc -l`
   - Check `#CHROM` line is present and has 8 columns (no FORMAT/sample since `--no-samples`)
   - Check data lines have exactly 8 tab-separated columns
   - Verify all chromosomes in data lines are `chrY`
   - Spot-check: print first 5 data lines and visually inspect REF/ALT/INFO

4. **Check for errors/warnings:**
   - If json2vcf crashes or raises exceptions, capture the traceback
   - Note any positions that produce empty or malformed INFO strings

### Success Criteria
- json2vcf exits with code 0
- Output has valid VCF header + `#CHROM` line + data lines
- Every data line has 8 tab-separated fields
- Chromosomes are all `chrY`
- INFO field keys are all declared in the header

### Potential Issues to Watch For
- The JSON `header` line format: Nirvana's first line is `{"header":{...},"positions":[` — parser must handle the trailing `,"positions":[`
- Shard 016 is known to be empty/broken — skip it
- The JSON has no `samples` field — json2vcf must handle `samples: None` gracefully
- Large shards may take significant time — monitor progress

---

## Phase 2 — Position-Level Concordance

### Objective
Verify that json2vcf produces the same set of variants (CHROM, POS, REF, ALT) as the reference VCF for a given genomic region.

### Steps

1. **Convert 2–3 shards spanning different chromosomes.** Good candidates:
   - Shard 099 (chrY — small, ~single shard)
   - Shard 090 (chr21 — small chromosome)
   - Shard 001 (chr1 first chunk — large, complex)

   ```bash
   json2vcf -i /tmp/dragen.anno_099.json.gz -o /tmp/json2vcf_099.vcf --no-samples
   json2vcf -i /tmp/dragen.anno_090.json.gz -o /tmp/json2vcf_090.vcf --no-samples
   json2vcf -i /tmp/dragen.anno_001.json.gz -o /tmp/json2vcf_001.vcf --no-samples
   ```

2. **Extract corresponding regions from the reference VCF on the remote VM.**
   Since `bcftools` is unavailable, use streaming extraction:
   ```bash
   # On remote VM — extract chrY lines
   ssh advapp@172.32.79.51 'zcat /gpfs/scratch/daily-process/externs/Muiz/joint-calling/3202_samples_cohort_annotation_only.vcf.gz | awk -F"\t" "\$1==\"chrY\"" | gzip > /tmp/ref_chrY.vcf.gz'
   scp advapp@172.32.79.51:/tmp/ref_chrY.vcf.gz /tmp/ref_chrY.vcf.gz
   ```
   **WARNING:** The reference VCF is 14 GB compressed. Streaming the whole file with `zcat` takes a long time. For chromosomes near the end of the file (chrY), this could take 30+ minutes. Strategies:
   - Run the extraction as a background job on the VM (`nohup ... &`)
   - Alternatively, extract only the first/last N positions for a quick check
   - Alternatively, use Python on the VM to seek more efficiently

3. **Compare variant sets** with a Python script:
   ```python
   # compare_positions.py
   # Reads two VCFs, extracts (CHROM, POS, REF, ALT) tuples, computes set overlap
   # Outputs: shared count, only-in-json2vcf count, only-in-reference count
   # Plus: first 20 examples from each exclusive set
   ```

### Success Criteria
- **>99% overlap** of `(CHROM, POS, REF, ALT)` tuples between json2vcf output and reference VCF for the same chromosome
- Any discrepancies should be explainable (e.g., allele normalization differences, multi-allelic splitting differences)

### Potential Issues to Watch For
- **Allele representation:** Nirvana may left-align or normalize differently than the in-house pipeline. A deletion might be `POS=100 REF=AC ALT=A` in one and `POS=101 REF=C ALT=.` (as a deletion) in the other.
- **Multi-allelic splitting:** Nirvana JSON has variants nested per-allele within a position. json2vcf groups them back into a single VCF row. The in-house VCF may or may not split multi-allelics into separate rows.
- **Shard boundary variants:** Variants at chromosome boundaries between shards should not be duplicated or lost.

---

## Phase 3 — Annotation Concordance

### Objective
Compare annotation values that both pipelines derive from the same underlying databases: gnomAD, ClinVar, dbSNP, REVEL, SIFT, PolyPhen, TOPMed.

### Steps

1. **Write a Python comparison tool** (`tests/compare_annotations.py` or `scripts/compare_vcfs.py`). The tool should:

   **a) Parse both VCFs streaming (line-by-line):**
   - Index by `(CHROM, POS, REF, ALT)` key
   - For multi-allelic rows, expand to per-allele entries
   - Handle both gzipped and plain text input

   **b) Compare these annotation pairs:**

   | Comparison | json2vcf Field | Reference Field | Method |
   |---|---|---|---|
   | gnomAD AF (all) | `gnomAD_AF` | `gnomadWGS_AF` | Numeric: abs diff < 0.01 or relative diff < 5% |
   | gnomAD AF (AFR) | `gnomAD_AFR_AF` | `gnomadWGS_AF_afr` | Numeric |
   | gnomAD AF (AMR) | `gnomAD_AMR_AF` | `gnomadWGS_AF_amr` | Numeric |
   | gnomAD AF (EAS) | `gnomAD_EAS_AF` | `gnomadWGS_AF_eas` | Numeric |
   | gnomAD AF (EUR/NFE) | `gnomAD_EUR_AF` | `gnomadWGS_AF_nfe` | Numeric |
   | gnomAD AF (SAS) | `gnomAD_SAS_AF` | `gnomadWGS_AF_sas` | Numeric |
   | ClinVar significance | `CLINVAR_SIG` | `CLNSIG` | Categorical exact match (after normalizing case/underscores) |
   | ClinVar review status | `CLINVAR_REVSTAT` | `CLNREVSTAT` | Categorical exact match (after normalizing) |
   | dbSNP rsID | `ID` column | `rs_ids` INFO field | Set overlap (variant may have multiple rsIDs) |
   | REVEL score | `REVEL` | `dbNSFP_REVEL_score` | Numeric: abs diff < 0.01 |
   | TOPMed AF | `TOPMed_AF` | `TOPMED` | Numeric: abs diff < 0.01 |
   | Gene symbol | CSQ `SYMBOL` | ANN gene field | Exact match per transcript |
   | SIFT | CSQ `SIFT` | `dbNSFP_SIFT_pred` | Categorical (D=deleterious, T=tolerated) |
   | PolyPhen | CSQ `PolyPhen` | `dbNSFP_Polyphen2_HDIV_pred` | Categorical (D/P/B) |

   **c) Produce a concordance report:**

   ```
   === ANNOTATION CONCORDANCE REPORT ===
   Total variants compared: 1,234,567

   POSITION MATCH
     Shared:          1,234,000 (99.95%)
     Only in json2vcf:      200 (0.02%)
     Only in reference:     367 (0.03%)

   gnomAD_AF vs gnomadWGS_AF
     Both present:    800,000
     Concordant:      795,000 (99.4%)
     Discordant:        5,000 (0.6%)
     Only in json2vcf:  10,000
     Only in reference:  5,000
     Mean abs diff:   0.00023
     Max abs diff:    0.15
     Correlation:     0.9998

   CLINVAR_SIG vs CLNSIG
     Both present:    12,000
     Exact match:     11,500 (95.8%)
     Mismatch:           500 (4.2%)
     Top mismatches:
       "Pathogenic" vs "Pathogenic/Likely_pathogenic": 120
       ...

   [... similar blocks for each field pair ...]

   === DISCORDANT EXAMPLES (first 20) ===
   chr1:12345 A>T
     gnomAD_AF: json2vcf=0.0012  ref=0.0010  (diff=0.0002)
     CLNSIG:    json2vcf=pathogenic  ref=Pathogenic  (case diff only)
   ...
   ```

   **d) Optionally produce TSV output for downstream analysis:**
   ```
   CHROM  POS  REF  ALT  FIELD  JSON2VCF_VALUE  REFERENCE_VALUE  MATCH  DIFF
   chr1   12345  A   T   gnomAD_AF  0.0012  0.0010  CONCORDANT  0.0002
   ...
   ```

2. **Run on chrY first** (smallest), then expand to chr21, then chr1 shard 001.

3. **Interpret results:**
   - **gnomAD AF:** Expect high concordance if both use gnomAD v4. If Nirvana bundles a different gnomAD version, expect systematic small differences across all variants.
   - **ClinVar:** Expect high concordance for significance, but exact wording may differ (e.g., `pathogenic` vs `Pathogenic`, `Likely_pathogenic` vs `likely_pathogenic`). Normalize before comparing.
   - **dbSNP rsID:** Should be exact match. Differences indicate version mismatch.
   - **REVEL:** Nirvana may parse REVEL differently (dict `{"score": X}` vs raw float). Check that json2vcf normalizes correctly.

### Success Criteria
- gnomAD AF correlation > 0.99 (differences are version-related, not conversion errors)
- ClinVar significance concordance > 90% (after normalization)
- dbSNP rsID concordance > 99%
- REVEL concordance > 95% where both have values
- No systematic conversion errors (e.g., all AFs being 0, or all ClinVar being empty)

### Key Implementation Notes for the Comparison Tool

**Parsing json2vcf output INFO field:**
```
INFO fields are semicolon-separated: key=value;key=value
Number=A fields: comma-separated per alt allele
CSQ: Allele|Consequence|SYMBOL|...,Allele|Consequence|SYMBOL|...
     (comma separates transcript entries; pipe separates subfields)
```

**Parsing reference VCF INFO field:**
```
gnomadWGS_AF=0.123        (single value)
CLNSIG=Pathogenic          (single value)
rs_ids=rs12345             (may have commas for multiple)
ANN=T|missense_variant|MODERATE|GENE|...  (pipe-separated, comma for multi-transcript)
dbNSFP_REVEL_score=0.85   (single value)
TOPMED=0.123               (single value)
```

**Multi-allelic handling:**
- json2vcf: single row with comma-separated ALTs and comma-separated Number=A values
- Reference VCF: may be single row or split into biallelic rows — need to handle both

**CSQ parsing (for gene/consequence comparison):**
CSQ field order (19 subfields, pipe-separated):
```
Allele|Consequence|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|CANONICAL|PolyPhen|SIFT
```

**ANN parsing (snpEff format, for gene/consequence comparison):**
```
Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos/cDNA.length|CDS.pos/CDS.length|AA.pos/AA.length|Distance|ERRORS_WARNINGS_INFO
```

---

## Phase 4 — Full-Scale Validation (Optional)

### Objective
Convert all 99 shards and do a genome-wide comparison.

### Steps

1. **Convert all shards (on the remote VM to avoid transfer):**
   - Install json2vcf on the VM: `pip install --user -e .` (after scp-ing the repo)
   - Or: `scp` the json2vcf package to the VM, run in a loop:
   ```bash
   for i in $(seq -w 1 99); do
     [ "$i" = "016" ] && continue  # skip broken shard
     json2vcf -i /gpfs/scratch/.../dragen.anno_${i}.json.gz \
              -o /tmp/json2vcf_${i}.vcf --no-samples
   done
   ```

2. **Concatenate into a single VCF:**
   ```bash
   # Take header from first output
   grep '^#' /tmp/json2vcf_001.vcf > /tmp/json2vcf_all.vcf
   # Append data lines from all shards (in order)
   for i in $(seq -w 1 99); do
     [ "$i" = "016" ] && continue
     grep -v '^#' /tmp/json2vcf_${i}.vcf >> /tmp/json2vcf_all.vcf
   done
   gzip /tmp/json2vcf_all.vcf
   ```

3. **Run the comparison tool genome-wide:**
   ```bash
   python3 compare_vcfs.py \
     --json2vcf /tmp/json2vcf_all.vcf.gz \
     --reference /gpfs/scratch/.../3202_samples_cohort_annotation_only.vcf.gz \
     --output /tmp/concordance_report.txt \
     --tsv /tmp/discordant_records.tsv
   ```

4. **Generate summary statistics and plots** (if matplotlib available):
   - Scatter: gnomAD AF (json2vcf) vs gnomAD AF (reference)
   - Histogram: AF difference distribution
   - Bar chart: variant type distribution in both
   - ClinVar confusion matrix
   - Per-chromosome concordance rates

### Disk Space Considerations
- 99 VCFs uncompressed could be very large (the 14 GB gzipped reference VCF is probably 100+ GB uncompressed)
- Write VCFs to `/tmp/` or a scratch directory on the VM
- Clean up intermediate files after comparison
- Consider streaming comparison (never fully materialize both VCFs) to avoid disk issues

---

## Comparison Tool Design Specification

### Input
```
compare_vcfs.py --json2vcf <path> --reference <path> [--output <report.txt>] [--tsv <records.tsv>] [--max-records N] [--chromosomes chr1,chr2,...]
```

### Architecture
Streaming, memory-efficient:
1. Read both VCFs simultaneously, sorted by (CHROM, POS)
2. Use a merge-join approach: advance the file that's "behind" until positions align
3. When positions match, compare annotations
4. Accumulate statistics in fixed-size counters (no per-variant storage in memory)
5. Collect first N discordant examples for the report

### Dependencies
- Pure Python (to match the zero-dependency philosophy of json2vcf)
- `gzip` module for `.gz` files
- Optional: `matplotlib` for plots (skip if not available)

### Output Files
1. **Console/text report**: Human-readable concordance summary (see Phase 3 example)
2. **TSV file**: Machine-readable discordant records for downstream analysis
3. **Plots** (optional): PNG files if matplotlib is available

---

## Quick-Reference: json2vcf CLI for This Test

```bash
# Basic conversion (annotation-only, no samples)
json2vcf -i input.json.gz -o output.vcf --no-samples

# CSQ-only mode (only VEP-style CSQ in INFO, no flat fields)
json2vcf -i input.json.gz -o output.vcf --no-samples --csq-only

# Force assembly (normally auto-detected from JSON header)
json2vcf -i input.json.gz -o output.vcf --no-samples --assembly GRCh38

# Output to stdout (pipe to gzip)
json2vcf -i input.json.gz --no-samples | gzip > output.vcf.gz
```

---

## Quick-Reference: Remote VM File Extraction

Since `bcftools` is not available, use these patterns:

```bash
# Extract VCF header from gzipped file
ssh advapp@172.32.79.51 "zcat FILE.vcf.gz | head -500 | grep '^##'"

# Extract #CHROM line
ssh advapp@172.32.79.51 "zcat FILE.vcf.gz | grep '^#CHROM'"

# Extract data lines for a specific chromosome (SLOW for late chromosomes)
ssh advapp@172.32.79.51 "zcat FILE.vcf.gz | awk -F'\t' '\$1==\"chrY\"' | head -100"

# Count data lines for a chromosome
ssh advapp@172.32.79.51 "zcat FILE.vcf.gz | awk -F'\t' '\$1==\"chrY\"' | wc -l"

# Extract first N data lines (any chromosome)
ssh advapp@172.32.79.51 "zcat FILE.vcf.gz | grep -v '^#' | head -1000"

# Count total data lines
ssh advapp@172.32.79.51 "zcat FILE.vcf.gz | grep -v '^#' | wc -l"
```

**WARNING:** The reference VCF is ~14 GB gzipped. `zcat` on the full file takes 30+ minutes. For testing, prefer:
- Extracting early chromosomes (chr1 data starts at the beginning of the file)
- Using `head` with large line counts to grab a prefix
- Running long extractions as background jobs (`nohup ... &`)

---

## Summary of Phases

| Phase | Scope | Time Estimate | Validates |
|---|---|---|---|
| 1: Smoke Test | 1 shard (chrY) | Quick | json2vcf runs, VCF is well-formed |
| 2: Position Concordance | 2–3 shards | Medium | Same variants in both files |
| 3: Annotation Concordance | 2–3 shards | Medium-Long | Annotations agree where expected |
| 4: Full-Scale | All 99 shards | Long | Genome-wide validation |
