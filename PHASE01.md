# Phase 01 — Smoke Test Results

**Date:** 2026-02-11
**Shard:** `dragen.anno_099.json.gz` (chrY — smallest shard)
**VM:** `advapp@172.32.79.51`
**json2vcf source:** `/home/advapp/json2vcf-src` (installed editable)

---

## Command

```bash
~/.local/bin/json2vcf \
  -i /gpfs/scratch/daily-process/externs/Muiz/joint-calling/1kg-dragen-nirvana/dragen.anno_099.json.gz \
  -o /tmp/test_099.vcf \
  --no-samples
```

---

## Result: PASS

json2vcf successfully converted the chrY shard (239,349 variants) to valid VCF 4.2 with zero errors.

---

## Validation Checks

| # | Check | Result | Detail |
|---|-------|--------|--------|
| 1 | Exit code is 0 | PASS | `EXIT_CODE=0` |
| 2 | Output file exists and non-empty | PASS | `/tmp/test_099.vcf` present |
| 3 | Total lines | PASS | 239,414 lines |
| 4 | Meta-header lines (`##`) | PASS | 64 lines |
| 5 | `#CHROM` line present | PASS | 1 line found |
| 6 | `#CHROM` line has 8 columns | PASS | 8 tab-separated columns (no FORMAT/sample) |
| 7 | Data line count | PASS | 239,349 data lines |
| 8 | All data lines have 8 columns | PASS | Every line has exactly 8 tab-separated fields |
| 9 | All CHROM values are `chrY` | PASS | Only `chrY` found |
| 10 | All INFO keys declared in header | PASS | 30 used keys, all among 37 declared |
| 11 | Spot-check first 5 data lines | PASS | See sample output below |

---

## Output Summary

| Metric | Value |
|--------|-------|
| Total lines | 239,414 |
| Header lines (`##`) | 64 |
| Column header (`#CHROM`) | 1 |
| Data lines | 239,349 |
| Columns per line | 8 (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) |

### Variant Breakdown

| Type | Count | % |
|------|-------|---|
| SNVs | 202,342 | 84.5% |
| Indels | 23,954 | 10.0% |
| Multi-allelic | 13,053 | 5.5% |
| **Total** | **239,349** | **100%** |

### Annotation Coverage

| Annotation | Variants with data | % of total |
|------------|-------------------|------------|
| gnomAD AF | 121,094 | 50.6% |
| CSQ (transcripts) | 72,248 | 30.2% |
| SpliceAI | 90 | <0.1% |
| REVEL | 82 | <0.1% |
| ClinVar | 5 | <0.1% |

### INFO Keys

**Declared in header (37):** CIEND, CIPOS, CLINVAR_ID, CLINVAR_REVSTAT, CLINVAR_SIG, CSQ, CytoBand, DANN, GERP, REVEL, SVEND, SVLEN, SVTYPE, SpliceAI_AG_DIST, SpliceAI_AG_SCORE, SpliceAI_AL_DIST, SpliceAI_AL_SCORE, SpliceAI_DG_DIST, SpliceAI_DG_SCORE, SpliceAI_DL_DIST, SpliceAI_DL_SCORE, TOPMed_AF, gnomAD_AC, gnomAD_AF, gnomAD_AFR_AF, gnomAD_AMR_AF, gnomAD_AN, gnomAD_EAS_AF, gnomAD_EUR_AF, gnomAD_SAS_AF, oneKG_AF, oneKG_AFR_AF, oneKG_AMR_AF, oneKG_EAS_AF, oneKG_EUR_AF, oneKG_SAS_AF, phyloP

**Used in data (30):** CLINVAR_ID, CLINVAR_REVSTAT, CLINVAR_SIG, CSQ, CytoBand, DANN, GERP, REVEL, SpliceAI_AG_DIST, SpliceAI_AG_SCORE, SpliceAI_AL_DIST, SpliceAI_AL_SCORE, SpliceAI_DG_DIST, SpliceAI_DG_SCORE, SpliceAI_DL_DIST, SpliceAI_DL_SCORE, gnomAD_AC, gnomAD_AF, gnomAD_AFR_AF, gnomAD_AMR_AF, gnomAD_AN, gnomAD_EAS_AF, gnomAD_SAS_AF, oneKG_AF, oneKG_AFR_AF, oneKG_AMR_AF, oneKG_EAS_AF, oneKG_EUR_AF, oneKG_SAS_AF, phyloP

**Declared but not used (7):** CIEND, CIPOS, SVEND, SVLEN, SVTYPE, TOPMed_AF, gnomAD_EUR_AF — expected since chrY has no structural variants in this shard and limited annotation coverage.

---

## Error/Warning Check

| Check | Result |
|-------|--------|
| Lines with empty INFO | 0 |
| Lines with INFO=`.` | 0 |
| Lines with malformed INFO (`;;`, leading/trailing `;`) | 0 |
| Lines with empty REF or ALT | 0 |
| Crash or exception | None |

---

## Sample Output (First 5 Data Lines)

```
chrY  2781554  .              G   A   .  .  CytoBand=Yp11.2;phyloP=0.1;DANN=0.84;gnomAD_AF=9.4e-05;gnomAD_AC=1;gnomAD_AN=10643;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_EAS_AF=0;gnomAD_SAS_AF=0;CSQ=A|upstream_gene_variant|RN...
chrY  2781605  .              CA  C   .  .  CytoBand=Yp11.2;CSQ=-|upstream_gene_variant|RNU6-1334P|ENSG00000251841|Transcript|ENST00000516032.1|snRNA||||||||||YES||
chrY  2781622  rs2051119814   C   T   .  .  CytoBand=Yp11.2;phyloP=-2.1;DANN=0.86;gnomAD_AF=5.7e-05;gnomAD_AC=1;gnomAD_AN=17499;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_EAS_AF=0;gnomAD_SAS_AF=0;CSQ=T|upstream_gen...
chrY  2781625  .              A   G   .  .  CytoBand=Yp11.2;phyloP=0.1;DANN=0.07;gnomAD_AF=0.000286;gnomAD_AC=5;gnomAD_AN=17480;gnomAD_AFR_AF=0.000468;gnomAD_AMR_AF=0;gnomAD_EAS_AF=0;gnomAD_SAS_AF=0;CSQ=G|upstream_gene_va...
chrY  2781639  .              A   G   .  .  CytoBand=Yp11.2;phyloP=0.1;DANN=0.07;gnomAD_AF=0.000296;gnomAD_AC=5;gnomAD_AN=16872;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0.003077;gnomAD_EAS_AF=0;gnomAD_SAS_AF=0;CSQ=G|upstream_gene_va...
```

(INFO fields truncated for readability — full lines have valid semicolon-separated key=value pairs.)

---

## VCF Header

```
##fileformat=VCFv4.2
##source=json2vcf (from Nirvana 3.18.1)
##INFO=<ID=CSQ,...>
... (64 meta-header lines total, including 37 INFO definitions and 25 contig lines)
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
```

---

## Issues Found

**None.** json2vcf processed all 239,349 positions in shard 099 without errors. Every output line is well-formed VCF with correctly declared INFO keys.

---

## Conclusion

**Phase 01: PASS**

json2vcf successfully converts a real Nirvana JSON shard to syntactically valid VCF 4.2 output. The tool handles annotation-only data (no samples) correctly, produces proper 8-column output, and all INFO field keys are properly declared in the header. Ready to proceed to Phase 2 (position-level concordance).
