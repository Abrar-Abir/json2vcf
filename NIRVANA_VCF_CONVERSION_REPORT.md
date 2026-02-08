# Nirvana JSON to VCF Conversion: Test Data & Validation Report

**Date:** February 8, 2026  
**Topic:** Finding test JSON-VCF pairs and validation methods for Nirvana format conversion

---

## Executive Summary

There are **multiple sources of test data and validation approaches** available for verifying Nirvana JSON to VCF conversion. While dedicated test pairs are limited, the Nirvana official repository includes test data, and comprehensive VCF validation tools exist to ensure output correctness.

---

## 1. Test Data Availability

### 1.1 Nirvana Official Resources

**Nirvana GitHub Repository**
- **URL:** https://github.com/Illumina/Nirvana
- **Contents:** Includes a `UnitTests` directory with test cases
- **Value:** Real test data for integration testing and format validation
- **Limitation:** Tests are primarily for VCF→JSON conversion, not reverse

**Illumina Connected Annotations** (Modern successor to Nirvana)
- **URL:** https://developer.illumina.com/illumina-connected-annotations
- **Documentation:** https://illumina.github.io/IlluminaConnectedAnnotationsDocumentation/
- **Status:** Active, maintained product
- **Caveat:** May require contact for JSON→VCF conversion capabilities

### 1.2 Public Reference Test Data

| Source | Type | Use Case |
|--------|------|----------|
| **1000 Genomes Project** | Real-world VCF files | Baseline test data, annotation examples |
| **GATK Resources** (https://github.com/broadinstitute/gatk) | Test variants in `src/test/resources/` | Complex variant scenarios |
| **ClinVar** | Clinical variants in VCF | Real annotation examples |
| **gnomAD** | Population variants | Large-scale test data |

### 1.3 Creating Custom Test Data

**Recommended Approach:**
1. Start with **minimal test cases** (3-5 simple variants)
2. Progress to **medium complexity** (multiple alleles, mixed types)
3. Add **edge cases** (missing values, structural variants, non-standard chromosomes)

**Minimal Test Case Example:**
- Single SNV variant
- Multi-allelic variant
- Variant with missing QUAL
- Variant with complex annotations

---

## 2. VCF Format Validation Tools

### 2.1 Professional VCF Validators

#### **vcf-validator** (EBI Variation Archive) ⭐ Recommended
- **Repository:** https://github.com/EBIvariation/vcf-validator
- **Installation:** `conda install -c bioconda vcf-validator`
- **Language:** C++14
- **Status:** Actively maintained

**Features:**
- Validates against VCF 4.1, 4.2, 4.3, 4.4 specifications
- Three validation levels:
  - `error` - Syntax only
  - `warning` - Full compliance checks
  - `stop` - Halt on first error
- Two main tools:
  - `vcf_validator` - Syntax and semantic validation
  - `vcf_assembly_checker` - Validates REF alleles against FASTA reference
- Multiple report formats: `summary`, `text`, `valid`

**Example Usage:**
```bash
vcf_validator -i output.vcf -l warning -r summary,text
vcf_assembly_checker -i output.vcf -f reference.fa
```

#### **HTSlib/BCFtools** (Samtools) ⭐ Standard Industry Tool
- **HTSlib Repository:** https://github.com/samtools/htslib
- **BCFtools Repository:** https://github.com/samtools/bcftools
- **Language:** C/Perl
- **Status:** Production-grade, widely used

**Features:**
- Fast VCF parsing and validation
- VCF merging, sorting, filtering
- Variant comparison capabilities
- Used in standard bioinformatics pipelines

**Example Usage:**
```bash
bcftools view -l 0 file.vcf
bcftools isec file1.vcf file2.vcf  # Compare two VCF files
```

#### **VCFtools**
- **Repository:** https://github.com/vcftools/vcftools
- **Website:** https://vcftools.github.io/
- **Language:** C++/Perl
- **Features:** Basic validation, filtering, comparison tools

---

### 2.2 Python-Based Validation

#### **PyVCF**
- **Repository:** https://github.com/jamescasbon/PyVCF
- **Language:** Python
- **VCF Support:** v4.0 and v4.1

**Capabilities:**
- Record-by-record parsing and validation
- Easy access to metadata, filters, INFO fields, FORMAT fields
- Custom validation script writing
- Integration with test frameworks

**Example:**
```python
import vcf

vcf_reader = vcf.Reader(filename='file.vcf')
for record in vcf_reader:
    assert record.CHROM is not None
    assert record.POS > 0
    assert record.REF is not None
```

---

## 3. VCF Comparison Methods

### 3.1 Tool-Based Comparison

**bcftools isec** (Recommended)
```bash
bcftools isec file1.vcf file2.vcf
# Output: Shows intersection, unique to file1, unique to file2
```

**VCFtools --diff**
```bash
vcftools --vcf original.vcf --diff converted.vcf --out comparison
# Output: Detailed comparison report
```

### 3.2 Custom Python Comparison

**Script-Based Approach:**
```python
import vcf

reader1 = vcf.Reader(filename='expected.vcf')
reader2 = vcf.Reader(filename='output.vcf')

differences = []
for i, (rec1, rec2) in enumerate(zip(reader1, reader2)):
    if rec1.CHROM != rec2.CHROM:
        differences.append(f"Record {i}: CHROM mismatch")
    if rec1.POS != rec2.POS:
        differences.append(f"Record {i}: POS mismatch")
    if rec1.REF != rec2.REF:
        differences.append(f"Record {i}: REF mismatch")
    if str(rec1.ALT) != str(rec2.ALT):
        differences.append(f"Record {i}: ALT mismatch")

if differences:
    print("Differences found:")
    for diff in differences:
        print(f"  - {diff}")
else:
    print("Files match!")
```

---

## 4. VCF Format Specification

### 4.1 Official Standards

| Version | URL | Status |
|---------|-----|--------|
| VCF v4.5 (latest) | http://samtools.github.io/hts-specs/VCFv4.5.pdf | Current |
| VCF v4.4 | http://samtools.github.io/hts-specs/VCFv4.4.pdf | Supported |
| VCF v4.3 | http://samtools.github.io/hts-specs/VCFv4.3.pdf | Supported |
| VCF v4.2 | http://samtools.github.io/hts-specs/VCFv4.2.pdf | Supported |
| VCF v4.1 | http://samtools.github.io/hts-specs/VCFv4.1.pdf | Supported |

**Maintained by:** SAM/BAM HTS Specs https://samtools.github.io/hts-specs/

---

## 5. Nirvana JSON Schema

### 5.1 Documentation Sources

- **Nirvana GitHub:** https://github.com/Illumina/Nirvana (with documentation)
- **Connected Annotations Docs:** https://illumina.github.io/IlluminaConnectedAnnotationsDocumentation/

### 5.2 Typical Nirvana JSON Structure

```json
{
  "header": {
    "annotationVersion": "1.0",
    "dataVersion": "2021-01-01"
  },
  "genes": [ ... ],
  "transcripts": [ ... ],
  "variants": [
    {
      "chromosome": "chr1",
      "position": 12345,
      "refAllele": "A",
      "altAlleles": ["T"],
      "variantId": "rs123",
      "variantType": "SNV",
      "phylopScore": 0.5,
      "annotationSet": [ ... ]
    }
  ]
}
```

### 5.3 Key Fields to Map for VCF Conversion

| Nirvana Field | VCF Field | Notes |
|---------------|-----------|-------|
| `chromosome` | CHROM | Must include "chr" prefix |
| `position` | POS | 1-based coordinate |
| `refAllele` | REF | Single or multi-base |
| `altAlleles` | ALT | Comma-separated list |
| `variantId` | ID | rs-number or "." if absent |
| `variantType` | INFO | Can store as custom field |
| Annotation fields | INFO | Must define in header ##INFO lines |

---

## 6. Validation Testing Strategy

### 6.1 Key Validation Checkpoints

#### **Header Validation**
- [ ] VCF format version specified (##fileformat=VCFv4.2)
- [ ] Contig lines match reference genome
- [ ] INFO field definitions present for all INFO fields used
- [ ] FORMAT field definitions present for genotypes (if applicable)
- [ ] Required ##INFO lines for Nirvana annotations

#### **Record-Level Validation**
- [ ] CHROM field valid (chr1-22, chrX, chrY, chrM, or contig name)
- [ ] POS is positive integer
- [ ] REF allele matches reference (if validating against FASTA)
- [ ] ALT alleles properly formatted (no invalid characters)
- [ ] QUAL field valid (numeric or ".")
- [ ] FILTER field valid (PASS, . or defined filter ID)

#### **Annotation Integrity**
- [ ] All Nirvana fields mapped to INFO
- [ ] No loss of numeric precision (float formatting)
- [ ] Nested annotations properly encoded (semicolons, commas, etc.)
- [ ] Empty/null values handled correctly

#### **Sample Data Preservation** (if multi-sample)
- [ ] Genotypes (GT) maintained
- [ ] All FORMAT fields complete
- [ ] Phasing information preserved (if present)
- [ ] Sample count matches input

---

### 6.2 Automated Testing Pipeline

**Recommended Workflow:**

```bash
# Step 1: Run basic format validation
vcf_validator -i output.vcf -l warning -r summary

# Step 2: Check assembly compatibility (if reference available)
vcf_assembly_checker -i output.vcf -f reference.fa

# Step 3: Compare with baseline/expected output
bcftools isec expected.vcf output.vcf

# Step 4: Detailed comparison
vcftools --vcf expected.vcf --diff output.vcf --out comparison_report
```

---

### 6.3 Testing Phases

| Phase | Tests | Tools |
|-------|-------|-------|
| **Phase 1: Unit** | Single variant conversion | PyVCF, custom Python |
| **Phase 2: Format** | VCF compliance | vcf-validator |
| **Phase 3: Comparison** | Variant matching | bcftools, vcftools |
| **Phase 4: Annotation** | Field preservation | Custom Python scripts |
| **Phase 5: Integration** | Real Nirvana outputs | All tools combined |
| **Phase 6: Edge Cases** | Missing values, complex variants | Custom test cases |

---

## 7. Recommended Implementation Approach

### Step-by-Step Plan

1. **Gather Test Data**
   - Download sample Nirvana outputs from official GitHub or create minimal cases
   - Keep baseline VCF files for comparison

2. **Set Up Validation Stack**
   - Install: `conda install -c bioconda vcf-validator bcftools vcftools`
   - Install Python: `pip install PyVCF`

3. **Create Test Suite**
   - Unit tests: Single variant conversions
   - Integration tests: Full Nirvana JSON files
   - Regression tests: Known edge cases

4. **Implement Custom Validation**
   - Python script using PyVCF to verify annotations
   - Record-by-record comparison
   - Field integrity checks

5. **Continuous Validation**
   - Run vcf-validator on all outputs
   - Compare with baseline using bcftools
   - Generate reports for each conversion batch

---

## 8. Tools Comparison Table

| Tool | Purpose | Language | Pros | Cons | Use Case |
|------|---------|----------|------|------|----------|
| **vcf-validator** | Format validation | C++ | Comprehensive, EBI standard | May require installation | Primary format check |
| **bcftools** | VCF manipulation | C/Perl | Fast, widely adopted | Steeper learning curve | Comparison, merging |
| **PyVCF** | Python parsing | Python | Easy integration, flexible | Limited to v4.1 | Custom validation scripts |
| **VCFtools** | Analysis/filtering | C++/Perl | Good comparison tools | Less modern | Comparison workflows |

---

## 9. Resources & References

### Official Documentation
- VCF Specification: https://samtools.github.io/hts-specs/
- Nirvana Repo: https://github.com/Illumina/Nirvana
- Connected Annotations: https://developer.illumina.com/illumina-connected-annotations

### Tool Documentation
- vcf-validator: https://github.com/EBIvariation/vcf-validator
- BCFtools: http://samtools.github.io/bcftools/
- PyVCF: https://github.com/jamescasbon/PyVCF

### Test Data Sources
- 1000 Genomes: https://www.internationalgenome.org/
- GATK: https://github.com/broadinstitute/gatk
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/

---

## 10. Conclusion

**Available Resources:**
- ✅ Test data available in Nirvana GitHub repository
- ✅ Multiple professional VCF validators ready to use
- ✅ Comparison tools available (bcftools, VCFtools)
- ✅ Python ecosystem supports custom validation

**Recommended Approach:**
1. Start with EBI's `vcf-validator` for format checking
2. Use `bcftools` for variant comparison
3. Write custom PyVCF scripts for annotation-specific validation
4. Build test cases incrementally from simple to complex

**Key Advantage:**
Combining multiple validation tools provides comprehensive coverage of both VCF format compliance and annotation integrity preservation.
