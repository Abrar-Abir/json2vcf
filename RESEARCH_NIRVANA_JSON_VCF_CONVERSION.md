# Research: Nirvana JSON to VCF Conversion Tools & Methods

**Date**: February 8, 2026  
**Research Goal**: Find and document tools, methods, and approaches for converting Nirvana JSON format (from DRAGEN annotation service) to VCF (Variant Call Format) files.

---

## Executive Summary

While there is limited dedicated tooling specifically for Nirvana JSON to VCF conversion, we identified:
- **1 direct converter tool** (eggd_nirvana2vcf)
- **Multiple general-purpose approaches** including custom Python scripts
- **Parsing libraries** for working with both JSON and VCF formats
- **Official Illumina tools** that support JSON output but no longer actively maintained as open source

---

## 1. DEDICATED CONVERSION TOOLS

### 1.1 eggd_nirvana2vcf

**Status**: ✅ Exists (but archived/older)

**Repository**: [eastgenomics/eggd_nirvana2vcf](https://github.com/eastgenomics/eggd_nirvana2vcf)

**Details**:
- **Language**: Python (67%) + Shell (33%)
- **Latest Release**: v1.1.1 (May 12, 2021) - **5 years old**
- **Stars**: 2 (very low adoption)
- **Forks**: 0
- **Contributors**: 3
- **Status**: No recent activity, likely not maintained
- **Last Commit**: 5 years ago

**Description**:
Converts Nirvana JSON into a VEP-style annotated VCF format. This is the direct tool for the conversion task.

**Key Characteristics**:
- Streams JSON input file (memory efficient, doesn't load entire file)
- Outputs VEP-style annotations in VCF format
- Was part of the East Genomics pipeline
- Apache-2.0 licensed (permissive)

**Pros**:
- Direct purpose-built tool for this exact conversion
- Handles streaming (good for large files)
- Open source and available to examine code

**Cons**:
- Minimal maintenance (last update 5 years ago)
- Very low adoption (2 stars, 0 forks)
- No recent issues or pull requests
- Risk of compatibility issues with modern tools
- Limited documentation

**Difficulty Level**: ⭐⭐ (Low - tool exists but may need debugging)

---

## 2. OFFICIAL ILLUMINA TOOLS

### 2.1 Nirvana (Legacy - No Longer Maintained)

**Repository**: [Illumina/Nirvana](https://github.com/Illumina/Nirvana)

**Details**:
- **Language**: C# (99.4%)
- **Latest Release**: v3.18.1 (June 14, 2022) - **3+ years old**
- **Stars**: 190 (reasonable community adoption)
- **Forks**: 50
- **Contributors**: 4
- **Status**: ⛔ **OFFICIALLY NO LONGER MAINTAINED**

**Important Notice** (from repository):
> Nirvana is no longer actively maintained as an open sourced tool. Please visit [Illumina Connected Annotations](https://developer.illumina.com/illumina-connected-annotations) for the latest version.

**Description**:
The original Nirvana annotator that generates JSON output from VCF input. It can optionally output VCF (limited annotation sources) alongside JSON, but the reverse (JSON to VCF) was not its primary function.

**Capabilities**:
- Takes VCF input → generates JSON output
- Can optionally output limited VCF alongside JSON
- Clinical-grade annotation of genomic variants
- Handles SNVs, MNVs, insertions, deletions, indels, and SVs

**Cons**:
- No longer maintained
- Not designed for reverse conversion (VCF from JSON)
- Complex C# codebase, not suitable for simple conversion task

**Difficulty Level**: ⭐⭐⭐⭐⭐ (Very High - would require substantial C# development)

---

### 2.2 Illumina Connected Annotations (Successor - Commercial)

**Website**: [developer.illumina.com/illumina-connected-annotations](https://developer.illumina.com/illumina-connected-annotations)

**Documentation**: [illumina.github.io/IlluminaConnectedAnnotationsDocumentation](https://illumina.github.io/IlluminaConnectedAnnotationsDocumentation/)

**Details**:
- **Status**: ✅ Actively maintained (latest updates 2024-2025)
- **Licensing**: Commercial/Proprietary with free basic tier
- **Support Tiers**: 
  - Basic tier (free) - 17 databases
  - Professional tier (licensed) - adds COSMIC, OMIM, PrimateAI-3D, SpliceAI, PromoterAI
- **Access**: Requires registration at annotation_support@illumina.com

**Capabilities**:
- Modern replacement for Nirvana
- **VCF Input** → **JSON Output** (same direction as Nirvana)
- New **VCF Output support** with limited annotation sources
- Integrates with DRAGEN pipeline
- Supports GRCh37 and GRCh38+ reference genomes
- Clinical-grade annotation

**Relevant for Reverse Conversion**:
⚠️ **Not directly applicable** - this tool also goes VCF → JSON, not JSON → VCF. However, documentation may contain insight into JSON format structure.

**Contact**: annotation_support@illumina.com (for inquiries about JSON to VCF conversion capabilities)

**Difficulty Level**: ⭐⭐ (Low - if proprietary tool has built-in reverse conversion)

---

## 3. PYTHON PARSING LIBRARIES

### 3.1 PyVCF

**Package**: [pyvcf](https://pypi.org/project/PyVCF/)

**Version**: 0.6.8 (Last released: March 18, 2016)

**Details**:
- **Status**: ⚠️ Older, but stable
- **Language**: Python
- **License**: MIT-equivalent
- **Documentation**: http://pyvcf.rtfd.org/

**Capabilities**:
- Parse VCF files (read and write)
- Inspect variant records, samples, genotypes
- Handle VCFv4.0 and v4.1 formats
- Fetch records from tabix-indexed files
- Create VCF writers with custom metadata

**Usage for JSON to VCF**:
This is essential for the **write-side** of conversion - you would use PyVCF to create properly formatted VCF output from parsed JSON data.

**Example**:
```python
import vcf
vcf_writer = vcf.Writer(open('output.vcf', 'w'), vcf_reader)
vcf_writer.write_record(record)
```

**Pros**:
- Well-established, works with standard Python
- Good for writing VCF files
- Handles metadata correctly

**Cons**:
- Last updated 2016
- Doesn't parse Nirvana JSON format
- Need complementary JSON parsing

**Difficulty Level**: ⭐ (Very Low - straightforward library usage)

---

### 3.2 Python Standard JSON Library

**Module**: `json` (built-in)

**Capabilities**:
- Parse JSON files into Python dictionaries
- Simple and efficient for reading Nirvana JSON output

**Example**:
```python
import json
with open('nirvana_output.json', 'r') as f:
    data = json.load(f)
```

**Difficulty Level**: ⭐ (Very Low - built-in, no installation needed)

---

## 4. GENERAL APPROACHES & CUSTOM SCRIPT METHODS

### 4.1 Manual Custom Script (Python-based)

**Approach**: Write custom Python code to:
1. Parse Nirvana JSON using `json` library
2. Extract variant information (CHROM, POS, REF, ALT, etc.)
3. Build VCF record objects using PyVCF
4. Write to VCF file

**Advantages**:
- Full control over mapping logic
- Can handle custom annotation preservation
- Can adapt to specific Nirvana JSON schema variations
- No dependency on unmaintained tools

**Disadvantages**:
- Requires understanding both JSON and VCF formats
- Need to map Nirvana-specific fields to standard VCF
- Testing needed to ensure correctness
- Maintenance responsibility

**Required Knowledge**:
- Nirvana JSON structure (need sample file to analyze)
- VCF specification (columns, metadata, annotation format)
- Python programming

**Difficulty Level**: ⭐⭐⭐ (Medium - requires format knowledge and coding)

---

### 4.2 Scripting Wrapper Around eggd_nirvana2vcf

**Approach**: 
1. Clone and test eggd_nirvana2vcf repository
2. Fix any compatibility issues with modern Python versions
3. Wrap as Docker container or package for distribution

**Advantages**:
- Existing codebase as reference
- Don't need to understand every detail of format mapping
- Can examine what fields are preserved

**Disadvantages**:
- Tool is abandoned - need to maintain fixes
- Unknown if it handles all Nirvana variants
- May not work with newer DRAGEN output formats

**Difficulty Level**: ⭐⭐⭐ (Medium - debugging and adaptation needed)

---

### 4.3 Processing Pipeline Approach

**Approach**: 
1. Use DRAGEN to output both JSON and VCF simultaneously
2. Use existing DRAGEN VCF output as base
3. Enrich with JSON annotations as needed
4. Only convert JSON to VCF when original VCF unavailable

**Advantages**:
- Most reliable if starting data is available
- DRAGEN VCF output is high quality
- Only enhances existing data

**Disadvantages**:
- Requires access to DRAGEN pipeline
- Not a direct JSON-to-VCF solution
- Adds processing complexity

**Difficulty Level**: ⭐⭐ (Low - operational approach rather than coding)

---

## 5. VCF FORMAT RESOURCES

### 5.1 VCF Specification

- **Official Spec**: VCF v4.2 (and v4.1)
- **Available at**: samtools.github.io/hts-specs/
- **Key Sections for Conversion**:
  - Header format (##INFO, ##FORMAT, ##contig)
  - Fixed columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
  - Sample genotype columns (optional)
  - Annotation encoding in INFO field (comma/semicolon-delimited)

### 5.2 VEP Annotation Format

The eggd_nirvana2vcf tool converts to **VEP-style annotations**, meaning:
- Uses CSQ (Consequence) field in INFO
- Pipes (|) separate consequence fields
- Standard format for variant effect prediction annotations

---

## 6. SIMILAR TOOLS & ECOSYSTEM

### 6.1 Existing VCF Conversion Tools (Different Purposes)

**Note**: The following tools exist but convert *between other formats* - not specifically for Nirvana JSON:

- **GTCtoVCF** (Illumina) - Converts GTC/BPM files to VCF
- **VEP plugins** - Various Perl plugins for VEP annotation enhancement
- **vcf2xls, vcf2prot** - Convert VCF to other formats (opposite direction)

These demonstrate ecosystem exists for VCF conversions but most convert *to* VCF from array/other sequencing outputs.

---

## 7. ALTERNATIVE APPROACHES

### 7.1 Contact Illumina Support

**Option**: Reach out directly to Illumina

**Contact**: 
- Email: annotation_support@illumina.com
- Portal: developer.illumina.com
- Purpose: Ask if:
  - JSON to VCF reverse conversion is available in Connected Annotations
  - They have internal tools for this
  - They can provide technical specifications for JSON format

**Difficulty Level**: ⭐ (Very Low - just an email)

---

### 7.2 East Genomics Support

**Option**: Contact original authors of eggd_nirvana2vcf

**Potential Contacts**: 
- @mattgarner
- @Yu-jinKim
- @woook

**Purpose**: Ask for:
- Current status of tool
- Known limitations
- Maintenance plans

**Difficulty Level**: ⭐⭐ (Low - requires finding contact info and reaching out)

---

## 8. COMPARISON TABLE

| Tool/Approach | Purpose | Maintenance | Difficulty | Python | Cost | Notes |
|---|---|---|---|---|---|---|
| **eggd_nirvana2vcf** | Direct JSON→VCF | ⛔ None (5y) | ⭐⭐ | ✅ | Free | Only direct tool; risk of obsolescence |
| **Custom Script** | JSON→VCF via PyVCF | ✅ Your team | ⭐⭐⭐ | ✅ | Free | Full control; best long-term |
| **Illumina Connected** | VCF→JSON (reverse) | ✅ Active | ⭐ | ❌ | $ Licensed | Contact for capabilities |
| **PyVCF Library** | VCF writing | ⚠️ Stable | ⭐ | ✅ | Free | Essential for custom scripts |
| **Wrapper eggd_nirvana2vcf** | JSON→VCF (forked) | ✅ Your team | ⭐⭐⭐ | ✅ | Free | Moderate maintenance burden |
| **Pipeline Approach** | Strategic | ✅ Varies | ⭐⭐ | N/A | Varies | Only if DRAGEN access |

---

## 9. RECOMMENDATIONS

### For Immediate Use (Next 1-3 months):

**Option A** (Fastest): 
1. Try **eggd_nirvana2vcf** first
2. Test with sample DRAGEN JSON files
3. Fix Python 3 compatibility issues if needed
4. Document any issues/workarounds

**Difficulty**: ⭐⭐ | **Time**: 2-4 hours | **Risk**: Medium

---

### For Medium-term (3-6 months):

**Option B** (Recommended): 
1. Analyze Nirvana JSON format (obtain sample file)
2. Build **custom Python conversion script** using PyVCF
3. Handle VEP-style annotation format (CSQ field)
4. Include comprehensive tests
5. Document mapping logic

**Difficulty**: ⭐⭐⭐ | **Time**: 4-8 hours | **Risk**: Low

**Benefits**: 
- Maintainable long-term
- Understandable codebase
- No dependency on abandoned tools
- Team ownership

---

### For Long-term (6+ months):

**Option C** (Strategic): 
1. Contact Illumina about capabilities in **Connected Annotations**
2. Evaluate if reverse conversion is now available
3. If available, migrate to supported tool
4. If not, maintain custom script

**Difficulty**: ⭐ | **Time**: 1-2 hours | **Risk**: Very Low

---

## 10. IMPLEMENTATION CHECKLIST

If building a custom converter:

- [ ] Obtain sample Nirvana JSON file from DRAGEN output
- [ ] Document JSON schema/structure
- [ ] Map JSON fields to VCF columns:
  - [ ] Chromosome
  - [ ] Position
  - [ ] ID (rsID or .)
  - [ ] Reference allele
  - [ ] Alternate alleles
  - [ ] Quality score
  - [ ] Filter status
  - [ ] Annotations (CSQ field for VEP format)
  - [ ] Sample genotypes (if present)
- [ ] Install PyVCF: `pip install pyvcf`
- [ ] Write conversion function
- [ ] Handle edge cases:
  - [ ] Multi-allelic variants
  - [ ] Complex annotations
  - [ ] Missing/null values
- [ ] Test with known variants
- [ ] Validate output against VCF v4.2 spec
- [ ] Document conversion logic
- [ ] Consider performance for large files

---

## 11. REQUIRED INFORMATION FOR NEXT STEPS

To proceed, you will need:

1. **Sample DRAGEN JSON file** - To understand exact structure
2. **Nirvana JSON schema documentation** - Official format specification
3. **Clarification on annotation preservation** - Which fields must be preserved in VCF
4. **Performance requirements** - File size expectations
5. **VCF version target** - VCFv4.1 or v4.2?
6. **Annotation style preference** - VEP CSQ format or standard INFO fields?

---

## 12. REFERENCES & LINKS

### Official Documentation
- Illumina Connected Annotations: https://illumina.github.io/IlluminaConnectedAnnotationsDocumentation/
- VCF Specification: https://samtools.github.io/hts-specs/VCFv4.2.pdf
- VEP Documentation: https://www.ensembl.org/info/docs/tools/vep/index.html

### Tools
- eggd_nirvana2vcf: https://github.com/eastgenomics/eggd_nirvana2vcf
- Nirvana (Legacy): https://github.com/Illumina/Nirvana
- PyVCF: https://pypi.org/project/PyVCF/ | https://github.com/jamiedowell/pyvcf

### Contacts
- Illumina Annotations Support: annotation_support@illumina.com
- East Genomics (eggd): https://github.com/eastgenomics

---

## 13. CONCLUSION

**Direct tooling is scarce but exists**. The only dedicated Nirvana JSON to VCF converter (eggd_nirvana2vcf) is no longer actively maintained but is available and functional. 

**Recommended path**: Evaluate eggd_nirvana2vcf for quick wins, then invest in a custom Python-based converter using PyVCF for long-term maintainability and flexibility. Contact Illumina to verify if their modern Connected Annotations platform has added reverse conversion capabilities.

---

**Document Version**: 1.0  
**Last Updated**: February 8, 2026  
**Status**: Complete Research Phase
