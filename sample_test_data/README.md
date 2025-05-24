# Canonical Edge Case VCFs – sample_test_data/

This directory contains canonical, version-controlled edge case VCF files for robust testing, normalization, and validation of VCF processing tools.

## Purpose
- Provide a reproducible set of minimal, strictly formatted VCFs covering common and rare edge cases encountered in genomics workflows.
- Ensure all files are bgzipped, tabix-indexed, and validated with bcftools.
- Support regression testing, normalization, and compliance validation.

## Files
- **minimal.vcf.gz**: Minimal valid VCF (single variant, standard header)
- **multiallelic.vcf.gz**: Site with multiple alternate alleles
- **empty_alt.vcf.gz**: Record with an empty ALT field
- **inconsistent_format.vcf.gz**: FORMAT field inconsistent with data columns
- **bad_info.vcf.gz**: Malformed INFO field
- **symbolic_allele.vcf.gz**: Symbolic ALT allele (e.g., <DEL>)
- **nonstandard_chrom.vcf.gz**: Nonstandard chromosome name (e.g., 'chrUn')
- **missing_header_field.vcf.gz**: Missing INFO header for a field used in the body (e.g., AF)

All files are bgzipped and have a corresponding `.tbi` index.

## Validation
To validate all files with bcftools:
```bash
for f in sample_test_data/*.vcf.gz; do bcftools view "$f"; done
```

## Usage
- Use these files for normalization, comparison, and compliance tool testing.
- Add new edge cases as needed, ensuring strict formatting and bcftools validation.

## Notes
- These files are intended for automated and manual regression testing.
- See the project README for more details on edge case testing strategy. 