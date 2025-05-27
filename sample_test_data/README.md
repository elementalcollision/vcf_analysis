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

## Large VCF for Performance Testing

Some performance tests, such as `tests/test_vcf_comparison_tool.py::test_vcf_comparison_large_files_performance`, require a large VCF file.
This file is not included in the repository due to its size.

To run these tests, you first need the base large VCF file:

- **Base File:** `chr22.1kg.phase3.v5a.vcf.gz`
- **Source:** 1000 Genomes Project (via University of Washington)
- **URL:** `https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr22.1kg.phase3.v5a.vcf.gz`
- **Size:** Approximately 129MB

Place the downloaded file into the `sample_test_data/` directory.

Next, you need to create a "test-ready" version of this file that has the correct contig header for chromosome `22` (matching the `22.fa` reference):

- **Test-Ready File:** `chr22.1kg.phase3.v5a.testready.vcf.gz`
- **Creation Command (run from workspace root after downloading the base file and ensuring `22.fa` and `22.fa.fai` are present):
  ```bash
  CONTIG_LINE="##contig=<ID=22,length=51304566>" # Length from 22.fa.fai
  { gunzip -c sample_test_data/chr22.1kg.phase3.v5a.vcf.gz | grep '^##' | grep -v '^##contig=<ID=GL'; \
    echo "$CONTIG_LINE"; \
    gunzip -c sample_test_data/chr22.1kg.phase3.v5a.vcf.gz | grep -v '^##'; \
  } | bgzip -c > sample_test_data/chr22.1kg.phase3.v5a.testready.vcf.gz && \
  bcftools index sample_test_data/chr22.1kg.phase3.v5a.testready.vcf.gz
  ```

The performance test will use `chr22.1kg.phase3.v5a.testready.vcf.gz` and `22.fa`.
If the `testready` file is not found, the relevant performance tests will be skipped. 