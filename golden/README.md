# Golden File Validation System

## Overview

The golden file validation system ensures that VCF processing outputs remain consistent across code changes and catch regressions in future development. Golden files contain known-good reference outputs from validated VCF processing operations.

## Directory Structure

```
golden/
├── bcftools/           # bcftools command outputs
│   ├── view/           # bcftools view header outputs
│   ├── stats/          # bcftools stats statistical summaries
│   ├── query/          # bcftools query formatted outputs
│   ├── filter/         # bcftools filter filtered results
│   └── norm/           # bcftools norm normalized outputs
├── validation/         # VCF validation results
│   ├── valid_files/    # Validation reports for valid VCF files
│   ├── invalid_files/  # Error reports for invalid VCF files
│   └── multi_file/     # Multi-file validation summaries
├── comparison/         # VCF comparison outputs
│   ├── basic/          # Basic VCF comparison results
│   ├── normalized/     # Normalized comparison results
│   └── multi_sample/   # Multi-sample comparison summaries
└── analysis/           # Agent analysis outputs
    ├── summary/        # VCF summary analysis reports
    ├── detailed/       # Detailed variant analysis results
    └── graph_integration/ # Graph database integration outputs
```

## Golden File Naming Convention

Golden files follow the pattern: `{operation}_{test_file}_{parameters}.{ext}`

Examples:
- `stats_minimal_vcf.txt` - bcftools stats output for minimal.vcf
- `view_header_sample1_vcf.txt` - bcftools view header for sample1.vcf
- `validation_invalid_malformed_vcf.json` - validation result for malformed VCF

## File Formats

- **Text outputs**: `.txt` files for bcftools command outputs
- **JSON outputs**: `.json` files for structured validation results
- **Analysis outputs**: `.md` files for agent analysis reports

## Usage

Golden files are used in automated tests to validate that current outputs match expected results:

```python
def test_bcftools_stats_golden():
    result = run_bcftools_stats("minimal.vcf")
    expected = load_golden_file("bcftools/stats/stats_minimal_vcf.txt")
    assert_golden_match(result, expected)
```

## Maintenance

### Generating Golden Files
1. Run validated operations on known-good test data
2. Manually review outputs for correctness
3. Store in appropriate directory with proper naming
4. Update test cases to use golden files

### Updating Golden Files
1. Review changes carefully when tests fail
2. Verify that changes are intentional improvements
3. Update golden files only after thorough validation
4. Document reasons for updates in commit messages

### Version Control
- Golden files are tracked in Git
- Large files may use Git LFS if needed
- Changes require code review and approval

## Test Integration

Golden file tests are integrated into the main test suite:
- Unit tests: Individual operation validation
- Integration tests: Workflow output validation
- E2E tests: Complete CLI output validation

## Quality Assurance

### Deterministic Outputs
- Timestamps are normalized or excluded
- Random elements are controlled or filtered
- Platform-specific differences are handled

### Tolerance Handling
- Floating-point comparisons use appropriate tolerance
- Non-critical variations are filtered out
- Strict matching for critical data integrity

### Performance
- Golden file comparisons are optimized for speed
- Large outputs are chunked or summarized
- CI/CD integration minimizes overhead

## Troubleshooting

### Test Failures
1. Check if output changes are intentional
2. Verify test data hasn't changed
3. Review for platform-specific differences
4. Update golden files if changes are valid

### Missing Golden Files
1. Generate missing files using validated operations
2. Review outputs manually before committing
3. Add corresponding test cases
4. Document generation process

---

*This system ensures consistent, high-quality VCF processing outputs across all development phases.* 