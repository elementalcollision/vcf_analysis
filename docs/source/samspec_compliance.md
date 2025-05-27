# SAMspec Compliance Validation

The VCF Analysis Agent includes comprehensive SAMspec compliance validation to ensure VCF files conform to the SAM/VCF specification standards (VCF 4.2/4.3).

## Overview

SAMspec compliance validation performs thorough checks of VCF files against the official specification, detecting:

- **Critical violations**: Specification violations that break compatibility
- **Warnings**: Recommended practices not followed
- **Info**: Minor issues or suggestions for improvement

## Features

### Comprehensive Validation Rules

The validator checks over 30 different compliance rules including:

#### File Format Validation
- `MISSING_FILEFORMAT`: Missing required ##fileformat header line
- `INVALID_FILEFORMAT`: Invalid ##fileformat line format
- `DUPLICATE_FILEFORMAT`: Multiple ##fileformat lines found
- `FILEFORMAT_NOT_FIRST`: ##fileformat should be the first line
- `UNSUPPORTED_VERSION`: VCF version may not be fully supported

#### Header Structure Validation
- `MISSING_CHROM_LINE`: Missing required #CHROM header line
- `INVALID_CHROM_LINE`: #CHROM line has insufficient fields
- `INVALID_CHROM_FIELD`: Invalid field in #CHROM line
- `MISSING_FORMAT_FIELD`: FORMAT field required when sample columns present
- `INVALID_FORMAT_FIELD`: Expected 'FORMAT' field not found

#### INFO Field Validation
- `INVALID_INFO_FORMAT`: Invalid ##INFO format
- `MISSING_INFO_FIELD`: Missing required field in INFO definition
- `INVALID_INFO_TYPE`: Invalid Type for INFO field
- `INVALID_INFO_NUMBER`: Invalid Number for INFO field
- `UNDEFINED_INFO`: INFO field not defined in header

#### FORMAT Field Validation
- `INVALID_FORMAT_FORMAT`: Invalid ##FORMAT format
- `MISSING_FORMAT_FIELD`: Missing required field in FORMAT definition
- `INVALID_FORMAT_TYPE`: Invalid Type for FORMAT field
- `INVALID_FORMAT_NUMBER`: Invalid Number for FORMAT field
- `UNDEFINED_FORMAT`: FORMAT field not defined in header

#### FILTER Field Validation
- `INVALID_FILTER_FORMAT`: Invalid ##FILTER format
- `MISSING_FILTER_FIELD`: Missing required field in FILTER definition
- `MISSING_PASS_FILTER`: Missing PASS filter definition
- `MISSING_FILTER_DESCRIPTION`: FILTER missing description
- `UNDEFINED_FILTER`: Filter not defined in header

#### Data Record Validation
- `INSUFFICIENT_FIELDS`: Data record has insufficient fields
- `INVALID_CHROM`: CHROM field cannot be empty or '.'
- `INVALID_POS_FORMAT`: POS field must be an integer
- `INVALID_POS`: POS field must be >= 1
- `INVALID_ID_FORMAT`: ID field contains invalid characters
- `INVALID_REF`: REF field cannot be empty or '.'
- `INVALID_REF_BASES`: REF field contains invalid bases
- `INVALID_ALT_BASES`: ALT field contains invalid allele
- `INVALID_QUAL_FORMAT`: QUAL field must be a number or '.'
- `NEGATIVE_QUAL`: QUAL field is negative
- `MISSING_SAMPLE_DATA`: FORMAT field present but no sample data
- `FORMAT_SAMPLE_MISMATCH`: Sample data doesn't match FORMAT fields

#### Contig Validation
- `INVALID_CONTIG_FORMAT`: Invalid ##contig format
- `MISSING_CONTIG_ID`: Missing ID field in contig definition
- `MISSING_CONTIG_LENGTH`: Contig missing length field

### File Format Support

- **Uncompressed VCF files** (`.vcf`)
- **Gzip-compressed VCF files** (`.vcf.gz`)
- **VCF versions**: 4.0, 4.1, 4.2, 4.3

## CLI Usage

### Basic Validation

```bash
# Validate a single VCF file
vcf-agent samspec validate sample.vcf

# Validate with verbose output
vcf-agent samspec validate sample.vcf --verbose

# Validate with quiet output (summary only)
vcf-agent samspec validate sample.vcf --quiet

# Treat warnings as failures (strict mode)
vcf-agent samspec validate sample.vcf --strict
```

### Output Formats

```bash
# Text output (default)
vcf-agent samspec validate sample.vcf --format text

# JSON output
vcf-agent samspec validate sample.vcf --format json

# Save report to file
vcf-agent samspec validate sample.vcf --output report.txt
vcf-agent samspec validate sample.vcf --format json --output report.json
```

### Batch Validation

```bash
# Validate multiple files
vcf-agent samspec batch-validate file1.vcf file2.vcf file3.vcf

# Batch validation with output directory
vcf-agent samspec batch-validate *.vcf --output-dir reports/

# Batch validation with summary report
vcf-agent samspec batch-validate *.vcf --summary

# Strict mode for batch validation
vcf-agent samspec batch-validate *.vcf --strict --summary
```

### Explain Violations

```bash
# Explain all violations in a file
vcf-agent samspec explain sample.vcf

# Filter by violation level
vcf-agent samspec explain sample.vcf --level critical
vcf-agent samspec explain sample.vcf --level warning
vcf-agent samspec explain sample.vcf --level info

# Explain specific rule
vcf-agent samspec explain sample.vcf --rule-id MISSING_FILEFORMAT
```

## Python API Usage

### Basic Validation

```python
from vcf_agent.samspec_compliance import validate_vcf_samspec_compliance

# Validate a VCF file
report = validate_vcf_samspec_compliance("sample.vcf")

# Check compliance status
if report.is_compliant:
    print("File is SAMspec compliant!")
else:
    print(f"File has {report.critical_count} critical violations")

# Access violation details
for violation in report.violations:
    print(f"{violation.rule_id}: {violation.message}")
    if violation.suggestion:
        print(f"  Suggestion: {violation.suggestion}")
```

### Advanced Usage

```python
from vcf_agent.samspec_compliance import SAMspecValidator, ComplianceLevel

# Create validator instance
validator = SAMspecValidator()

# Validate file
report = validator.validate_file("sample.vcf")

# Filter violations by level
critical_violations = [
    v for v in report.violations 
    if v.level == ComplianceLevel.CRITICAL
]

# Convert report to dictionary for JSON serialization
report_dict = report.to_dict()

# Access detailed information
print(f"VCF Version: {report.vcf_version}")
print(f"Total violations: {report.total_violations}")
print(f"Critical: {report.critical_count}")
print(f"Warnings: {report.warning_count}")
print(f"Info: {report.info_count}")
```

## Report Structure

### ComplianceReport

The validation returns a `ComplianceReport` object with the following attributes:

- `file_path`: Path to the validated VCF file
- `vcf_version`: Detected VCF version (e.g., "4.2")
- `total_violations`: Total number of violations found
- `critical_count`: Number of critical violations
- `warning_count`: Number of warning violations
- `info_count`: Number of info violations
- `violations`: List of `ComplianceViolation` objects
- `is_compliant`: Boolean indicating if file is compliant (no critical violations)

### ComplianceViolation

Each violation includes:

- `level`: Violation severity (`CRITICAL`, `WARNING`, `INFO`)
- `rule_id`: Unique identifier for the validation rule
- `message`: Human-readable description of the violation
- `line_number`: Line number where violation occurred (if applicable)
- `field`: Field name related to the violation (if applicable)
- `value`: Problematic value (if applicable)
- `suggestion`: Suggested fix for the violation (if applicable)

## Example Reports

### Text Format

```
SAMspec Compliance Report
==================================================
File: sample.vcf
VCF Version: 4.2
Compliant: No

Violation Summary:
  Critical: 2
  Warnings: 1
  Info: 0
  Total: 3

Detailed Violations:
------------------------------

1. MISSING_FILEFORMAT (CRITICAL)
   Missing required ##fileformat header line
   Suggestion: Add ##fileformat=VCFv4.2 or ##fileformat=VCFv4.3

2. INVALID_POS_FORMAT (CRITICAL)
   POS field must be an integer, got 'invalid_pos'
   Line: 15
   Field: POS

3. MISSING_PASS_FILTER (WARNING)
   Missing PASS filter definition
   Suggestion: Add ##FILTER=<ID=PASS,Description="All filters passed">
```

### JSON Format

```json
{
  "file_path": "sample.vcf",
  "vcf_version": "4.2",
  "total_violations": 3,
  "critical_count": 2,
  "warning_count": 1,
  "info_count": 0,
  "is_compliant": false,
  "violations": [
    {
      "level": "critical",
      "rule_id": "MISSING_FILEFORMAT",
      "message": "Missing required ##fileformat header line",
      "line_number": null,
      "field": null,
      "value": null,
      "suggestion": "Add ##fileformat=VCFv4.2 or ##fileformat=VCFv4.3"
    },
    {
      "level": "critical",
      "rule_id": "INVALID_POS_FORMAT",
      "message": "POS field must be an integer, got 'invalid_pos'",
      "line_number": 15,
      "field": "POS",
      "value": "invalid_pos",
      "suggestion": null
    },
    {
      "level": "warning",
      "rule_id": "MISSING_PASS_FILTER",
      "message": "Missing PASS filter definition",
      "line_number": null,
      "field": null,
      "value": null,
      "suggestion": "Add ##FILTER=<ID=PASS,Description=\"All filters passed\">"
    }
  ]
}
```

## Integration with CI/CD

### Exit Codes

The CLI commands use standard exit codes:

- `0`: Validation passed (no critical violations)
- `1`: Validation failed (critical violations found, or warnings in strict mode)

### Example GitHub Actions Workflow

```yaml
name: VCF Compliance Check

on: [push, pull_request]

jobs:
  validate-vcf:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v2
    
    - name: Set up Python
      uses: actions/setup-python@v2
      with:
        python-version: '3.9'
    
    - name: Install VCF Agent
      run: pip install vcf-analysis-agent
    
    - name: Validate VCF files
      run: |
        vcf-agent samspec batch-validate data/*.vcf \
          --strict \
          --summary \
          --output-dir reports/ \
          --format json
    
    - name: Upload reports
      uses: actions/upload-artifact@v2
      with:
        name: compliance-reports
        path: reports/
```

## Performance Considerations

- **Large files**: The validator processes files line-by-line for memory efficiency
- **Compressed files**: Gzip-compressed files are supported with minimal overhead
- **Validation speed**: Typical validation speed is 10,000+ variants per second
- **Memory usage**: Memory usage scales with header complexity, not file size

## Best Practices

1. **Regular validation**: Include SAMspec validation in your VCF processing pipeline
2. **Strict mode**: Use strict mode in production environments to catch all issues
3. **Batch processing**: Use batch validation for multiple files to get summary reports
4. **JSON output**: Use JSON output for programmatic processing and integration
5. **Fix suggestions**: Follow the provided suggestions to resolve violations
6. **Version consistency**: Ensure all VCF files in a project use the same version

## Troubleshooting

### Common Issues

1. **File not found**: Ensure the VCF file path is correct and accessible
2. **Permission errors**: Check file permissions for read access
3. **Encoding issues**: Ensure VCF files use UTF-8 encoding
4. **Large file timeouts**: For very large files, consider splitting or using streaming validation

### Getting Help

For issues with SAMspec compliance validation:

1. Check the violation suggestions for specific fixes
2. Refer to the VCF specification documentation
3. Use the `explain` command for detailed violation information
4. Review the test files in the repository for examples

## Related Features

- **VCF Validation**: Basic VCF file structure validation
- **bcftools Integration**: Use bcftools for additional validation
- **Golden File Testing**: Compare validation results against reference outputs
- **Metrics and Monitoring**: Track compliance metrics over time 