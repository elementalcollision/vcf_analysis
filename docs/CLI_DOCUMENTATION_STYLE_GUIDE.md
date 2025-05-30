# VCF Analysis Agent CLI Documentation Style Guide

## Overview

This document provides comprehensive standards, templates, and validation tools for documenting CLI command handlers in the VCF Analysis Agent project. Following these standards ensures consistent, high-quality documentation that serves both developers and end users.

## Implementation Status

✅ **COMPLETED**: CLI Documentation Style Guide (CLI-DOC-001)
- Comprehensive documentation standards established
- Multiple command-specific templates created
- Validation framework with automated checking
- Real-world examples for common CLI patterns
- Integration with argparse help system

## Related Tools

This style guide works alongside other CLI validation tools in the `scripts/` directory:

### `scripts/validate_cli_documentation.py`
- **Purpose**: Validates CLI module-level documentation completeness
- **Scope**: Ensures all implemented commands are documented in the module docstring
- **Usage**: `python scripts/validate_cli_documentation.py --verbose`
- **Focus**: Command presence/absence, example accuracy, overall CLI documentation

### `scripts/cli_documentation_style_guide.py` (This Tool)
- **Purpose**: Provides standards and validates individual function docstrings
- **Scope**: Ensures each command handler function has comprehensive documentation
- **Usage**: `python scripts/cli_documentation_style_guide.py`
- **Focus**: Docstring quality, style compliance, template guidance

### Complementary Workflow
1. **Development**: Use style guide templates for new command handlers
2. **Individual Validation**: Use style guide validator for function-level documentation
3. **Module Validation**: Use `validate_cli_documentation.py` for overall CLI completeness
4. **CI/CD**: Run both validators to ensure comprehensive documentation coverage

## Core Principles

### 1. Dual-Purpose Documentation
Every CLI command handler docstring must serve two purposes:
- **Developer Documentation**: Comprehensive technical details for code maintenance
- **User Documentation**: Clear usage instructions accessible via `--help`

### 2. Standardized Structure
All CLI command docstrings must include these required sections:
- **Summary**: One-line imperative description ending with period
- **Description**: Detailed explanation of purpose and context
- **Arguments**: All required parameters with types and constraints
- **Options**: All optional parameters with defaults and behavior
- **Returns**: Exit code information
- **Raises**: All possible exceptions with conditions
- **Examples**: Real CLI usage examples with shell syntax
- **Exit Codes**: Numerical codes and their meanings

### 3. Quality Standards
- **95% Coverage Minimum**: All functions must achieve 95% documentation coverage score
- **PEP 257 Compliance**: Follow Python docstring conventions
- **CLI-Specific Format**: Include proper command invocation examples
- **Imperative Mood**: Use "Execute...", "Process...", not "Executes..."

## Available Templates

### General CLI Commands
```python
from scripts.cli_documentation_style_guide import cli_command_template
template = cli_command_template()
```

### LanceDB Operations
```python
from scripts.cli_documentation_style_guide import lancedb_command_template
template = lancedb_command_template()
```

### VCF Processing
```python
from scripts.cli_documentation_style_guide import vcf_processing_template
template = vcf_processing_template()
```

### Compliance Validation
```python
from scripts.cli_documentation_style_guide import compliance_command_template
template = compliance_command_template()
```

## Template Customization

Use the template generation utility for specific commands:

```python
from scripts.cli_documentation_style_guide import generate_command_docstring

# Generate a custom LanceDB command docstring
docstring = generate_command_docstring(
    'lancedb',
    command_purpose='initialize variant database',
    database_operation='table creation and schema setup',
    additional_args='--schema-file (str): Path to custom schema definition'
)
```

## Validation

### Automated Validation
```python
from scripts.cli_documentation_style_guide import CLIDocstringValidator

validator = CLIDocstringValidator()
result = validator.validate_docstring(your_cli_function)

if not result.is_valid:
    print(f"Coverage Score: {result.coverage_score:.1f}%")
    for issue in result.issues:
        print(f"Issue: {issue}")
    for suggestion in result.suggestions:
        print(f"Suggestion: {suggestion}")
```

### Coverage Requirements
- **Minimum Score**: 95%
- **Required Sections**: All 8 core sections must be present
- **Format Compliance**: PEP 257 standards
- **CLI Integration**: Examples must show proper invocation format

## Real-World Examples

### Example 1: LanceDB Initialization
```python
def init_lancedb_command():
    """Initialize LanceDB table for variant storage.

    Creates a new LanceDB database and table structure optimized for storing
    genomic variant data with vector embeddings. The table schema includes
    variant identifiers, genomic coordinates, reference/alternate alleles,
    clinical significance annotations, and high-dimensional embedding vectors
    for similarity searches.

    Arguments:
        --db-path (str): Path to LanceDB database directory.
            Directory will be created if it doesn't exist. Must have
            write permissions for database creation and management.
        --table-name (str, optional): Name for the variants table.
            Default: 'variants'. Must be valid table name.

    Options:
        --schema-version (str, optional): Table schema version to use.
            Default: 'latest'. Specify version for compatibility.
        --embedding-dim (int, optional): Dimension of embedding vectors.
            Default: 384. Must match embedding model output dimensions.

    Returns:
        int: Exit code (0 for success, non-zero for errors).

    Raises:
        FileNotFoundError: When parent directory for db-path doesn't exist.
        PermissionError: When lacking write permissions for database creation.
        ValueError: When table-name contains invalid characters.
        lancedb.Error: When database initialization fails.

    Examples:
        Initialize with default settings:
            $ python -m vcf_agent.cli init-lancedb --db-path ./variants_db
        
        Custom table name:
            $ python -m vcf_agent.cli init-lancedb --db-path ./my_db \\
                --table-name my_variants
        
        Specify embedding dimensions:
            $ python -m vcf_agent.cli init-lancedb --db-path ./db \\
                --embedding-dim 512

    Exit Codes:
        0: Database initialized successfully
        1: General initialization error
        2: Invalid arguments provided
        3: Parent directory not found
        4: Permission denied for database creation
        5: Invalid table name format
        6: Database system error
        7: Insufficient disk space

    Performance:
        - Initial database creation is fast (< 1 second)
        - Directory creation includes metadata and schema files
        - No significant memory or CPU requirements for initialization

    See Also:
        - add-variant: Add individual variants to initialized database
        - ingest-vcf: Bulk load variants from VCF files
        - create-lancedb-index: Optimize database performance
    """
```

### Example 2: VCF Processing
```python
def ingest_vcf_command():
    """Ingest VCF file into both LanceDB and Kuzu databases.

    Processes a VCF file through the complete ingestion pipeline, extracting
    variant information, generating embeddings, and storing data in both the
    LanceDB vector database (for similarity searches) and Kuzu graph database
    (for relationship analysis).

    Arguments:
        --vcf-file (str): Path to input VCF file (.vcf, .vcf.gz, .vcf.bgz).
            File must be valid VCF 4.0+ format and accessible for reading.

    Options:
        --lancedb-path (str, optional): LanceDB database directory.
            Default: './lancedb'. Database will be created if it doesn't exist.
        --kuzu-path (str, optional): Kuzu database directory.
            Default: './kuzu_db'. Database schema will be initialized if needed.
        --batch-size (int, optional): Variants to process per batch.
            Default: 1000. Larger values use more memory but may be faster.

    Returns:
        int: Exit code (0 for success, non-zero for errors).

    Raises:
        FileNotFoundError: When VCF file or reference files cannot be found.
        ValueError: When VCF format validation fails or arguments are invalid.
        PermissionError: When database directories cannot be accessed.

    Examples:
        Basic ingestion:
            $ python -m vcf_agent.cli ingest-vcf --vcf-file sample.vcf.gz
        
        Custom database paths:
            $ python -m vcf_agent.cli ingest-vcf --vcf-file data.vcf \\
                --lancedb-path ./vector_db --kuzu-path ./graph_db

    Exit Codes:
        0: Ingestion completed successfully
        1: General processing error
        2: Invalid command arguments
        3: VCF file not found
        4: Invalid VCF format
        5: Database access error

    See Also:
        - init-lancedb: Initialize vector database
        - populate-kuzu-from-vcf: Graph database only ingestion
        - samspec validate: VCF format validation
    """
```

## Integration with Development Workflow

### Pre-commit Hooks
The style guide integrates with pre-commit hooks for immediate validation:
```yaml
# .pre-commit-config.yaml
repos:
  - repo: local
    hooks:
      - id: cli-docstring-validation
        name: CLI Docstring Validation
        entry: python scripts/cli_documentation_style_guide.py
        language: python
        files: ^src/vcf_agent/cli.*\.py$
```

### IDE Integration
The validation framework works with popular IDEs:
- **VS Code**: Automatic docstring validation on save
- **PyCharm**: Integration with code quality tools
- **Vim/Neovim**: Linting integration with ale or coc

### CI/CD Pipeline
Automated validation in continuous integration:
```yaml
# .github/workflows/documentation.yml
- name: Validate CLI Documentation
  run: |
    python -c "
    import sys
    sys.path.append('scripts')
    from cli_documentation_style_guide import validate_module_documentation
    results = validate_module_documentation('src/vcf_agent/cli.py')
    failing = [name for name, result in results.items() if not result.is_valid]
    if failing:
        print(f'Documentation validation failed for: {failing}')
        exit(1)
    "
```

## Common Patterns and Best Practices

### 1. Command Naming
- Use descriptive, action-oriented names
- Follow kebab-case convention: `init-lancedb`, `ingest-vcf`
- Group related commands: `lancedb-*`, `vcf-*`, `compliance-*`

### 2. Argument Documentation
- Always specify type and constraints
- Include examples of valid values
- Document default behavior clearly
- Note file format requirements

### 3. Error Handling
- Document all possible exceptions
- Map exceptions to specific exit codes
- Provide remediation guidance
- Include troubleshooting notes

### 4. Examples Section
- Show common usage patterns first
- Include complex scenarios with line continuations
- Use realistic file paths and parameters
- Demonstrate error recovery patterns

### 5. Performance Notes
- Document memory requirements
- Note processing time expectations
- Include scaling characteristics
- Suggest optimization strategies

## Migration Guide

### For Existing Commands
1. **Audit Current Documentation**: Use the validation framework to identify gaps
2. **Apply Templates**: Choose appropriate template for command type
3. **Customize Content**: Fill in command-specific details
4. **Validate**: Ensure 95% coverage score
5. **Test Integration**: Verify argparse help output

### For New Commands
1. **Start with Template**: Choose appropriate command template
2. **Customize During Development**: Fill sections as implementation progresses
3. **Validate Early**: Check coverage score regularly
4. **Review Examples**: Ensure examples work as documented

## Quality Gates

### Development Phase
- [ ] Docstring template selected and customized
- [ ] All required sections present
- [ ] Examples tested and verified
- [ ] Coverage score ≥ 95%

### Code Review Phase
- [ ] Reviewer validates documentation completeness
- [ ] Examples tested in realistic scenarios
- [ ] Integration with help system verified
- [ ] Performance notes reviewed for accuracy

### Pre-Release Phase
- [ ] Full documentation validation suite passes
- [ ] User testing of help output completed
- [ ] Documentation renders correctly in all contexts
- [ ] No validation errors in CI/CD pipeline

## Troubleshooting

### Common Issues

**Low Coverage Score**
- Missing required sections (most common)
- Insufficient detail in existing sections
- No CLI-specific examples provided
- Missing exit code documentation

**Format Validation Errors**
- Summary doesn't end with period
- Lines exceed 72 character limit
- Missing blank line after summary
- Non-imperative mood in summary

**Integration Issues**
- Examples don't match actual CLI interface
- Documented arguments don't exist in implementation
- Exit codes don't match actual command behavior
- Help text truncated or malformed

### Getting Help

1. **Validation Framework**: Use built-in validator for detailed feedback
2. **Template Examples**: Reference provided real-world examples
3. **Team Review**: Request documentation review during code review
4. **Style Guide Reference**: Consult this document for clarification

## Future Enhancements

The style guide will evolve to include:
- Additional command type templates
- Enhanced validation rules
- Integration with documentation generators
- Automated example testing
- Performance benchmark integration

## Related Documentation

- [CLI Module Architecture](./CLI_ARCHITECTURE.md)
- [Development Workflow](./DEVELOPMENT_WORKFLOW.md)
- [Code Review Guidelines](./CODE_REVIEW_GUIDELINES.md)
- [Testing Standards](./TESTING_STANDARDS.md)

---

**Version**: 1.0  
**Last Updated**: 2025-01-05  
**Status**: ✅ Implemented and Active 