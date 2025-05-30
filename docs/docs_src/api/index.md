# API Reference

Complete API documentation for the CLI Enhanced Validation Engine, automatically generated from source code docstrings.

## Quick Navigation

=== "Core Components"
    - **[CLI Enhanced Validation Engine](../../api/index.html){target=_blank}** - Main validation engine and tools
    - **[Module Documentation](../../api/_modules/cli_enhanced_validation.html){target=_blank}** - Source code with syntax highlighting
    - **[Python Module Index](../../api/py-modindex.html){target=_blank}** - Complete module index
    - **[General Index](../../api/genindex.html){target=_blank}** - Alphabetical index of all items

=== "Configuration"
    - **[Validation Configuration Guide](../user_guide/configuration.md)** - Configuration management documentation
    - **[CLI Reference](../user_guide/cli_reference.md)** - Command-line interface reference
    - **[Examples](../user_guide/examples.md)** - Usage examples and patterns

=== "Integration"
    - **[GitHub Actions](../integration/github_actions.md)** - CI/CD integration
    - **[Pre-commit Hooks](../integration/pre_commit_hooks.md)** - Development workflow integration
    - **[IDE Integration](../integration/ide_integration.md)** - Editor and IDE setup

## Complete API Documentation

For complete API documentation with all classes, methods, and detailed docstrings, visit:

**[📚 Complete CLI Enhanced Validation Engine API Reference](../../api/index.html){target=_blank}**

This comprehensive reference includes:

- **All Classes and Methods**: Complete coverage of the validation engine
- **Detailed Docstrings**: Parameter descriptions, return values, and examples
- **Type Annotations**: Full type information for all functions and methods
- **Source Code Links**: Direct links to implementation code
- **Cross-references**: Links between related components

## Architecture Overview

```mermaid
graph TB
    CLI[CLI Interface] --> Validator[EnhancedCLIValidator]
    Validator --> AST[ASTAnalyzer]
    Validator --> Parser[MultiFormatDocstringParser]
    Validator --> Cache[CacheManager]
    Validator --> Git[GitIntegration]
    
    AST --> Definitions[CodeDefinition]
    Parser --> Parsed[ParsedDocstring]
    Cache --> Performance[PerformanceOptimizer]
    
    Validator --> Report[ValidationReport]
    Report --> Output[OutputFormat]
    
    subgraph "Core Analysis"
        AST
        Parser
        Definitions
        Parsed
    end
    
    subgraph "Performance"
        Cache
        Performance
        Git
    end
    
    subgraph "Results"
        Report
        Output
    end
```

## Key Components

### EnhancedCLIValidator

The main validation engine that orchestrates all validation activities.

**Key Methods**:
- `validate_comprehensive()` - Complete validation with full analysis
- `validate_incremental()` - Git-based incremental validation
- `validate_staged_files()` - Pre-commit validation
- `generate_ci_report()` - CI/CD integration reporting

### ASTAnalyzer

Advanced AST-based code analysis for extracting function and class definitions.

**Key Methods**:
- `extract_all_definitions()` - Extract all code definitions from files
- `analyze_docstring_coverage()` - Calculate documentation coverage
- `extract_function_signature()` - Extract detailed function signatures

### MultiFormatDocstringParser

Intelligent docstring parsing supporting multiple formats (Google, NumPy, Sphinx).

**Key Methods**:
- `parse_docstring()` - Parse docstring into structured format
- `validate_structure()` - Validate docstring completeness and accuracy

### CacheManager

High-performance caching system for validation results.

**Key Methods**:
- `get()` / `set()` - Cache operations
- `get_stats()` - Cache performance metrics
- `invalidate()` - Cache invalidation

## Usage Examples

### Basic Validation

```python
from cli_enhanced_validation import EnhancedCLIValidator

# Initialize validator
validator = EnhancedCLIValidator()

# Validate files
result = validator.validate_comprehensive(['src/'])

# Check results
print(f"Coverage: {result.coverage_percentage}%")
print(f"Issues: {len(result.issues)}")
```

### Configuration

```python
from cli_enhanced_validation import ValidationConfig, ValidationMode

# Custom configuration
config = ValidationConfig(
    mode=ValidationMode.STRICT,
    coverage_threshold=95.0,
    parallel_workers=4
)

validator = EnhancedCLIValidator(config=config)
```

### Incremental Validation

```python
# Git-based incremental validation
result = validator.validate_incremental(
    base_branch='main',
    include_staged=True
)
```

## Integration with User Guides

This API reference complements the user-focused documentation:

- **[Quick Start Guide](../user_guide/quick_start.md)** - Get started with the CLI tool
- **[Configuration Guide](../user_guide/configuration.md)** - Configure validation behavior
- **[Examples](../user_guide/examples.md)** - Real-world usage examples
- **[Troubleshooting](../user_guide/troubleshooting.md)** - Common issues and solutions

## Development and Extension

For developers looking to extend or integrate the CLI Enhanced Validation Engine:

1. **Study the Core Classes**: Start with `EnhancedCLIValidator` and `ASTAnalyzer`
2. **Understand Data Models**: Review `CodeDefinition` and `ValidationReport`
3. **Explore Configuration**: Use `ValidationConfig` for customization
4. **Performance Optimization**: Leverage `CacheManager` and `PerformanceOptimizer`

## Support and Contributing

- **Issues**: Report bugs and feature requests on GitHub
- **Documentation**: Contribute to improving API documentation
- **Code**: Submit pull requests for enhancements

The API documentation is automatically generated from source code docstrings, ensuring it stays current with the implementation. 