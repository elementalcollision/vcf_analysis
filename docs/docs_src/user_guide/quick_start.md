# Quick Start

Get up and running with the CLI Enhanced Validation Engine in under 5 minutes! This guide will take you from installation to your first successful validation.

## Prerequisites

Before getting started, ensure you have:

- **Python 3.9 or higher** - Check with `python --version`
- **Git** (optional, for source installation and Git integration features)
- **Terminal/Command Prompt** access

## Installation

Choose your preferred installation method:

=== "PyPI (Recommended)"

    ```bash
    # Install the latest stable version
    pip install vcf-cli-validation
    
    # Verify installation
    vcf-validate --version
    ```

=== "From Source"

    ```bash
    # Clone the repository
    git clone https://github.com/vcf-agent/cli-enhanced-validation.git
    cd cli-enhanced-validation
    
    # Install in development mode
    pip install -e .
    
    # Verify installation
    python scripts/cli_enhanced_validation.py --version
    ```

=== "Virtual Environment"

    ```bash
    # Create and activate virtual environment
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate
    
    # Install the CLI tool
    pip install vcf-cli-validation
    
    # Verify installation
    vcf-validate --version
    ```

!!! success "Installation Complete"
    If you see version information, you're ready to proceed!

## Your First Validation

Let's run your first validation to see the CLI Enhanced Validation Engine in action.

### Self-Validation

The best way to understand the tool is to validate itself:

```bash
# Validate the CLI Enhanced Validation Engine
python scripts/cli_enhanced_validation.py scripts/cli_enhanced_validation.py
```

**Expected Output:**
```
🔍 Analyzing CLI Enhanced Validation Engine...
✅ SUCCESS - 121 definitions, 88.4% coverage, 0.44s execution
📊 Performance: Cache hit rate: 0.0%, Parallel workers: 4, Git integration: available

Summary:
- Total definitions: 121
- Documented: 107
- Missing documentation: 14
- Coverage: 88.4%
```

### Validate a Python File

Now let's validate a Python file in your project:

```bash
# Validate a single Python file
python scripts/cli_enhanced_validation.py src/your_module.py

# Example with a sample file
python scripts/cli_enhanced_validation.py src/vcf_agent/cli/__init__.py
```

### Validate a Directory

For comprehensive validation of a Python package:

```bash
# Validate entire directory
python scripts/cli_enhanced_validation.py src/

# Validate with quick mode (faster for large projects)
python scripts/cli_enhanced_validation.py src/ --mode quick
```

## Understanding the Output

The CLI Enhanced Validation Engine provides rich, informative output:

### Success Output
```
✅ SUCCESS - 25 definitions, 96.0% coverage, 0.15s execution
📊 Performance: Cache hit rate: 75.2%, Parallel workers: 4
```

- **Definitions**: Total functions, methods, and classes found
- **Coverage**: Percentage of definitions with docstrings
- **Execution Time**: How long the validation took
- **Cache Hit Rate**: Efficiency of the caching system

### Detailed Analysis
```
📁 Analyzing 5 files...
🔍 Found 25 definitions
📊 Coverage: 96.0% (24/25 definitions documented)
ℹ️  12 info messages
⏱️  Execution time: 0.15s

Files analyzed:
  ✅ src/my_module.py (100% coverage)
  ⚠️  src/utils.py (80% coverage - 2 missing docstrings)
```

### Warning Examples
```
⚠️  Function 'process_data' missing docstring at line 45
ℹ️  Function 'helper_func' has minimal docstring - consider adding parameters and returns
```

## Basic Usage Patterns

### Common Commands

```bash
# Quick validation (recommended for CI/CD)
python scripts/cli_enhanced_validation.py src/ --mode quick

# Comprehensive validation with detailed analysis
python scripts/cli_enhanced_validation.py src/ --mode comprehensive

# JSON output for automation
python scripts/cli_enhanced_validation.py src/ --json

# Verbose output for debugging
python scripts/cli_enhanced_validation.py src/ --verbose

# Help information
python scripts/cli_enhanced_validation.py --help
```

### Validation Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| `quick` | Fast validation, basic checks | CI/CD pipelines, pre-commit hooks |
| `comprehensive` | Full analysis with detailed reporting | Code reviews, quality gates |
| `strict` | Strictest validation rules | Release preparation, compliance |

### Exit Codes

The CLI tool returns standard exit codes for automation:

- **0**: All validations passed ✅
- **1**: Validation issues found ⚠️
- **2**: Script execution error ❌

## Configuration Basics

Create a simple configuration file to customize behavior:

=== "YAML Configuration"

    Create `.cli-validation.yml` in your project root:

    ```yaml
    validation:
      coverage_threshold: 90
      modes: [google, numpy, sphinx]
    
    cache:
      enabled: true
      ttl_hours: 24
    
    performance:
      parallel_workers: 4
    ```

=== "TOML Configuration"

    Create `.cli-validation.toml` in your project root:

    ```toml
    [validation]
    coverage_threshold = 90
    modes = ["google", "numpy", "sphinx"]
    
    [cache]
    enabled = true
    ttl_hours = 24
    
    [performance]
    parallel_workers = 4
    ```

## Performance Features

### Caching

The CLI tool automatically caches validation results:

```bash
# First run (no cache)
python scripts/cli_enhanced_validation.py src/
# ⏱️ Execution time: 2.1s, Cache hit rate: 0.0%

# Second run (with cache)
python scripts/cli_enhanced_validation.py src/
# ⏱️ Execution time: 0.8s, Cache hit rate: 85.2%
```

### Incremental Validation

For large projects, use Git integration for faster validation:

```bash
# Only validate changed files since main branch
python scripts/cli_enhanced_validation.py --incremental --base-branch main

# Only validate staged files
python scripts/cli_enhanced_validation.py --staged-files
```

## Next Steps

Congratulations! You've successfully run your first validations. Here's what to explore next:

### Immediate Next Steps
1. **[Configuration](configuration.md)** - Customize validation rules and behavior
2. **[Examples](examples.md)** - See real-world usage patterns and workflows
3. **[CLI Reference](cli_reference.md)** - Master all available command options

### Integration
4. **[GitHub Actions](../integration/github_actions.md)** - Set up automated validation in CI/CD
5. **[Pre-commit Hooks](../integration/pre_commit_hooks.md)** - Add validation to your development workflow
6. **[IDE Integration](../integration/ide_integration.md)** - Enable real-time validation in your editor

### Advanced Usage
7. **[Performance Tuning](../advanced/performance_tuning.md)** - Optimize for large codebases
8. **[Extending](../advanced/extending.md)** - Add custom validation rules
9. **[Contributing](../advanced/contributing.md)** - Help improve the tool

## Troubleshooting

### Common Issues

!!! warning "Python Version"
    **Error**: `ModuleNotFoundError` or syntax errors
    
    **Solution**: Ensure you're using Python 3.9 or higher:
    ```bash
    python --version  # Should show 3.9+
    ```

!!! warning "Installation Issues"
    **Error**: `pip install` fails
    
    **Solutions**:
    ```bash
    # Update pip first
    pip install --upgrade pip
    
    # Install with user flag if permissions issue
    pip install --user vcf-cli-validation
    ```

!!! warning "No Files Found"
    **Error**: `No Python files found for validation`
    
    **Solution**: Check your file paths and ensure you're in the correct directory:
    ```bash
    # Check current directory
    pwd
    
    # List Python files
    find . -name "*.py" | head -10
    ```

### Getting Help

If you encounter issues:

1. **Check [Troubleshooting](troubleshooting.md)** for detailed solutions
2. **Run with verbose output** for more information:
   ```bash
   python scripts/cli_enhanced_validation.py src/ --verbose
   ```
3. **Check the logs** for detailed error information

---

**Ready to dive deeper?** Check out our [Configuration Guide](configuration.md) to customize the validation engine for your specific needs! 