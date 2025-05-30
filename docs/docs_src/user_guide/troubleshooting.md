# Troubleshooting

Common issues and solutions for the CLI Enhanced Validation Engine.

## Quick Diagnostics

### Health Check Command

```bash
# Check system health and configuration
python scripts/cli_enhanced_validation.py --stats
```

**Expected Output**:
```
📊 Performance Statistics:
Cache: enabled, Parallel workers: 4, Git integration: available
Cache hit rate: 87.5%, Cache entries: 45, Cache size: 12.3 MB
```

### Common First Steps

```bash
# 1. Check version and basic functionality
python scripts/cli_enhanced_validation.py --help

# 2. Test with verbose output
python scripts/cli_enhanced_validation.py --verbose src/vcf_agent/cli/

# 3. Clear cache if issues persist
python scripts/cli_enhanced_validation.py --clear-cache

# 4. Run without cache for baseline
python scripts/cli_enhanced_validation.py --no-cache src/vcf_agent/cli/
```

## Common Issues and Solutions

### 1. Import and Dependency Issues

#### Problem: `ModuleNotFoundError: No module named 'docstring_parser'`

```
Error: ModuleNotFoundError: No module named 'docstring_parser'
```

**Solution**:
```bash
# Install missing dependency
pip install docstring_parser

# Or install all optional dependencies
pip install "vcf-analysis-agent[cli-validation]"

# Verify installation
python -c "import docstring_parser; print('✅ docstring_parser available')"
```

#### Problem: `ModuleNotFoundError: No module named 'pydantic_settings'`

```
Warning: pydantic_settings not available, configuration features limited
```

**Solution**:
```bash
# Install pydantic settings (Python 3.8+)
pip install pydantic-settings

# For older Python versions
pip install "pydantic<2.0"

# Alternative: Use environment variables instead of config files
export CLI_VALIDATION_MODE=comprehensive
```

#### Problem: `ImportError: cannot import name 'tomllib'`

```
Warning: TOML support not available, .toml config files will be ignored
```

**Solution**:
```bash
# For Python < 3.11
pip install tomli

# Or use YAML configuration instead
cp config.toml config.yml  # Convert to YAML format
```

### 2. Performance Issues

#### Problem: Validation is Very Slow

**Symptoms**:
```
⏱️  Execution time: 45.32s
💾 Cache hit rate: 12%
```

**Diagnosis**:
```bash
# Check cache status
python scripts/cli_enhanced_validation.py --stats

# Profile with verbose output
python scripts/cli_enhanced_validation.py --verbose --no-cache src/
```

**Solutions**:

=== "Enable Caching"
    ```bash
    # Ensure cache is enabled
    python scripts/cli_enhanced_validation.py --stats
    
    # Clear and rebuild cache
    python scripts/cli_enhanced_validation.py --clear-cache
    python scripts/cli_enhanced_validation.py src/  # Rebuild cache
    ```

=== "Increase Parallelism"
    ```bash
    # Use more workers for large codebases
    python scripts/cli_enhanced_validation.py --parallel-workers 8 src/
    
    # Find optimal worker count
    for workers in 2 4 6 8; do
        echo "Testing $workers workers:"
        time python scripts/cli_enhanced_validation.py --parallel-workers $workers src/
    done
    ```

=== "Use Incremental Validation"
    ```bash
    # Only validate changed files
    python scripts/cli_enhanced_validation.py --incremental
    
    # Use quick mode for development
    python scripts/cli_enhanced_validation.py --mode quick src/
    ```

#### Problem: High Memory Usage

**Symptoms**:
```bash
# Memory usage monitoring
ps aux | grep cli_enhanced_validation
# Shows high RSS memory usage
```

**Solutions**:
```bash
# Reduce parallel workers
python scripts/cli_enhanced_validation.py --parallel-workers 2 src/

# Process files in smaller batches
python scripts/cli_enhanced_validation.py src/vcf_agent/cli/
python scripts/cli_enhanced_validation.py tests/

# Clear cache to reduce memory footprint
python scripts/cli_enhanced_validation.py --clear-cache
```

### 3. Git Integration Issues

#### Problem: `git: command not found` or Git errors

```
Error: Git integration failed: git: command not found
```

**Solutions**:

=== "Install Git"
    ```bash
    # Ubuntu/Debian
    sudo apt-get install git
    
    # macOS
    brew install git
    
    # Or use without git integration
    python scripts/cli_enhanced_validation.py --no-git src/
    ```

=== "Check Git Repository"
    ```bash
    # Ensure you're in a git repository
    git status
    
    # Initialize if needed
    git init
    git add .
    git commit -m "Initial commit"
    ```

#### Problem: Incremental validation not finding changes

```
🔄 Incremental validation against 'main'
📁 Changed files: 0, Files analyzed: 0
```

**Diagnosis**:
```bash
# Check git diff manually
git diff --name-only main

# Check current branch
git branch --show-current

# Check git history
git log --oneline -5
```

**Solutions**:
```bash
# Specify correct base branch
python scripts/cli_enhanced_validation.py --incremental --base-branch develop

# Use different comparison target
python scripts/cli_enhanced_validation.py --incremental --base-branch HEAD~1

# Fall back to full validation
python scripts/cli_enhanced_validation.py src/
```

### 4. Configuration Issues

#### Problem: Configuration file not found or ignored

```
Warning: Configuration file not found, using defaults
```

**Diagnosis**:
```bash
# Check configuration discovery
python scripts/cli_enhanced_validation.py --verbose src/ 2>&1 | grep -i config

# List configuration files in current directory
ls -la .cli-validation.* pyproject.toml
```

**Solutions**:
```bash
# Create basic configuration
cat > .cli-validation.yml << EOF
validation:
  mode: comprehensive
  coverage_threshold: 90.0
cache:
  enabled: true
EOF

# Specify configuration explicitly
python scripts/cli_enhanced_validation.py --config custom-config.yml src/

# Validate configuration syntax
python -c "import yaml; yaml.safe_load(open('.cli-validation.yml'))"
```

#### Problem: Invalid configuration values

```
Error: Invalid configuration: coverage_threshold must be between 0 and 100
```

**Solutions**:
```yaml
# Fix configuration values
validation:
  mode: comprehensive  # Must be: quick, comprehensive, strict
  coverage_threshold: 95.0  # Must be 0.0-100.0
  parallel_workers: 4  # Must be 1-16

cache:
  enabled: true  # Must be boolean
  ttl_hours: 24  # Must be positive integer
```

### 5. Validation Errors and Warnings

#### Problem: Low documentation coverage

```
❌ 712 definitions, 45.2% coverage (below 95% threshold)
```

**Analysis**:
```bash
# Get detailed coverage report
python scripts/cli_enhanced_validation.py --verbose src/ | grep -A 10 "Coverage Analysis"

# Focus on specific modules
python scripts/cli_enhanced_validation.py --verbose src/vcf_agent/cli/
```

**Solutions**:

=== "Lower Threshold Temporarily"
    ```yaml
    # .cli-validation.yml
    validation:
      coverage_threshold: 75.0  # Temporarily lower
    ```

=== "Focus on Critical Files"
    ```bash
    # Validate only CLI-specific files
    python scripts/cli_enhanced_validation.py src/vcf_agent/cli/
    
    # Skip test files
    python scripts/cli_enhanced_validation.py src/ --exclude tests/
    ```

=== "Improve Documentation"
    ```python
    # Add missing docstrings
    def validate_cli_command(command: str) -> bool:
        """
        Validate CLI command format and structure.
        
        Args:
            command: CLI command string to validate
            
        Returns:
            True if command is valid, False otherwise
            
        Raises:
            ValueError: If command format is invalid
        """
        pass
    ```

#### Problem: Structured parsing errors

```
Warning: unknown_docstring_format - Unable to parse docstring format
```

**Solutions**:

=== "Fix Docstring Format"
    ```python
    # Before (problematic)
    def process_data(data):
        """Process data
        data -- input data
        returns processed result
        """
        
    # After (Google format)
    def process_data(data):
        """Process data according to specified rules.
        
        Args:
            data: Input data to process
            
        Returns:
            Processed result data
        """
    ```

=== "Use Different Format"
    ```python
    # NumPy format
    def process_data(data):
        """
        Process data according to specified rules.
        
        Parameters
        ----------
        data : str
            Input data to process
            
        Returns
        -------
        str
            Processed result data
        """
    ```

### 6. CI/CD Integration Issues

#### Problem: GitHub Actions workflow fails

```
Error: Process completed with exit code 1
```

**Diagnosis in CI**:
```yaml
# Add debugging to workflow
- name: Debug CLI Validation
  run: |
    python scripts/cli_enhanced_validation.py --stats
    python scripts/cli_enhanced_validation.py --verbose --mode quick src/vcf_agent/cli/
```

**Common Solutions**:

=== "Missing Dependencies"
    ```yaml
    # Install all dependencies in CI
    - name: Install dependencies
      run: |
        pip install -r requirements.txt
        pip install docstring_parser pydantic-settings
    ```

=== "Git History Issues"
    ```yaml
    # Ensure full git history for incremental validation
    - uses: actions/checkout@v3
      with:
        fetch-depth: 0  # Important for git diff
    ```

=== "Permission Issues"
    ```yaml
    # Fix cache permissions
    - name: Setup cache
      run: |
        mkdir -p ~/.cache/cli-validation
        chmod 755 ~/.cache/cli-validation
    ```

#### Problem: Pre-commit hook failures

```
CLI Documentation Validation............................Failed
- hook id: cli-doc-validation
- exit code: 1
```

**Solutions**:
```bash
# Test hook manually
pre-commit run cli-doc-validation --all-files

# Debug specific files
python scripts/cli_enhanced_validation.py --staged-files --verbose

# Temporarily bypass for urgent fixes
git commit --no-verify -m "Emergency fix"

# Fix hook configuration
# .pre-commit-config.yaml
repos:
  - repo: local
    hooks:
      - id: cli-doc-validation
        name: CLI Documentation Validation
        entry: python scripts/cli_enhanced_validation.py --staged-files --mode quick
        language: system
        types: [python]
        pass_filenames: false  # Important!
```

## Performance Troubleshooting

### Identifying Bottlenecks

```bash
# Profile execution time
time python scripts/cli_enhanced_validation.py --verbose src/

# Check file-by-file performance
python scripts/cli_enhanced_validation.py --verbose src/ 2>&1 | grep "Processing"

# Monitor system resources
# Terminal 1: Run validation
python scripts/cli_enhanced_validation.py src/
# Terminal 2: Monitor resources
top -p $(pgrep -f cli_enhanced_validation)
```

### Cache Troubleshooting

```bash
# Check cache status
python scripts/cli_enhanced_validation.py --stats

# Manual cache inspection
ls -la ~/.cache/cli-validation/

# Cache performance test
echo "First run (cold cache):"
time python scripts/cli_enhanced_validation.py --no-cache src/vcf_agent/cli/
echo "Second run (warm cache):"
time python scripts/cli_enhanced_validation.py src/vcf_agent/cli/

# Cache corruption recovery
python scripts/cli_enhanced_validation.py --clear-cache
rm -rf ~/.cache/cli-validation/
```

## Error Message Reference

### Exit Codes

| Code | Meaning | Troubleshooting |
|------|---------|----------------|
| `0` | Success | No action needed |
| `1` | Validation failures | Check coverage and fix documentation |
| `2` | Script error | Check dependencies and configuration |

### Common Error Patterns

#### Import Errors
```python
# Error pattern
ModuleNotFoundError: No module named 'X'

# Solution template
pip install X
# or
pip install "package[extra]"
```

#### Configuration Errors
```yaml
# Error pattern
Error: Invalid configuration: field 'X' validation failed

# Solution: Check configuration schema
validation:
  mode: "comprehensive"  # String, not bare word
  coverage_threshold: 95.0  # Float, not integer
```

#### Git Errors
```bash
# Error pattern
fatal: not a git repository

# Solution
git init
# or run without git integration
python scripts/cli_enhanced_validation.py --no-git
```

## Advanced Debugging

### Enable Debug Logging

```python
# Add to script or config
import logging
logging.basicConfig(level=logging.DEBUG)

# Or use environment variable
export CLI_VALIDATION_LOG_LEVEL=DEBUG
python scripts/cli_enhanced_validation.py src/
```

### Custom Debugging

```python
# Create debug script
#!/usr/bin/env python3
import sys
sys.path.insert(0, '.')

from scripts.cli_enhanced_validation import EnhancedCLIValidator

# Debug specific file
validator = EnhancedCLIValidator()
result = validator.validate_file('src/vcf_agent/cli/main.py')
print(f"Definitions: {len(result.definitions)}")
print(f"Issues: {len(result.issues)}")
for issue in result.issues:
    print(f"  {issue.severity}: {issue.message}")
```

### Memory Profiling

```bash
# Install memory profiler
pip install memory_profiler

# Profile memory usage
python -m memory_profiler scripts/cli_enhanced_validation.py src/

# Monitor memory during execution
mprof run python scripts/cli_enhanced_validation.py src/
mprof plot
```

## Getting Help

### Information Gathering

Before reporting issues, gather this information:

```bash
# System information
python --version
pip list | grep -E "(docstring|pydantic|yaml)"

# Tool information
python scripts/cli_enhanced_validation.py --stats

# Error reproduction
python scripts/cli_enhanced_validation.py --verbose src/ > debug.log 2>&1
```

### Useful Commands for Support

```bash
# Minimal reproduction case
python scripts/cli_enhanced_validation.py --mode quick --verbose single_file.py

# Configuration debugging
python scripts/cli_enhanced_validation.py --config /dev/null --verbose src/

# Clean environment test
python scripts/cli_enhanced_validation.py --no-cache --parallel-workers 1 src/
```

### Reporting Issues

Include this information when reporting issues:

1. **Command used**: Full command line
2. **Error output**: Complete error message
3. **Environment**: Python version, OS, dependencies
4. **Configuration**: Content of config files
5. **Sample file**: Minimal example that reproduces the issue

**Example issue report**:
```
Command: python scripts/cli_enhanced_validation.py --mode comprehensive src/
Error: ModuleNotFoundError: No module named 'tomli'
Python: 3.9.7
OS: Ubuntu 20.04
Config: Default (.cli-validation.yml not present)
Sample: Any Python file with docstrings
```

This troubleshooting guide covers the most common issues users encounter, providing both quick fixes and detailed diagnostic procedures for complex problems.
