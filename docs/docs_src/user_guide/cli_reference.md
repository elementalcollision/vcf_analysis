# CLI Reference

Complete command-line interface reference for the CLI Enhanced Validation Engine.

## Basic Usage

```bash
python scripts/cli_enhanced_validation.py [OPTIONS] [PATHS...]
```

The CLI Enhanced Validation Engine analyzes Python files for CLI documentation quality, providing comprehensive validation with caching, Git integration, and multiple output formats.

## Arguments

### Positional Arguments

#### `PATHS`
File or directory paths to validate.

- **Type**: Multiple paths
- **Default**: Current directory (`.`)
- **Examples**:
```bash
# Validate current directory
python scripts/cli_enhanced_validation.py

# Validate specific file
python scripts/cli_enhanced_validation.py src/vcf_agent/cli/main.py

# Validate multiple paths
python scripts/cli_enhanced_validation.py src/vcf_agent/cli/ tests/unit/
```

## Validation Mode Options

### `--mode {quick,comprehensive,strict}`
Validation thoroughness level.

- **Default**: `comprehensive`
- **Choices**:
  - **`quick`**: Fast validation with basic checks
  - **`comprehensive`**: Complete validation with all features
  - **`strict`**: Comprehensive validation with stricter requirements

=== "Quick Mode"
    ```bash
    # Fast validation for development workflow
    python scripts/cli_enhanced_validation.py --mode quick src/vcf_agent/cli/
    ```
    
    **Output**: Basic docstring coverage and structural validation only.

=== "Comprehensive Mode"
    ```bash
    # Complete analysis (default)
    python scripts/cli_enhanced_validation.py --mode comprehensive src/vcf_agent/
    ```
    
    **Output**: Full AST analysis, structured docstring parsing, and performance metrics.

=== "Strict Mode"
    ```bash
    # Strictest validation for production
    python scripts/cli_enhanced_validation.py --mode strict src/vcf_agent/cli/
    ```
    
    **Output**: All comprehensive features plus stricter coverage requirements.

## Configuration Options

### `--config PATH`
Path to configuration file.

- **Default**: Auto-discovery (`.cli-validation.yml`, `.cli-validation.yaml`, `pyproject.toml`)
- **Format**: YAML or TOML

```bash
# Use custom configuration
python scripts/cli_enhanced_validation.py --config custom-validation.yml src/

# Configuration discovery order:
# 1. Specified with --config
# 2. .cli-validation.yml
# 3. .cli-validation.yaml  
# 4. pyproject.toml [tool.cli-validation] section
```

### `--no-cache`
Disable caching for this run.

```bash
# Force fresh analysis without cache
python scripts/cli_enhanced_validation.py --no-cache src/vcf_agent/cli/
```

**Use Cases**:
- Development testing
- CI/CD clean builds
- Debugging cache issues

## Git Integration Options

### `--incremental`
Validate only files changed since base branch.

```bash
# Validate changes against main branch
python scripts/cli_enhanced_validation.py --incremental

# Validate with custom base branch
python scripts/cli_enhanced_validation.py --incremental --base-branch develop
```

**Performance**: Typically 70-90% faster than full validation.

### `--base-branch BRANCH`
Base branch for incremental validation.

- **Default**: `main`
- **Requires**: `--incremental` flag

```bash
# Compare against develop branch
python scripts/cli_enhanced_validation.py --incremental --base-branch develop

# Compare against specific commit
python scripts/cli_enhanced_validation.py --incremental --base-branch HEAD~3
```

### `--staged-files`
Validate only staged files (pre-commit integration).

```bash
# Pre-commit hook usage
python scripts/cli_enhanced_validation.py --staged-files --mode quick

# Combined with JSON output for automation
python scripts/cli_enhanced_validation.py --staged-files --json
```

**Integration**: Perfect for Git pre-commit hooks and automated workflows.

## Output Options

### `--json`
Output results in JSON format for CI/CD integration.

```bash
# JSON output for automation
python scripts/cli_enhanced_validation.py --json src/vcf_agent/cli/

# Combine with other options
python scripts/cli_enhanced_validation.py --incremental --json --github-annotations
```

**JSON Structure**:
```json
{
  "validation_results": {
    "success": true,
    "total_files": 6,
    "total_definitions": 121,
    "coverage_percentage": 88.4
  },
  "issues": [],
  "metadata": {
    "execution_time": 0.44,
    "cache_stats": {"hit_rate_percent": 100.0}
  },
  "github_annotations": []
}
```

### `--verbose`, `-v`
Enable detailed output with debug information.

```bash
# Verbose console output
python scripts/cli_enhanced_validation.py -v src/vcf_agent/cli/

# Verbose JSON output  
python scripts/cli_enhanced_validation.py --json --verbose src/vcf_agent/
```

**Additional Information**:
- File-by-file analysis details
- Performance metrics per file
- Cache hit/miss statistics
- Detailed issue descriptions

### `--github-annotations`
Include GitHub Actions annotations in JSON output.

```bash
# CI/CD integration with inline comments
python scripts/cli_enhanced_validation.py --json --github-annotations src/
```

**GitHub Actions Integration**:
```yaml
- name: Validate CLI Documentation
  run: |
    python scripts/cli_enhanced_validation.py --json --github-annotations . > validation.json
    # Annotations appear as inline PR comments
```

## Utility Options

### `--clear-cache`
Clear validation cache and exit.

```bash
# Clear cache
python scripts/cli_enhanced_validation.py --clear-cache

# Clear cache with JSON output
python scripts/cli_enhanced_validation.py --clear-cache --json
```

**Output**:
```
✅ Cache cleared: 24 entries removed
Cache directory: /Users/dave/.cache/cli-validation
```

### `--stats`
Show performance statistics and exit.

```bash
# Performance overview
python scripts/cli_enhanced_validation.py --stats

# JSON format for monitoring
python scripts/cli_enhanced_validation.py --stats --json
```

**Statistics Include**:
- Cache configuration and hit rates
- Git integration availability
- Performance settings
- Cache size and entries

## Performance Options

### `--parallel-workers N`
Number of parallel workers for validation.

- **Default**: Auto-detected (CPU cores)
- **Range**: 1-16

```bash
# Single-threaded processing
python scripts/cli_enhanced_validation.py --parallel-workers 1 src/

# Maximum parallelism
python scripts/cli_enhanced_validation.py --parallel-workers 8 src/
```

**Performance Impact**: 2-4x speedup on multi-file validation.

### `--timeout N`
Timeout in seconds for validation operations.

- **Default**: 300 seconds
- **Minimum**: 10 seconds

```bash
# Short timeout for CI/CD
python scripts/cli_enhanced_validation.py --timeout 60 src/

# Extended timeout for large codebases
python scripts/cli_enhanced_validation.py --timeout 600 src/
```

## Exit Codes

| Code | Meaning | Description |
|------|---------|-------------|
| `0` | Success | All validation passed |
| `1` | Validation Failures | Issues found below threshold |
| `2` | Execution Error | Script failed to run |

```bash
# Check exit code in scripts
python scripts/cli_enhanced_validation.py src/vcf_agent/cli/
echo "Exit code: $?"

# CI/CD usage
if python scripts/cli_enhanced_validation.py --json src/; then
    echo "✅ Validation passed"
else
    echo "❌ Validation failed with code $?"
fi
```

## Real-World Examples

### Development Workflow

```bash
# Quick check during development
python scripts/cli_enhanced_validation.py --mode quick src/vcf_agent/cli/

# Pre-commit validation
python scripts/cli_enhanced_validation.py --staged-files --mode quick
```

### CI/CD Integration

```bash
# GitHub Actions validation
python scripts/cli_enhanced_validation.py \
  --incremental \
  --json \
  --github-annotations \
  --mode comprehensive \
  src/

# Performance monitoring
python scripts/cli_enhanced_validation.py \
  --stats \
  --json > validation-metrics.json
```

### Enterprise Usage

```bash
# Complete codebase analysis
python scripts/cli_enhanced_validation.py \
  --mode strict \
  --config enterprise-validation.yml \
  --parallel-workers 4 \
  --verbose \
  src/

# Incremental PR validation
python scripts/cli_enhanced_validation.py \
  --incremental \
  --base-branch main \
  --json \
  --github-annotations \
  --timeout 300
```

### Performance Optimization

```bash
# Cache management
python scripts/cli_enhanced_validation.py --stats
python scripts/cli_enhanced_validation.py --clear-cache

# Incremental validation workflow
python scripts/cli_enhanced_validation.py --incremental --mode quick
python scripts/cli_enhanced_validation.py --incremental --mode comprehensive
```

## Common Patterns

### Pre-commit Hook Setup

```bash
# .pre-commit-config.yaml
repos:
  - repo: local
    hooks:
      - id: cli-doc-validation
        name: CLI Documentation Validation
        entry: python scripts/cli_enhanced_validation.py --staged-files --mode quick
        language: system
        types: [python]
```

### Makefile Integration

```makefile
# Makefile
.PHONY: validate-docs validate-docs-quick validate-docs-ci

validate-docs:
	python scripts/cli_enhanced_validation.py --mode comprehensive src/

validate-docs-quick:
	python scripts/cli_enhanced_validation.py --mode quick src/vcf_agent/cli/

validate-docs-ci:
	python scripts/cli_enhanced_validation.py --incremental --json --github-annotations
```

### Configuration File Examples

=== "YAML Configuration"
    ```yaml
    # .cli-validation.yml
    validation:
      mode: comprehensive
      coverage_threshold: 95.0
      parallel_workers: 4
      
    cache:
      enabled: true
      ttl_hours: 24
      
    output:
      json: false
      verbose: false
      github_annotations: false
    ```

=== "TOML Configuration"
    ```toml
    # pyproject.toml
    [tool.cli-validation]
    mode = "comprehensive"
    coverage_threshold = 95.0
    parallel_workers = 4
    
    [tool.cli-validation.cache]
    enabled = true
    ttl_hours = 24
    ```

## Troubleshooting Commands

```bash
# Debug cache issues
python scripts/cli_enhanced_validation.py --no-cache --verbose src/

# Performance debugging
python scripts/cli_enhanced_validation.py --stats
python scripts/cli_enhanced_validation.py --parallel-workers 1 --verbose src/

# Configuration debugging
python scripts/cli_enhanced_validation.py --config /dev/null --verbose src/
```

## Integration Examples

See the [Integration Guide](../integration/index.md) for detailed examples of:
- [GitHub Actions workflows](../integration/github_actions.md)
- [Pre-commit hook setup](../integration/pre_commit_hooks.md)
- [IDE integration](../integration/ide_integration.md)
- [Development workflow](../integration/development_workflow.md) 