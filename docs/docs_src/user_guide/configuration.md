# Configuration

The CLI Enhanced Validation Engine provides flexible configuration through multiple sources with a clear precedence hierarchy. This guide covers all configuration options, formats, and best practices.

## Configuration Sources

The validation engine supports configuration from multiple sources in order of precedence:

1. **Command Line Arguments** (highest priority)
2. **Environment Variables**
3. **YAML Configuration Files** (`.cli-validation.yml`)
4. **TOML Configuration Files** (`.cli-validation.toml`)
5. **Default Values** (lowest priority)

## Configuration Files

### YAML Configuration

Create `.cli-validation.yml` in your project root:

```yaml
# CLI Enhanced Validation Engine Configuration
# Complete configuration reference

# Validation settings
validation:
  coverage_threshold: 95.0           # Minimum coverage percentage
  modes: [google, numpy, sphinx]     # Supported docstring formats
  require_examples: false            # Require examples in docstrings
  validate_parameter_types: true     # Validate parameter type annotations
  validate_return_types: true        # Validate return type annotations
  strict_mode: false                 # Enable strict validation rules

# Cache configuration
cache:
  enabled: true                      # Enable caching system
  directory: ".validation-cache"     # Cache directory location
  ttl_hours: 24                     # Time-to-live in hours
  max_size_mb: 100                  # Maximum cache size in MB
  compression: true                  # Enable cache compression

# Performance settings
performance:
  parallel_workers: 4                # Number of parallel workers
  chunk_size: 10                    # Files per worker chunk
  max_memory_mb: 512                # Maximum memory usage
  timeout_seconds: 300              # Maximum execution timeout

# Git integration
git:
  enabled: true                      # Enable Git integration
  base_branch: "main"               # Default base branch
  ignore_staged: false              # Include staged files in incremental mode
  exclude_patterns:                 # Files to exclude from Git analysis
    - "*.pyc"
    - "__pycache__/*"
    - "*.egg-info/*"

# Reporting options
reporting:
  format: "console"                 # Output format: console, json
  verbose: false                    # Enable verbose output
  colors: true                      # Enable colored output
  show_progress: true               # Show progress indicators
  github_annotations: false        # Generate GitHub Actions annotations

# File processing
files:
  include_patterns:                 # File patterns to include
    - "*.py"
  exclude_patterns:                 # File patterns to exclude
    - "test_*.py"
    - "*_test.py"
    - "tests/*"
  max_file_size_kb: 1024           # Maximum file size to process
```

### TOML Configuration

Create `.cli-validation.toml` in your project root:

```toml
# CLI Enhanced Validation Engine Configuration
# Complete configuration reference

[validation]
coverage_threshold = 95.0
modes = ["google", "numpy", "sphinx"]
require_examples = false
validate_parameter_types = true
validate_return_types = true
strict_mode = false

[cache]
enabled = true
directory = ".validation-cache"
ttl_hours = 24
max_size_mb = 100
compression = true

[performance]
parallel_workers = 4
chunk_size = 10
max_memory_mb = 512
timeout_seconds = 300

[git]
enabled = true
base_branch = "main"
ignore_staged = false
exclude_patterns = ["*.pyc", "__pycache__/*", "*.egg-info/*"]

[reporting]
format = "console"
verbose = false
colors = true
show_progress = true
github_annotations = false

[files]
include_patterns = ["*.py"]
exclude_patterns = ["test_*.py", "*_test.py", "tests/*"]
max_file_size_kb = 1024
```

## Environment Variables

All configuration options can be set via environment variables using the prefix `CLI_VALIDATION_`:

```bash
# Validation settings
export CLI_VALIDATION_COVERAGE_THRESHOLD=90.0
export CLI_VALIDATION_MODES="google,numpy"
export CLI_VALIDATION_STRICT_MODE=true

# Cache settings
export CLI_VALIDATION_CACHE_ENABLED=true
export CLI_VALIDATION_CACHE_TTL_HOURS=48
export CLI_VALIDATION_CACHE_DIRECTORY="/tmp/validation-cache"

# Performance settings
export CLI_VALIDATION_PARALLEL_WORKERS=8
export CLI_VALIDATION_MAX_MEMORY_MB=1024

# Git settings
export CLI_VALIDATION_GIT_BASE_BRANCH="develop"
export CLI_VALIDATION_GIT_IGNORE_STAGED=true

# Reporting settings
export CLI_VALIDATION_REPORTING_FORMAT="json"
export CLI_VALIDATION_REPORTING_VERBOSE=true
```

## Configuration Options Reference

### Validation Settings

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `coverage_threshold` | float | 80.0 | Minimum docstring coverage percentage |
| `modes` | list | ["google"] | Supported docstring formats |
| `require_examples` | bool | false | Require examples in docstrings |
| `validate_parameter_types` | bool | true | Validate parameter type annotations |
| `validate_return_types` | bool | true | Validate return type annotations |
| `strict_mode` | bool | false | Enable strictest validation rules |

### Cache Settings

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `enabled` | bool | true | Enable caching system |
| `directory` | string | ".validation-cache" | Cache directory location |
| `ttl_hours` | int | 24 | Time-to-live in hours |
| `max_size_mb` | int | 100 | Maximum cache size in MB |
| `compression` | bool | true | Enable cache compression |

### Performance Settings

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `parallel_workers` | int | 4 | Number of parallel workers |
| `chunk_size` | int | 10 | Files per worker chunk |
| `max_memory_mb` | int | 512 | Maximum memory usage |
| `timeout_seconds` | int | 300 | Maximum execution timeout |

### Git Integration

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `enabled` | bool | true | Enable Git integration |
| `base_branch` | string | "main" | Default base branch for comparisons |
| `ignore_staged` | bool | false | Include staged files in incremental mode |
| `exclude_patterns` | list | ["*.pyc"] | File patterns to exclude |

### Reporting Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `format` | string | "console" | Output format: console, json |
| `verbose` | bool | false | Enable verbose output |
| `colors` | bool | true | Enable colored output |
| `show_progress` | bool | true | Show progress indicators |
| `github_annotations` | bool | false | Generate GitHub Actions annotations |

### File Processing

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `include_patterns` | list | ["*.py"] | File patterns to include |
| `exclude_patterns` | list | ["test_*.py"] | File patterns to exclude |
| `max_file_size_kb` | int | 1024 | Maximum file size to process |

## Configuration Examples

### Development Workflow

For day-to-day development with fast feedback:

```yaml
validation:
  coverage_threshold: 85.0
  modes: [google]
  strict_mode: false

cache:
  enabled: true
  ttl_hours: 12

performance:
  parallel_workers: 2
  timeout_seconds: 60

reporting:
  verbose: false
  show_progress: true
```

### CI/CD Pipeline

For automated CI/CD validation:

```yaml
validation:
  coverage_threshold: 95.0
  modes: [google, numpy, sphinx]
  strict_mode: true
  require_examples: true

cache:
  enabled: true
  ttl_hours: 168  # 1 week

performance:
  parallel_workers: 8
  timeout_seconds: 600

git:
  enabled: true
  base_branch: "main"

reporting:
  format: "json"
  github_annotations: true
  verbose: true
```

### Enterprise/Compliance

For enterprise environments with strict requirements:

```yaml
validation:
  coverage_threshold: 98.0
  modes: [google, numpy, sphinx]
  strict_mode: true
  require_examples: true
  validate_parameter_types: true
  validate_return_types: true

cache:
  enabled: true
  directory: "/opt/validation-cache"
  ttl_hours: 720  # 30 days
  max_size_mb: 500

performance:
  parallel_workers: 16
  max_memory_mb: 2048
  timeout_seconds: 1800

files:
  exclude_patterns:
    - "test_*.py"
    - "*_test.py"
    - "tests/*"
    - "*.pyc"
    - "__pycache__/*"
    - "migrations/*"

reporting:
  format: "json"
  verbose: true
  colors: false
```

## Configuration Validation

The CLI Enhanced Validation Engine validates configuration files on startup:

```bash
# Test configuration validity
python scripts/cli_enhanced_validation.py --config-test

# Show current configuration
python scripts/cli_enhanced_validation.py --show-config
```

### Common Configuration Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `Invalid coverage_threshold` | Value outside 0.0-100.0 range | Set value between 0.0 and 100.0 |
| `Unknown validation mode` | Unsupported docstring format | Use: google, numpy, or sphinx |
| `Invalid parallel_workers` | Non-positive integer | Set to positive integer (1-16) |
| `Cache directory not writable` | Permission issues | Check directory permissions |

## Configuration Precedence

When the same option is specified in multiple sources, the CLI Enhanced Validation Engine follows this precedence order:

1. **Command Line Arguments** (highest)
   ```bash
   python scripts/cli_enhanced_validation.py --coverage-threshold 90.0
   ```

2. **Environment Variables**
   ```bash
   export CLI_VALIDATION_COVERAGE_THRESHOLD=85.0
   ```

3. **YAML Configuration File**
   ```yaml
   validation:
     coverage_threshold: 80.0
   ```

4. **TOML Configuration File**
   ```toml
   [validation]
   coverage_threshold = 75.0
   ```

5. **Default Values** (lowest)

### Example: Precedence in Action

Given this configuration scenario:

- **Config file**: `coverage_threshold: 80.0`
- **Environment**: `CLI_VALIDATION_COVERAGE_THRESHOLD=85.0`
- **Command line**: `--coverage-threshold 90.0`

The final value used will be **90.0** (command line takes precedence).

## Advanced Configuration

### Custom File Patterns

Include specific file types and exclude test files:

```yaml
files:
  include_patterns:
    - "src/**/*.py"
    - "lib/**/*.py"
    - "apps/**/*.py"
  exclude_patterns:
    - "test_*.py"
    - "*_test.py"
    - "tests/**"
    - "conftest.py"
    - "**/migrations/*.py"
```

### Performance Tuning

Optimize for large codebases:

```yaml
performance:
  parallel_workers: 8           # 2x CPU cores
  chunk_size: 20               # Larger chunks for big projects
  max_memory_mb: 1024          # Higher memory limit
  timeout_seconds: 1800        # 30 minutes for large projects

cache:
  max_size_mb: 500            # Larger cache for big projects
  compression: true           # Save disk space
```

### Multi-Environment Configuration

Use different configurations for different environments:

```bash
# Development
cp .cli-validation.dev.yml .cli-validation.yml

# Staging
cp .cli-validation.staging.yml .cli-validation.yml

# Production
cp .cli-validation.prod.yml .cli-validation.yml
```

## Troubleshooting Configuration

### Debug Configuration Loading

```bash
# Show which configuration files are being loaded
python scripts/cli_enhanced_validation.py --verbose

# Test specific configuration file
python scripts/cli_enhanced_validation.py --config custom-config.yml --config-test
```

### Common Issues

!!! warning "Configuration File Not Found"
    **Issue**: Configuration file not loaded
    
    **Solutions**:
    - Check file exists in current directory
    - Verify file permissions (readable)
    - Check YAML/TOML syntax validity

!!! warning "Environment Variables Not Working"
    **Issue**: Environment variables ignored
    
    **Solutions**:
    - Use correct prefix: `CLI_VALIDATION_`
    - Check variable names match configuration keys
    - Export variables in current shell session

!!! warning "Performance Issues"
    **Issue**: Validation too slow
    
    **Solutions**:
    - Increase `parallel_workers`
    - Enable `cache.enabled`
    - Reduce `coverage_threshold` for development
    - Use `--mode quick` for faster validation

## Configuration Best Practices

### 1. Environment-Specific Configurations

- **Development**: Lower thresholds, faster execution
- **CI/CD**: Comprehensive validation, JSON output
- **Production**: Strict rules, detailed reporting

### 2. Cache Management

- Enable caching for repeated validations
- Set appropriate TTL based on code change frequency
- Monitor cache size and clean periodically

### 3. Git Integration

- Use incremental validation for large projects
- Set appropriate base branch for comparisons
- Exclude irrelevant file patterns

### 4. Performance Optimization

- Tune parallel workers based on CPU cores
- Adjust chunk size for project size
- Set memory limits to prevent system issues

---

Ready to customize your validation? Start with the [Examples](examples.md) guide to see configuration in action, or check the [CLI Reference](cli_reference.md) for command-line options. 