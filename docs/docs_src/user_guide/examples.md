# Examples

Real-world usage examples and workflows for the CLI Enhanced Validation Engine.

## Quick Start Examples

### Basic Validation

```bash
# Validate current directory
python scripts/cli_enhanced_validation.py

# Validate specific CLI module
python scripts/cli_enhanced_validation.py src/vcf_agent/cli/
```

**Expected Output**:
```
============================================================
CLI Documentation Validation Report
============================================================
📊 Comprehensive validation
📁 Files analyzed: 6
🔍 Definitions found: 121
📝 Coverage: 88.4%
⏱️  Execution time: 0.44s
💾 Cache hit rate: 100.0%

✅ SUCCESS - Validation completed with minor warnings
```

### Development Workflow

```bash
# Quick check during development
python scripts/cli_enhanced_validation.py --mode quick src/vcf_agent/cli/main.py

# Verbose output for debugging
python scripts/cli_enhanced_validation.py --verbose src/vcf_agent/cli/
```

**Use Case**: Fast feedback loop during active development.

## Git Integration Workflows

### Pull Request Validation

```bash
# Validate only changed files in PR
python scripts/cli_enhanced_validation.py --incremental --base-branch main

# With GitHub Actions integration
python scripts/cli_enhanced_validation.py \
  --incremental \
  --json \
  --github-annotations \
  --base-branch main
```

**GitHub Actions Workflow**:
```yaml
name: CLI Documentation Validation
on: [pull_request]

jobs:
  validate-docs:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
      with:
        fetch-depth: 0  # Need full history for git diff
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    - name: Install dependencies
      run: pip install -r requirements.txt
    - name: Validate CLI Documentation
      run: |
        python scripts/cli_enhanced_validation.py \
          --incremental \
          --json \
          --github-annotations \
          --base-branch origin/main > validation.json
        cat validation.json
```

### Pre-commit Hook Integration

```bash
# Validate staged files only
python scripts/cli_enhanced_validation.py --staged-files --mode quick
```

**.pre-commit-config.yaml**:
```yaml
repos:
  - repo: local
    hooks:
      - id: cli-doc-validation
        name: CLI Documentation Validation
        entry: python scripts/cli_enhanced_validation.py --staged-files --mode quick
        language: system
        types: [python]
        pass_filenames: false
        stages: [commit]
```

**Pre-commit Setup**:
```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install

# Test the hook
pre-commit run cli-doc-validation --all-files
```

## CI/CD Integration Examples

### GitHub Actions Complete Workflow

```yaml
name: Documentation Quality Gate
on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

jobs:
  validate-documentation:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ['3.9', '3.10', '3.11']
        
    steps:
    - name: Checkout code
      uses: actions/checkout@v3
      with:
        fetch-depth: 0
        
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}
        
    - name: Cache dependencies
      uses: actions/cache@v3
      with:
        path: ~/.cache/pip
        key: ${{ runner.os }}-pip-${{ hashFiles('**/requirements.txt') }}
        
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install docstring_parser  # For enhanced validation
        
    - name: Full Documentation Validation
      if: github.event_name == 'push'
      run: |
        python scripts/cli_enhanced_validation.py \
          --mode comprehensive \
          --json \
          --github-annotations \
          src/ > validation-full.json
        echo "Full validation results:"
        cat validation-full.json
        
    - name: Incremental Validation (PR)
      if: github.event_name == 'pull_request'
      run: |
        python scripts/cli_enhanced_validation.py \
          --incremental \
          --base-branch origin/main \
          --json \
          --github-annotations > validation-incremental.json
        echo "Incremental validation results:"
        cat validation-incremental.json
        
    - name: Performance Monitoring
      run: |
        python scripts/cli_enhanced_validation.py --stats --json > performance.json
        echo "Performance metrics:"
        cat performance.json
        
    - name: Upload validation results
      uses: actions/upload-artifact@v3
      with:
        name: validation-results-${{ matrix.python-version }}
        path: |
          validation-*.json
          performance.json
```

### Jenkins Pipeline

```groovy
pipeline {
    agent any
    
    environment {
        PYTHON_ENV = "${WORKSPACE}/.venv"
    }
    
    stages {
        stage('Setup') {
            steps {
                sh '''
                    python -m venv ${PYTHON_ENV}
                    . ${PYTHON_ENV}/bin/activate
                    pip install -r requirements.txt
                '''
            }
        }
        
        stage('CLI Documentation Validation') {
            parallel {
                stage('Quick Validation') {
                    steps {
                        sh '''
                            . ${PYTHON_ENV}/bin/activate
                            python scripts/cli_enhanced_validation.py \
                              --mode quick \
                              --json \
                              src/vcf_agent/cli/ > quick-validation.json
                        '''
                        archiveArtifacts artifacts: 'quick-validation.json'
                    }
                }
                
                stage('Comprehensive Validation') {
                    steps {
                        sh '''
                            . ${PYTHON_ENV}/bin/activate
                            python scripts/cli_enhanced_validation.py \
                              --mode comprehensive \
                              --json \
                              --parallel-workers 4 \
                              src/ > comprehensive-validation.json
                        '''
                        archiveArtifacts artifacts: 'comprehensive-validation.json'
                    }
                }
            }
        }
        
        stage('Performance Metrics') {
            steps {
                sh '''
                    . ${PYTHON_ENV}/bin/activate
                    python scripts/cli_enhanced_validation.py \
                      --stats \
                      --json > performance-metrics.json
                '''
                publishHTML([
                    allowMissing: false,
                    alwaysLinkToLastBuild: true,
                    keepAll: true,
                    reportDir: '.',
                    reportFiles: 'performance-metrics.json',
                    reportName: 'CLI Validation Performance'
                ])
            }
        }
    }
    
    post {
        always {
            sh '. ${PYTHON_ENV}/bin/activate && python scripts/cli_enhanced_validation.py --clear-cache'
        }
    }
}
```

## Development Team Workflows

### Code Review Process

```bash
# Reviewer validates PR changes
python scripts/cli_enhanced_validation.py \
  --incremental \
  --base-branch origin/main \
  --verbose

# Check specific files mentioned in PR
python scripts/cli_enhanced_validation.py \
  --mode comprehensive \
  --verbose \
  src/vcf_agent/cli/commands/validate.py
```

**Code Review Checklist**:
1. Run incremental validation
2. Check coverage percentage
3. Review new documentation quality
4. Verify no regression in existing docs

### Release Preparation

```bash
# Complete validation before release
python scripts/cli_enhanced_validation.py \
  --mode strict \
  --parallel-workers 4 \
  --verbose \
  src/

# Generate documentation quality report
python scripts/cli_enhanced_validation.py \
  --json \
  --mode comprehensive \
  src/ > release-docs-report.json
```

**Release Quality Gates**:
- ✅ 95%+ documentation coverage
- ✅ Zero critical issues
- ✅ All CLI components documented
- ✅ Performance within acceptable limits

## Performance Optimization Examples

### Cache Management

```bash
# Check cache performance
python scripts/cli_enhanced_validation.py --stats

# Clear cache if performance degrades
python scripts/cli_enhanced_validation.py --clear-cache

# Run without cache for baseline
python scripts/cli_enhanced_validation.py --no-cache src/vcf_agent/cli/
```

**Performance Comparison**:
```bash
# First run (cold cache)
time python scripts/cli_enhanced_validation.py src/
# Result: 2.34s

# Second run (warm cache)  
time python scripts/cli_enhanced_validation.py src/
# Result: 0.89s (62% improvement)
```

### Large Codebase Optimization

```bash
# Parallel processing for large repositories
python scripts/cli_enhanced_validation.py \
  --parallel-workers 8 \
  --timeout 600 \
  --mode comprehensive \
  src/

# Incremental validation for regular checks
python scripts/cli_enhanced_validation.py \
  --incremental \
  --mode quick \
  --parallel-workers 4
```

**Scaling Recommendations**:
- **Small projects** (<50 files): Default settings
- **Medium projects** (50-200 files): `--parallel-workers 4`
- **Large projects** (200+ files): `--parallel-workers 8` + incremental

## Configuration Examples

### Basic Configuration

**.cli-validation.yml**:
```yaml
# Basic configuration for small teams
validation:
  mode: comprehensive
  coverage_threshold: 90.0
  
cache:
  enabled: true
  ttl_hours: 24
  
output:
  json: false
  verbose: false
```

### Enterprise Configuration

**.cli-validation.yml**:
```yaml
# Enterprise configuration with strict requirements
validation:
  mode: strict
  coverage_threshold: 95.0
  parallel_workers: 8
  timeout_seconds: 300
  
cache:
  enabled: true
  ttl_hours: 168  # 1 week
  max_size_mb: 500
  
output:
  json: true
  verbose: true
  github_annotations: true
  
rules:
  require_examples: true
  validate_types: true
  check_parameter_docs: true
  
performance:
  incremental_by_default: true
  cache_optimization: true
```

### Project-specific Configuration

**pyproject.toml**:
```toml
[tool.cli-validation]
mode = "comprehensive"
coverage_threshold = 92.0

[tool.cli-validation.paths]
include = ["src/vcf_agent/cli/"]
exclude = ["tests/", "scripts/"]

[tool.cli-validation.rules]
require_return_docs = true
require_parameter_docs = true
allow_missing_examples = false

[tool.cli-validation.cache]
enabled = true
directory = ".validation-cache"
```

## Integration with IDEs

### VS Code Integration

**.vscode/tasks.json**:
```json
{
    "version": "2.0.0",
    "tasks": [
        {
            "label": "Validate CLI Documentation",
            "type": "shell",
            "command": "python",
            "args": [
                "scripts/cli_enhanced_validation.py",
                "--mode", "quick",
                "--verbose",
                "${file}"
            ],
            "group": "test",
            "presentation": {
                "echo": true,
                "reveal": "always",
                "focus": false,
                "panel": "shared"
            },
            "problemMatcher": {
                "owner": "cli-validation",
                "fileLocation": ["relative", "${workspaceFolder}"],
                "pattern": {
                    "regexp": "^(.*):(\\d+):(\\d+):\\s+(warning|error):\\s+(.*)$",
                    "file": 1,
                    "line": 2,
                    "column": 3,
                    "severity": 4,
                    "message": 5
                }
            }
        },
        {
            "label": "Full CLI Validation",
            "type": "shell",
            "command": "python",
            "args": [
                "scripts/cli_enhanced_validation.py",
                "--mode", "comprehensive",
                "src/vcf_agent/cli/"
            ],
            "group": "test"
        }
    ]
}
```

### PyCharm Integration

**External Tool Configuration**:
- **Name**: CLI Doc Validation
- **Program**: `python`
- **Arguments**: `scripts/cli_enhanced_validation.py --mode quick --verbose $FilePath$`
- **Working Directory**: `$ProjectFileDir$`

## Monitoring and Reporting

### Continuous Monitoring

```bash
# Daily documentation quality check
#!/bin/bash
# daily-doc-check.sh

REPORT_FILE="daily-docs-$(date +%Y%m%d).json"

python scripts/cli_enhanced_validation.py \
  --mode comprehensive \
  --json \
  src/ > "$REPORT_FILE"

# Extract key metrics
COVERAGE=$(cat "$REPORT_FILE" | jq '.validation_results.coverage_percentage')
ISSUES=$(cat "$REPORT_FILE" | jq '.issues | length')

echo "Daily Documentation Report - $(date)"
echo "Coverage: ${COVERAGE}%"
echo "Issues: ${ISSUES}"

# Alert if coverage drops below threshold
if (( $(echo "$COVERAGE < 90.0" | bc -l) )); then
    echo "⚠️  WARNING: Documentation coverage below 90%"
    # Send alert to team channel
fi
```

### Team Dashboard Integration

```bash
# Generate team metrics
python scripts/cli_enhanced_validation.py \
  --stats \
  --json > team-metrics.json

# Extract metrics for dashboard
jq '{
  coverage: .validation_results.coverage_percentage,
  total_files: .validation_results.total_files,
  execution_time: .metadata.execution_time,
  cache_efficiency: .metadata.cache_stats.hit_rate_percent
}' team-metrics.json > dashboard-data.json
```

## Testing and Validation

### Automated Testing Integration

```python
# test_documentation_quality.py
import subprocess
import json
import pytest

def test_cli_documentation_coverage():
    """Test that CLI documentation meets coverage requirements."""
    result = subprocess.run([
        'python', 'scripts/cli_enhanced_validation.py',
        '--mode', 'comprehensive',
        '--json',
        'src/vcf_agent/cli/'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Validation failed: {result.stderr}"
    
    report = json.loads(result.stdout)
    coverage = report['validation_results']['coverage_percentage']
    
    assert coverage >= 90.0, f"Coverage {coverage}% below 90% threshold"

def test_cli_documentation_performance():
    """Test that validation completes within performance limits."""
    result = subprocess.run([
        'python', 'scripts/cli_enhanced_validation.py',
        '--mode', 'quick',
        '--json',
        'src/vcf_agent/cli/'
    ], capture_output=True, text=True)
    
    report = json.loads(result.stdout)
    execution_time = report['metadata']['execution_time']
    
    assert execution_time < 5.0, f"Validation took {execution_time}s (>5s limit)"

if __name__ == "__main__":
    pytest.main([__file__])
```

## Advanced Use Cases

### Multi-repository Validation

```bash
# Validate multiple repositories
#!/bin/bash
# multi-repo-validation.sh

REPOS=("vcf-agent" "cli-tools" "validation-engine")
RESULTS_DIR="validation-results-$(date +%Y%m%d)"

mkdir -p "$RESULTS_DIR"

for repo in "${REPOS[@]}"; do
    echo "Validating $repo..."
    cd "$repo"
    
    python scripts/cli_enhanced_validation.py \
      --mode comprehensive \
      --json \
      src/ > "../$RESULTS_DIR/$repo-validation.json"
    
    cd ..
done

# Aggregate results
python tools/aggregate-validation-results.py "$RESULTS_DIR"/*.json
```

### Custom Validation Rules

```python
# custom-validation.py
from scripts.cli_enhanced_validation import EnhancedCLIValidator

class CustomCLIValidator(EnhancedCLIValidator):
    def validate_custom_rules(self, code_definition):
        """Add custom validation rules specific to your project."""
        issues = []
        
        # Custom rule: CLI commands must have examples
        if code_definition.is_cli_handler:
            if not self.has_examples(code_definition.docstring):
                issues.append({
                    'type': 'missing_examples',
                    'severity': 'warning',
                    'message': 'CLI command missing usage examples'
                })
        
        return issues

# Usage
validator = CustomCLIValidator()
report = validator.validate_comprehensive(['src/'])
```

This comprehensive examples section provides real-world scenarios that developers and teams encounter, from basic usage through enterprise deployment patterns. Each example includes expected outputs, integration code, and practical implementation details.
