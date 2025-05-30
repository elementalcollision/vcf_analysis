# GitHub Actions Integration

Integrate the CLI Enhanced Validation Engine with GitHub Actions for automated documentation validation in your CI/CD pipeline. This guide provides complete workflow templates and best practices for different use cases.

## 🚀 Quick Setup

### Basic Workflow

Create `.github/workflows/validate-docs.yml` in your repository:

```yaml
name: Documentation Validation

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main ]

jobs:
  validate-docs:
    runs-on: ubuntu-latest
    
    steps:
    - name: Checkout code
      uses: actions/checkout@v4
      
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.9'
        
    - name: Install CLI Enhanced Validation Engine
      run: |
        pip install --upgrade pip
        pip install vcf-cli-validation
        
    - name: Validate documentation
      run: |
        python scripts/cli_enhanced_validation.py src/ \
          --mode comprehensive \
          --json \
          --github-annotations
```

### Expected Output

```bash
✅ Documentation validation passed
📊 Coverage: 95.2% (345/362 definitions documented)
🔍 Analyzed 28 files in 2.14s
💾 Cache hit rate: 73.5%
```

## 📋 Workflow Templates

### 1. **Standard Validation Workflow**

For most projects with standard validation requirements:

```yaml
name: Documentation Validation

on:
  push:
    branches: [ main, develop ]
    paths: 
      - '**.py'
      - '.cli-validation.yml'
      - '.github/workflows/validate-docs.yml'
  pull_request:
    branches: [ main ]
    paths:
      - '**.py'

jobs:
  validate-docs:
    runs-on: ubuntu-latest
    
    steps:
    - name: Checkout code
      uses: actions/checkout@v4
      with:
        fetch-depth: 0  # Needed for incremental validation
        
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.9'
        cache: 'pip'
        
    - name: Install dependencies
      run: |
        pip install --upgrade pip
        pip install vcf-cli-validation
        
    - name: Validate documentation (Incremental)
      if: github.event_name == 'pull_request'
      run: |
        python scripts/cli_enhanced_validation.py \
          --incremental \
          --base-branch origin/${{ github.base_ref }} \
          --json \
          --github-annotations
          
    - name: Validate documentation (Full)
      if: github.event_name == 'push'
      run: |
        python scripts/cli_enhanced_validation.py src/ \
          --mode comprehensive \
          --coverage-threshold 90 \
          --json
```

### 2. **Multi-Environment Workflow**

For projects with different validation rules per environment:

```yaml
name: Documentation Validation

on:
  push:
    branches: [ main, develop, 'feature/*' ]
  pull_request:
    branches: [ main, develop ]

jobs:
  validate-feature:
    if: startsWith(github.ref, 'refs/heads/feature/') || github.event_name == 'pull_request'
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v4
    - uses: actions/setup-python@v4
      with:
        python-version: '3.9'
        
    - name: Install CLI validation
      run: pip install vcf-cli-validation
      
    - name: Feature validation (Relaxed)
      run: |
        python scripts/cli_enhanced_validation.py src/ \
          --mode quick \
          --coverage-threshold 75 \
          --json
          
  validate-develop:
    if: github.ref == 'refs/heads/develop'
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v4
    - uses: actions/setup-python@v4
      with:
        python-version: '3.9'
        
    - name: Install CLI validation
      run: pip install vcf-cli-validation
      
    - name: Development validation (Standard)
      run: |
        python scripts/cli_enhanced_validation.py src/ \
          --mode comprehensive \
          --coverage-threshold 85 \
          --json
          
  validate-main:
    if: github.ref == 'refs/heads/main'
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v4
    - uses: actions/setup-python@v4
      with:
        python-version: '3.9'
        
    - name: Install CLI validation
      run: pip install vcf-cli-validation
      
    - name: Production validation (Strict)
      run: |
        python scripts/cli_enhanced_validation.py src/ \
          --mode strict \
          --coverage-threshold 95 \
          --require-examples \
          --json \
          --github-annotations
```

### 3. **Matrix Testing Workflow**

For testing across multiple Python versions:

```yaml
name: Documentation Validation Matrix

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  validate-matrix:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, windows-latest, macos-latest]
        python-version: ['3.9', '3.10', '3.11', '3.12']
        
    steps:
    - uses: actions/checkout@v4
    
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}
        
    - name: Install CLI validation
      run: |
        pip install --upgrade pip
        pip install vcf-cli-validation
        
    - name: Validate documentation
      run: |
        python scripts/cli_enhanced_validation.py src/ \
          --mode comprehensive \
          --json
          
    - name: Upload validation results
      uses: actions/upload-artifact@v3
      if: failure()
      with:
        name: validation-results-${{ matrix.os }}-${{ matrix.python-version }}
        path: validation-report.json
```

## 🎯 Advanced Features

### 1. **Caching for Performance**

Enable caching to speed up repeated validations:

```yaml
- name: Cache validation results
  uses: actions/cache@v3
  with:
    path: .validation-cache
    key: validation-${{ runner.os }}-${{ hashFiles('**/*.py') }}
    restore-keys: |
      validation-${{ runner.os }}-
      
- name: Validate with caching
  run: |
    python scripts/cli_enhanced_validation.py src/ \
      --cache-enabled \
      --cache-directory .validation-cache \
      --json
```

### 2. **GitHub Annotations**

Add inline comments to pull requests:

```yaml
- name: Validate with annotations
  run: |
    python scripts/cli_enhanced_validation.py src/ \
      --json \
      --github-annotations \
      --output-file validation-results.json
      
- name: Process annotations
  if: always()
  uses: actions/github-script@v6
  with:
    script: |
      const fs = require('fs');
      try {
        const results = JSON.parse(fs.readFileSync('validation-results.json', 'utf8'));
        if (results.github_annotations) {
          for (const annotation of results.github_annotations) {
            core.error(annotation.message, {
              file: annotation.file,
              startLine: annotation.line,
              title: annotation.title
            });
          }
        }
      } catch (error) {
        core.warning('Could not process validation annotations');
      }
```

### 3. **Conditional Validation**

Run different validation based on changed files:

```yaml
- name: Get changed files
  id: changed-files
  uses: tj-actions/changed-files@v40
  with:
    files: |
      **.py
      
- name: Validate changed files only
  if: steps.changed-files.outputs.any_changed == 'true'
  run: |
    echo "Changed Python files: ${{ steps.changed-files.outputs.all_changed_files }}"
    python scripts/cli_enhanced_validation.py \
      ${{ steps.changed-files.outputs.all_changed_files }} \
      --mode quick \
      --json
```

### 4. **Performance Monitoring**

Track validation performance over time:

```yaml
- name: Validate with metrics
  run: |
    python scripts/cli_enhanced_validation.py src/ \
      --stats \
      --json \
      --output-file metrics.json
      
- name: Upload metrics
  uses: actions/upload-artifact@v3
  with:
    name: validation-metrics
    path: metrics.json
    
- name: Report performance
  run: |
    echo "## Validation Performance" >> $GITHUB_STEP_SUMMARY
    echo "$(jq -r '.performance_metrics | to_entries[] | "- \(.key): \(.value)"' metrics.json)" >> $GITHUB_STEP_SUMMARY
```

## 🔧 Configuration Management

### Environment-Specific Configurations

Use different configurations for different environments:

```yaml
- name: Select configuration
  run: |
    if [ "${{ github.ref }}" == "refs/heads/main" ]; then
      cp .cli-validation.prod.yml .cli-validation.yml
    elif [ "${{ github.ref }}" == "refs/heads/develop" ]; then
      cp .cli-validation.staging.yml .cli-validation.yml
    else
      cp .cli-validation.dev.yml .cli-validation.yml
    fi
    
- name: Validate with environment config
  run: |
    python scripts/cli_enhanced_validation.py src/ --json
```

### Using GitHub Secrets

Store sensitive configuration in GitHub Secrets:

```yaml
- name: Configure validation
  env:
    CLI_VALIDATION_COVERAGE_THRESHOLD: ${{ secrets.COVERAGE_THRESHOLD }}
    CLI_VALIDATION_MODES: ${{ secrets.VALIDATION_MODES }}
  run: |
    python scripts/cli_enhanced_validation.py src/ --json
```

## 📊 Reporting and Artifacts

### 1. **Detailed Reporting**

Generate comprehensive reports:

```yaml
- name: Generate detailed report
  run: |
    python scripts/cli_enhanced_validation.py src/ \
      --mode comprehensive \
      --json \
      --output-file detailed-report.json
      
- name: Create summary report
  run: |
    echo "# Documentation Validation Report" > report.md
    echo "" >> report.md
    jq -r '
      "## Summary",
      "- **Coverage**: \(.validation_results.coverage_percentage)%",
      "- **Files analyzed**: \(.metadata.files_analyzed)",
      "- **Definitions found**: \(.validation_results.total_definitions)",
      "- **Execution time**: \(.metadata.execution_time_seconds)s",
      "",
      "## Issues Found",
      if (.validation_results.issues | length) > 0 then
        (.validation_results.issues[] | "- **\(.file):\(.line)**: \(.message)")
      else
        "✅ No issues found!"
      end
    ' detailed-report.json >> report.md
    
- name: Upload report
  uses: actions/upload-artifact@v3
  with:
    name: validation-report
    path: |
      detailed-report.json
      report.md
```

### 2. **Pull Request Comments**

Add validation results as PR comments:

```yaml
- name: Comment PR
  if: github.event_name == 'pull_request'
  uses: actions/github-script@v6
  with:
    script: |
      const fs = require('fs');
      const results = JSON.parse(fs.readFileSync('detailed-report.json', 'utf8'));
      
      const coverage = results.validation_results.coverage_percentage;
      const coverageEmoji = coverage >= 95 ? '✅' : coverage >= 85 ? '⚠️' : '❌';
      
      const comment = `## ${coverageEmoji} Documentation Validation Results
      
      **Coverage**: ${coverage}% (${results.validation_results.documented_definitions}/${results.validation_results.total_definitions} definitions)
      **Execution time**: ${results.metadata.execution_time_seconds}s
      **Cache hit rate**: ${results.metadata.cache_stats?.hit_rate_percent || 0}%
      
      ${results.validation_results.issues.length > 0 ? 
        `### Issues Found (${results.validation_results.issues.length})
        ${results.validation_results.issues.slice(0, 10).map(issue => 
          `- \`${issue.file}:${issue.line}\`: ${issue.message}`
        ).join('\n')}
        ${results.validation_results.issues.length > 10 ? '\n...and more' : ''}` :
        '✅ No issues found!'
      }`;
      
      github.rest.issues.createComment({
        issue_number: context.issue.number,
        owner: context.repo.owner,
        repo: context.repo.repo,
        body: comment
      });
```

## 🛡️ Security Best Practices

### 1. **Minimal Permissions**

Use minimal required permissions:

```yaml
permissions:
  contents: read        # Read repository contents
  pull-requests: write  # Comment on PRs
  checks: write        # Write check status
```

### 2. **Secure Dependencies**

Pin dependency versions and use trusted sources:

```yaml
- name: Install CLI validation (pinned)
  run: |
    pip install vcf-cli-validation==1.0.0 \
      --only-binary=all \
      --require-hashes
```

### 3. **Environment Isolation**

Use containers for consistent environments:

```yaml
jobs:
  validate-docs:
    runs-on: ubuntu-latest
    container: python:3.9-slim
    
    steps:
    - name: Install system dependencies
      run: |
        apt-get update && apt-get install -y git
        
    - uses: actions/checkout@v4
    
    - name: Install CLI validation
      run: pip install vcf-cli-validation
      
    - name: Validate documentation
      run: python scripts/cli_enhanced_validation.py src/ --json
```

## 🚨 Troubleshooting

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| **Permission denied** | Insufficient workflow permissions | Add required permissions to workflow |
| **Cache not working** | Incorrect cache key or path | Verify cache path and key patterns |
| **Annotations not appearing** | Missing github-annotations flag | Add `--github-annotations` to command |
| **Validation too slow** | Large codebase, no incremental mode | Use `--incremental` for PRs |

### Debug Mode

Enable verbose output for troubleshooting:

```yaml
- name: Debug validation
  run: |
    python scripts/cli_enhanced_validation.py src/ \
      --verbose \
      --debug \
      --json
```

### Performance Optimization

For large repositories:

```yaml
- name: Optimized validation
  run: |
    python scripts/cli_enhanced_validation.py src/ \
      --parallel-workers 8 \
      --cache-enabled \
      --incremental \
      --json
```

## 📈 Monitoring and Analytics

### Workflow Metrics

Track validation trends over time:

```yaml
- name: Send metrics to monitoring
  if: always()
  run: |
    curl -X POST "${{ secrets.MONITORING_WEBHOOK }}" \
      -H "Content-Type: application/json" \
      -d @metrics.json
```

### Success/Failure Notifications

Set up notifications for validation results:

```yaml
- name: Notify on failure
  if: failure()
  uses: 8398a7/action-slack@v3
  with:
    status: failure
    channel: '#dev-alerts'
    webhook_url: ${{ secrets.SLACK_WEBHOOK }}
```

## 🎉 Complete Example

Here's a production-ready workflow combining all best practices:

```yaml
name: Documentation Validation

on:
  push:
    branches: [ main, develop ]
    paths: [ '**.py', '.cli-validation.yml' ]
  pull_request:
    branches: [ main ]
    paths: [ '**.py' ]

permissions:
  contents: read
  pull-requests: write
  checks: write

jobs:
  validate-docs:
    runs-on: ubuntu-latest
    
    steps:
    - name: Checkout code
      uses: actions/checkout@v4
      with:
        fetch-depth: 0
        
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.9'
        cache: 'pip'
        
    - name: Cache validation results
      uses: actions/cache@v3
      with:
        path: .validation-cache
        key: validation-${{ runner.os }}-${{ hashFiles('**/*.py') }}
        
    - name: Install CLI Enhanced Validation Engine
      run: |
        pip install --upgrade pip
        pip install vcf-cli-validation==1.0.0
        
    - name: Validate documentation (PR)
      if: github.event_name == 'pull_request'
      run: |
        python scripts/cli_enhanced_validation.py \
          --incremental \
          --base-branch origin/${{ github.base_ref }} \
          --cache-enabled \
          --json \
          --github-annotations \
          --output-file results.json
          
    - name: Validate documentation (Push)
      if: github.event_name == 'push'
      run: |
        python scripts/cli_enhanced_validation.py src/ \
          --mode comprehensive \
          --coverage-threshold 90 \
          --cache-enabled \
          --json \
          --output-file results.json
          
    - name: Comment PR
      if: github.event_name == 'pull_request' && always()
      uses: actions/github-script@v6
      with:
        script: |
          const fs = require('fs');
          try {
            const results = JSON.parse(fs.readFileSync('results.json', 'utf8'));
            const coverage = results.validation_results.coverage_percentage;
            const emoji = coverage >= 95 ? '✅' : coverage >= 85 ? '⚠️' : '❌';
            
            const comment = `## ${emoji} Documentation Validation
            **Coverage**: ${coverage}% | **Time**: ${results.metadata.execution_time_seconds}s`;
            
            github.rest.issues.createComment({
              issue_number: context.issue.number,
              owner: context.repo.owner,
              repo: context.repo.repo,
              body: comment
            });
          } catch (error) {
            console.log('Could not create PR comment:', error);
          }
          
    - name: Upload results
      if: always()
      uses: actions/upload-artifact@v3
      with:
        name: validation-results
        path: results.json
```

---

Ready to set up GitHub Actions? Start with the basic workflow and gradually add advanced features as needed. Check the [Pre-commit Hooks](pre_commit_hooks.md) guide for local validation before pushing to GitHub! 