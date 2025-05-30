# Integration

Welcome to the CLI Enhanced Validation Engine Integration Guide! This section provides comprehensive documentation for integrating the validation engine into your development workflows, CI/CD pipelines, and development environments.

## 🔗 Integration Options

The CLI Enhanced Validation Engine is designed to integrate seamlessly with modern development workflows:

### **Continuous Integration**
- **[GitHub Actions](github_actions.md)** - Complete workflow templates for GitHub repositories
- **Jenkins** - Pipeline integration with Jenkins CI/CD
- **GitLab CI** - YAML configuration for GitLab pipelines
- **Azure DevOps** - Integration with Azure Pipelines

### **Development Workflows**
- **[Pre-commit Hooks](pre_commit_hooks.md)** - Automatic validation before commits
- **[Development Workflow](development_workflow.md)** - Team development integration patterns
- **Code Review** - Pull request and merge request validation

### **IDE and Editor Integration**
- **[IDE Integration](ide_integration.md)** - VS Code, PyCharm, and other editor support
- **Language Server** - Real-time validation in supported editors
- **Plugins** - Editor-specific plugins and extensions

## 🚀 Quick Integration

Get started with the most common integration patterns:

### 1. GitHub Actions (Recommended)
Add automated validation to your GitHub repository in 2 minutes:

```yaml
# .github/workflows/validate-docs.yml
name: Documentation Validation
on: [push, pull_request]
jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v4
        with:
          python-version: '3.9'
      - run: pip install vcf-cli-validation
      - run: vcf-validate src/ --json --github-annotations
```

### 2. Pre-commit Hook
Enable validation before every commit:

```yaml
# .pre-commit-config.yaml
repos:
  - repo: https://github.com/vcf-agent/cli-enhanced-validation
    rev: v1.0.0
    hooks:
      - id: cli-validation
        args: [--mode, quick]
```

### 3. Development Script
Add to your development workflow:

```bash
#!/bin/bash
# scripts/validate.sh
echo "🔍 Validating documentation..."
python scripts/cli_enhanced_validation.py src/ --mode comprehensive
```

## 📋 Integration Strategies

### For Small Teams (1-5 developers)
- **Pre-commit hooks** for immediate feedback
- **GitHub Actions** for pull request validation
- **Simple configuration** with default settings

### For Medium Teams (5-20 developers)
- **Multi-stage validation** (quick in pre-commit, comprehensive in CI)
- **Branch-specific rules** with configuration management
- **IDE integration** for real-time feedback

### For Large Teams/Enterprise (20+ developers)
- **Advanced CI/CD integration** with multiple environments
- **Custom validation rules** and enterprise configuration
- **Monitoring and reporting** with detailed analytics

## 🛠️ Integration Components

### Core Integration Features

| Feature | Description | Use Case |
|---------|-------------|----------|
| **Exit Codes** | Standard success/failure codes | CI/CD pipeline control |
| **JSON Output** | Machine-readable results | Automation and reporting |
| **GitHub Annotations** | Inline PR comments | Code review integration |
| **Incremental Validation** | Changed files only | Performance optimization |
| **Configuration Management** | Environment-specific settings | Multi-environment deployment |

### Advanced Integration Features

| Feature | Description | Use Case |
|---------|-------------|----------|
| **Cache Management** | Persistent validation cache | Performance across runs |
| **Git Integration** | Branch comparison and staging | Smart validation targeting |
| **Performance Monitoring** | Execution metrics and timing | Optimization and debugging |
| **Custom Reporting** | Flexible output formats | Enterprise reporting needs |

## 📊 Integration Patterns

### 1. **Fail-Fast Pattern**
Quick validation in development, comprehensive in CI:

```yaml
# Development (pre-commit)
validation:
  coverage_threshold: 80.0
  mode: quick
  timeout_seconds: 30

# CI/CD (GitHub Actions)  
validation:
  coverage_threshold: 95.0
  mode: comprehensive
  timeout_seconds: 300
```

### 2. **Progressive Validation Pattern**
Different rules for different stages:

```bash
# Feature branch: Relaxed rules
vcf-validate src/ --coverage-threshold 75 --mode quick

# Main branch: Standard rules
vcf-validate src/ --coverage-threshold 90 --mode comprehensive

# Release: Strict rules
vcf-validate src/ --coverage-threshold 98 --mode strict --require-examples
```

### 3. **Microservice Pattern**
Independent validation per service:

```yaml
# docker-compose.yml
services:
  validation:
    image: vcf-cli-validation:latest
    volumes:
      - ./src:/workspace/src
    command: vcf-validate /workspace/src --json
```

## 🔧 Configuration for Integration

### Environment-Specific Configuration

=== "Development"

    ```yaml
    # .cli-validation.dev.yml
    validation:
      coverage_threshold: 80.0
      modes: [google]
      strict_mode: false

    performance:
      parallel_workers: 2
      timeout_seconds: 60

    reporting:
      verbose: false
      show_progress: true
    ```

=== "CI/CD"

    ```yaml
    # .cli-validation.ci.yml
    validation:
      coverage_threshold: 95.0
      modes: [google, numpy, sphinx]
      strict_mode: true

    git:
      enabled: true
      base_branch: "main"

    reporting:
      format: "json"
      github_annotations: true
    ```

=== "Production"

    ```yaml
    # .cli-validation.prod.yml
    validation:
      coverage_threshold: 98.0
      modes: [google, numpy, sphinx]
      strict_mode: true
      require_examples: true

    cache:
      enabled: true
      ttl_hours: 168  # 1 week

    reporting:
      format: "json"
      verbose: true
    ```

## 🎯 Integration Best Practices

### 1. **Start Simple, Scale Up**
- Begin with basic GitHub Actions integration
- Add pre-commit hooks for immediate feedback
- Gradually increase coverage thresholds and rules

### 2. **Use Incremental Validation**
- Validate only changed files in development
- Full validation on main branch and releases
- Cache results for repeated validations

### 3. **Configure by Environment**
- Development: Fast feedback, lower thresholds
- CI/CD: Comprehensive validation, detailed reporting
- Production: Strict rules, enterprise compliance

### 4. **Monitor and Optimize**
- Track validation execution times
- Monitor cache hit rates
- Optimize worker counts and timeouts

## 📈 Integration Monitoring

### Key Metrics to Track

| Metric | Target | Purpose |
|--------|--------|---------|
| **Execution Time** | <30s for quick, <5min for comprehensive | Performance monitoring |
| **Cache Hit Rate** | >70% | Cache effectiveness |
| **Coverage Trend** | Increasing over time | Documentation quality |
| **Failure Rate** | <5% in development | Developer experience |

### Monitoring Setup

```bash
# Add to CI/CD pipeline for monitoring
vcf-validate src/ --stats --json > validation-metrics.json

# Parse metrics for monitoring systems
jq '.performance_metrics' validation-metrics.json
```

## 🔗 Quick Links

Ready to integrate? Choose your integration path:

<div class="grid cards" markdown>

-   :fontawesome-brands-github:{ .lg .middle } **[GitHub Actions](github_actions.md)**

    ---

    Complete CI/CD integration with GitHub

-   :material-git:{ .lg .middle } **[Pre-commit Hooks](pre_commit_hooks.md)**

    ---

    Validate before every commit

-   :material-code-braces:{ .lg .middle } **[Development Workflow](development_workflow.md)**

    ---

    Team development integration

-   :material-application:{ .lg .middle } **[IDE Integration](ide_integration.md)**

    ---

    Real-time validation in editors

</div>

## 🆘 Integration Support

Need help with integration?

- **Check specific integration guides** for detailed setup instructions
- **Review [Examples](../user_guide/examples.md)** for common patterns
- **See [Troubleshooting](../user_guide/troubleshooting.md)** for common issues
- **Visit [Configuration](../user_guide/configuration.md)** for detailed options

---

*The CLI Enhanced Validation Engine is designed to integrate seamlessly with your existing development workflows. Start with one integration and expand as needed.* 