# CLI Enhanced Validation Engine

[![Documentation Status](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://vcf-agent.github.io/cli-validation-docs)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Production-ready CLI documentation validation engine with comprehensive testing and automation features.**

The CLI Enhanced Validation Engine is a sophisticated tool that integrates existing CLI documentation validation capabilities with advanced AST parsing, multi-format docstring support, and CI/CD integration features. Built from extensive research using industry best practices, it provides enterprise-grade validation for Python CLI applications.

## ✨ Key Features

### 🔍 **Comprehensive Validation**
- **AST Analysis**: Deep code analysis with function, method, and class extraction
- **Multi-format Docstrings**: Support for Google, NumPy, and Sphinx formats
- **Coverage Analysis**: Detailed docstring coverage metrics and reporting
- **Structured Validation**: Parameter, return, and exception validation against function signatures

### ⚡ **Performance Optimized**
- **Incremental Validation**: Git integration for changed files only
- **Intelligent Caching**: Content-hash based caching with TTL support
- **Parallel Processing**: Multi-threaded analysis with configurable worker pools
- **Performance Monitoring**: Real-time metrics and cache statistics

### 🛠️ **Developer Experience**
- **Multiple Validation Modes**: Quick, comprehensive, and strict validation levels
- **Rich CLI Output**: Progress indicators, colored output, and detailed reporting
- **Configuration Management**: YAML, TOML, and environment variable support
- **IDE Integration**: Support for VS Code, PyCharm, and other editors

### 🚀 **CI/CD Ready**
- **GitHub Actions**: Pre-built workflow templates
- **Pre-commit Hooks**: Easy integration with development workflows
- **JSON Output**: Machine-readable results for automation
- **Exit Codes**: Proper status codes for pipeline integration

## 🚀 Quick Start

Get up and running in under 5 minutes:

### Installation

```bash
# Install from PyPI (recommended)
pip install vcf-cli-validation

# Or install from source
git clone https://github.com/vcf-agent/cli-enhanced-validation.git
cd cli-enhanced-validation
pip install -e .
```

### Basic Usage

```bash
# Validate a single file
python scripts/cli_enhanced_validation.py src/my_module.py

# Validate entire directory
python scripts/cli_enhanced_validation.py src/

# Quick validation mode
python scripts/cli_enhanced_validation.py src/ --mode quick

# Generate JSON report
python scripts/cli_enhanced_validation.py src/ --json
```

### First Validation

Run your first validation to see the power of the engine:

```bash
# Self-validation - validate the CLI tool itself
python scripts/cli_enhanced_validation.py scripts/cli_enhanced_validation.py

# Example output:
✅ SUCCESS - 6 definitions, 100.0% coverage, 0.08s execution
📊 Performance: Cache hit rate: 100.0%, Parallel workers: 4
```

## 📖 Documentation Sections

<div class="grid cards" markdown>

-   :fontawesome-solid-rocket:{ .lg .middle } **User Guide**

    ---

    Complete guides from installation to advanced usage

    [:octicons-arrow-right-24: Get Started](user_guide/quick_start.md)

-   :fontawesome-solid-cogs:{ .lg .middle } **Integration**

    ---

    CI/CD workflows, GitHub Actions, and development integration

    [:octicons-arrow-right-24: Setup Integration](integration/github_actions.md)

-   :fontawesome-solid-code:{ .lg .middle } **API Reference**

    ---

    Complete API documentation with examples

    [:octicons-arrow-right-24: Browse API](api/index.md)

-   :fontawesome-solid-tools:{ .lg .middle } **Advanced**

    ---

    Performance tuning, extending, and contributing

    [:octicons-arrow-right-24: Advanced Topics](advanced/performance_tuning.md)

</div>

## 🎯 Use Cases

### For Individual Developers
- **Code Quality**: Ensure comprehensive documentation coverage
- **Best Practices**: Follow Google, NumPy, and Sphinx docstring standards
- **IDE Integration**: Real-time validation in your development environment

### For Teams
- **Standardization**: Consistent documentation standards across projects
- **Code Reviews**: Automated validation in pull request workflows
- **Quality Gates**: Enforce documentation requirements in CI/CD pipelines

### For Enterprise
- **Compliance**: Meet documentation requirements for enterprise software
- **Automation**: Integrate with existing quality assurance processes
- **Reporting**: Generate detailed reports for stakeholders and audits

## 📊 Validation Results

The CLI Enhanced Validation Engine provides detailed insights into your codebase:

```bash
# Comprehensive analysis example
📁 Analyzing 28 files...
🔍 Found 712 definitions
📊 Coverage: 82.4% (587/712 definitions documented)
⚠️  889 warnings, 732 info messages
⏱️  Execution time: 2.06s

# Cache performance
💾 Cache statistics:
   Hit rate: 85.2%
   Entries: 234
   Size: 15.3 MB
```

## 🔗 Quick Links

- [Installation Guide](user_guide/quick_start.md#installation)
- [Configuration Reference](user_guide/configuration.md)
- [GitHub Actions Setup](integration/github_actions.md)
- [API Documentation](api/index.md)
- [Contributing Guide](advanced/contributing.md)

## 🏆 Why Choose CLI Enhanced Validation Engine?

Built from comprehensive research of industry best practices and real-world CLI tools, this engine provides:

- **Battle-tested Architecture**: Based on proven patterns from successful projects
- **Enterprise Ready**: Production-grade performance and reliability
- **Developer Friendly**: Intuitive interface with comprehensive documentation
- **Future Proof**: Extensible design supporting evolving validation needs

Ready to get started? Check out our [Quick Start Guide](user_guide/quick_start.md) or explore the [API Reference](api/index.md).

---

*Built with ❤️ by the VCF Agent Development Team* 