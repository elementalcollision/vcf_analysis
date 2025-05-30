# VCF Analysis Agent User Guide

## Overview

The VCF Analysis Agent is a powerful, AI-driven tool for analyzing Variant Call Format (VCF) files with advanced genomic insights. This guide will help you get started with installation, configuration, and common usage patterns.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Configuration](#configuration)
4. [Basic Usage](#basic-usage)
5. [Advanced Features](#advanced-features)
6. [Data Management](#data-management)
7. [AI-Powered Analysis](#ai-powered-analysis)
8. [Performance Optimization](#performance-optimization)
9. [Troubleshooting](#troubleshooting)
10. [Best Practices](#best-practices)

## Quick Start

### 5-Minute Setup

1. **Clone and Setup**:
   ```bash
   git clone https://github.com/your-org/vcf-analysis-agent.git
   cd vcf-analysis-agent
   pip install -e .
   ```

2. **Configure API Keys**:
   ```bash
   export OPENAI_API_KEY="your-openai-key"
   export CEREBRAS_API_KEY="your-cerebras-key"  # Optional
   ```

3. **Analyze Your First VCF**:
   ```bash
   vcf-agent analyze sample_data/example.vcf.gz
   ```

4. **Start Interactive Analysis**:
   ```bash
   vcf-agent ask "What are the most significant variants in this file?"
   ```

## Installation

### Prerequisites

- **Python**: 3.11 or higher
- **Operating System**: macOS, Linux, or Windows with WSL
- **Memory**: 4GB RAM minimum, 8GB recommended
- **Storage**: 2GB free space for databases and temporary files

### Installation Methods

#### Method 1: pip Install (Recommended)

```bash
pip install vcf-analysis-agent
```

#### Method 2: From Source

```bash
git clone https://github.com/your-org/vcf-analysis-agent.git
cd vcf-analysis-agent
pip install -e .
```

#### Method 3: Docker

```bash
docker pull vcf-agent:latest
docker run -it --rm -v $(pwd):/data vcf-agent:latest
```

### Verify Installation

```bash
vcf-agent --version
vcf-agent --help
```

## Configuration

### Environment Variables

Create a `.env` file in your working directory:

```bash
# AI Model Configuration
OPENAI_API_KEY=sk-your-openai-key-here
CEREBRAS_API_KEY=your-cerebras-key-here  # Optional
OLLAMA_BASE_URL=http://localhost:11434   # For local models

# Database Configuration
LANCEDB_PATH=./data/lancedb
KUZU_PATH=./data/kuzu_db
LANCEDB_TABLE_NAME=vcf_variants

# Performance Settings
MAX_WORKERS=4
BATCH_SIZE=1000
ENABLE_PERFORMANCE_MONITORING=true

# Logging Configuration
LOG_LEVEL=INFO
STRUCTURED_LOGGING=true
```

### Configuration File

Create `config/vcf_agent_config.yaml`:

```yaml
# VCF Agent Configuration
session:
  model_provider: "openai"  # openai, cerebras, ollama
  default_model: "gpt-4"
  temperature: 0.1
  max_tokens: 4000

data_stores:
  lancedb_path: "./data/lancedb"
  kuzu_path: "./data/kuzu_db"
  enable_sync: true
  batch_size: 1000

performance:
  max_workers: 4
  enable_monitoring: true
  cache_embeddings: true

security:
  mask_sensitive_data: true
  enable_audit_logging: true
```

### API Key Management

#### Option 1: Environment Variables
```bash
export OPENAI_API_KEY="your-key"
export CEREBRAS_API_KEY="your-key"
```

#### Option 2: Credentials File
Create `credentials.json`:
```json
{
  "openai_api_key": "your-openai-key",
  "cerebras_api_key": "your-cerebras-key"
}
```

#### Option 3: Secure Vault (Production)
```bash
# Using HashiCorp Vault
vault kv put secret/vcf-agent openai_api_key="your-key"

# Using Kubernetes Secrets
kubectl create secret generic vcf-agent-secrets \
  --from-literal=openai-api-key="your-key"
```

## Basic Usage

### Command Line Interface

#### Analyze VCF Files

```bash
# Basic analysis
vcf-agent analyze sample.vcf.gz

# With specific output directory
vcf-agent analyze sample.vcf.gz --output-dir ./results

# Multiple files
vcf-agent analyze *.vcf.gz --batch-mode

# With custom filters
vcf-agent analyze sample.vcf.gz --filter "QUAL > 30 && INFO/DP > 10"
```

#### Interactive AI Analysis

```bash
# Start interactive session
vcf-agent ask "Summarize the variants in sample.vcf.gz"

# Specific questions
vcf-agent ask "What are the pathogenic variants in BRCA1?"
vcf-agent ask "Compare variant quality between samples A and B"
vcf-agent ask "Identify potential compound heterozygotes"
```

#### Data Management

```bash
# Ingest VCF into databases
vcf-agent ingest-vcf sample.vcf.gz --sample-id SAMPLE_001

# Search variants
vcf-agent search "pathogenic BRCA1 variants"
vcf-agent search --gene BRCA1 --clinical-significance pathogenic

# Export results
vcf-agent export --format json --output results.json
vcf-agent export --format vcf --filter "gene=BRCA1" --output brca1_variants.vcf
```

### Python API

#### Basic Usage

```python
from vcf_agent import VCFAnalysisAgent
from vcf_agent.config import SessionConfig

# Initialize agent
config = SessionConfig(
    model_provider="openai",
    credentials_file="credentials.json"
)
agent = VCFAnalysisAgent(config)

# Analyze VCF file
results = agent.analyze_vcf("sample.vcf.gz")
print(f"Found {results['variant_count']} variants")

# Ask questions
response = agent.ask("What are the most significant variants?")
print(response.content)
```

#### Advanced Usage

```python
from vcf_agent.data_store_manager import create_data_store_manager
from vcf_agent.lancedb_integration import VCFVariant

# Create data store manager
manager = create_data_store_manager(
    lancedb_path="./data/lancedb",
    kuzu_path="./data/kuzu_db"
)

# Add sample data
sample_data = {
    "id": "SAMPLE_001",
    "name": "Patient Sample",
    "type": "germline"
}

variants_data = [
    {
        "id": "chr1-123456-A-G",
        "chr": "1",
        "pos": 123456,
        "ref": "A",
        "alt": "G",
        "gene_symbol": "BRCA1",
        "clinical_significance": "Pathogenic"
    }
]

# Ingest data
result = manager.add_sample_with_variants(
    sample_data=sample_data,
    variants_data=variants_data
)

# Search variants
search_results = manager.search_variants(
    query="pathogenic variants in BRCA1",
    limit=10
)

# Get sample analysis
analysis = manager.get_sample_analysis("SAMPLE_001")
```

## Advanced Features

### Dual-Database Architecture

The VCF Agent uses a sophisticated dual-database approach:

#### LanceDB (Vector Database)
- **Purpose**: Semantic similarity search using AI embeddings
- **Use Cases**: Finding similar variants, AI-powered search
- **Performance**: <100ms query response time

```python
# Vector similarity search
similar_variants = manager.search_similar_variants(
    reference_variant_id="chr1-123456-A-G",
    limit=10
)
```

#### Kuzu (Graph Database)
- **Purpose**: Complex genomic relationships and graph queries
- **Use Cases**: Sample-variant-gene relationships, pathway analysis
- **Performance**: <500ms for complex queries

```python
# Graph relationship queries
sample_variants = manager.get_sample_analysis("SAMPLE_001")
gene_associations = manager.find_variant_genes("chr1-123456-A-G")
```

### AI-Powered Embeddings

The system generates 1536-dimensional embeddings for semantic search:

```python
from vcf_agent.lancedb_integration import VariantEmbeddingService

# Initialize embedding service
embedding_service = VariantEmbeddingService(
    provider="openai",  # or "ollama"
    model="text-embedding-3-large"
)

# Generate embeddings
variant_description = "Pathogenic variant in BRCA1 gene"
embedding = embedding_service.generate_embedding_sync(variant_description)
```

### Batch Processing

For large-scale analysis:

```python
# Batch variant processing
batch_results = manager.batch_add_vcf_variants(
    variants_data=large_variant_list,
    batch_size=1000,
    max_workers=4
)

print(f"Processed {batch_results} variants")
```

### Performance Monitoring

Built-in performance tracking:

```python
# Get performance statistics
stats = manager.get_comprehensive_statistics()

print("Performance Metrics:")
for operation, metrics in stats['performance_stats'].items():
    print(f"- {operation}: {metrics['avg_duration']:.3f}s avg")
```

## Data Management

### Supported File Formats

- **VCF**: Variant Call Format (v4.0-4.3)
- **BCF**: Binary Variant Call Format
- **Compressed**: .vcf.gz, .bcf files with tabix indexes

### Data Ingestion Workflow

```bash
# Step 1: Validate VCF file
vcf-agent validate sample.vcf.gz

# Step 2: Ingest into databases
vcf-agent ingest-vcf sample.vcf.gz \
  --sample-id SAMPLE_001 \
  --sample-name "Patient Sample" \
  --sample-type germline

# Step 3: Verify ingestion
vcf-agent stats --sample-id SAMPLE_001
```

### Data Export Options

```bash
# Export to different formats
vcf-agent export --format json --output results.json
vcf-agent export --format csv --output results.csv
vcf-agent export --format vcf --output filtered.vcf

# Export with filters
vcf-agent export \
  --format json \
  --filter "gene=BRCA1 OR gene=BRCA2" \
  --clinical-significance pathogenic \
  --output brca_pathogenic.json
```

### Database Management

```bash
# Database statistics
vcf-agent db-stats

# Database cleanup
vcf-agent db-cleanup --older-than 30d

# Database backup
vcf-agent db-backup --output backup.tar.gz

# Database restore
vcf-agent db-restore backup.tar.gz
```

## AI-Powered Analysis

### Supported AI Models

#### OpenAI Models
- **GPT-4**: Best for complex analysis and reasoning
- **GPT-3.5-turbo**: Fast and cost-effective
- **Text-embedding-3-large**: For semantic embeddings

#### Cerebras Models
- **Llama-3.1-8B**: Fast inference for basic analysis
- **Llama-3.1-70B**: Advanced reasoning capabilities

#### Local Models (Ollama)
- **Llama3**: Privacy-focused local analysis
- **CodeLlama**: Code and technical analysis
- **Mistral**: Multilingual support

### Analysis Types

#### Variant Interpretation

```bash
vcf-agent ask "Interpret the clinical significance of variants in sample.vcf.gz"
```

#### Comparative Analysis

```bash
vcf-agent ask "Compare variant profiles between sample1.vcf and sample2.vcf"
```

#### Pathway Analysis

```bash
vcf-agent ask "What biological pathways are affected by these variants?"
```

#### Quality Assessment

```bash
vcf-agent ask "Assess the quality and reliability of variants in this file"
```

### Custom Prompts

Create custom analysis prompts:

```python
from vcf_agent.prompt_templates import create_custom_prompt

# Custom analysis prompt
custom_prompt = create_custom_prompt(
    template="""
    Analyze the following variants for {analysis_type}:
    
    Variants: {variant_data}
    
    Focus on:
    1. Clinical significance
    2. Population frequency
    3. Functional impact
    
    Provide a structured summary.
    """,
    variables=["analysis_type", "variant_data"]
)

# Use custom prompt
response = agent.ask_with_prompt(
    prompt=custom_prompt,
    analysis_type="cancer predisposition",
    variant_data=variant_list
)
```

## Performance Optimization

### Hardware Recommendations

#### Minimum Requirements
- **CPU**: 2 cores, 2.0 GHz
- **RAM**: 4GB
- **Storage**: 10GB SSD

#### Recommended Configuration
- **CPU**: 4+ cores, 3.0+ GHz
- **RAM**: 16GB+
- **Storage**: 50GB+ NVMe SSD
- **Network**: High-speed internet for AI API calls

### Performance Tuning

#### Batch Size Optimization

```python
# Optimize batch size based on available memory
optimal_batch_size = min(
    1000,  # Default maximum
    available_memory_gb * 100  # Scale with memory
)

manager = create_data_store_manager(
    batch_size=optimal_batch_size
)
```

#### Worker Thread Configuration

```python
import os

# Set optimal worker count
max_workers = min(
    os.cpu_count(),  # Don't exceed CPU cores
    4  # Reasonable maximum for I/O bound tasks
)

config = DataStoreConfig(max_workers=max_workers)
```

#### Caching Strategies

```python
# Enable embedding caching
embedding_service = VariantEmbeddingService(
    cache_embeddings=True,
    cache_size=10000  # Cache up to 10k embeddings
)
```

### Monitoring Performance

```bash
# Real-time performance monitoring
vcf-agent monitor --real-time

# Performance report
vcf-agent performance-report --output performance.html

# Resource usage
vcf-agent resource-usage --interval 5s
```

## Troubleshooting

### Common Issues

#### 1. API Key Errors

**Problem**: `Authentication failed for OpenAI API`

**Solution**:
```bash
# Verify API key
echo $OPENAI_API_KEY

# Test API connectivity
vcf-agent test-api --provider openai
```

#### 2. Memory Issues

**Problem**: `Out of memory during batch processing`

**Solutions**:
```bash
# Reduce batch size
vcf-agent ingest-vcf sample.vcf.gz --batch-size 500

# Use streaming mode
vcf-agent ingest-vcf sample.vcf.gz --streaming
```

#### 3. Database Connection Issues

**Problem**: `Failed to connect to database`

**Solutions**:
```bash
# Check database paths
vcf-agent db-status

# Reset databases
vcf-agent db-reset --confirm

# Repair databases
vcf-agent db-repair
```

#### 4. Performance Issues

**Problem**: Slow query performance

**Solutions**:
```bash
# Optimize databases
vcf-agent db-optimize

# Check indexes
vcf-agent db-indexes --status

# Rebuild indexes
vcf-agent db-indexes --rebuild
```

### Debug Mode

Enable detailed logging:

```bash
# Enable debug logging
export LOG_LEVEL=DEBUG
vcf-agent analyze sample.vcf.gz --debug

# Trace mode for detailed execution
vcf-agent analyze sample.vcf.gz --trace
```

### Log Analysis

```bash
# View recent logs
vcf-agent logs --tail 100

# Search logs
vcf-agent logs --search "error" --since 1h

# Export logs
vcf-agent logs --export logs.json --since 24h
```

## Best Practices

### Data Organization

#### Directory Structure
```
project/
├── data/
│   ├── raw/           # Original VCF files
│   ├── processed/     # Processed results
│   └── databases/     # LanceDB and Kuzu data
├── config/
│   ├── credentials.json
│   └── vcf_agent_config.yaml
├── results/
│   ├── analyses/      # Analysis results
│   └── exports/       # Exported data
└── logs/              # Application logs
```

#### File Naming Conventions
```bash
# VCF files
{sample_id}_{date}_{type}.vcf.gz
# Example: SAMPLE_001_20250105_germline.vcf.gz

# Analysis results
{sample_id}_{analysis_type}_{date}.json
# Example: SAMPLE_001_pathogenic_analysis_20250105.json
```

### Security Best Practices

#### API Key Management
```bash
# Use environment variables (not hardcoded)
export OPENAI_API_KEY="$(cat ~/.secrets/openai_key)"

# Rotate keys regularly
vcf-agent rotate-keys --provider openai

# Monitor API usage
vcf-agent api-usage --provider openai --since 30d
```

#### Data Protection
```bash
# Encrypt sensitive data
vcf-agent encrypt-data --input sensitive.vcf --output encrypted.vcf.enc

# Anonymize data
vcf-agent anonymize --input patient.vcf --output anonymized.vcf
```

### Workflow Automation

#### Batch Processing Script
```bash
#!/bin/bash
# batch_analysis.sh

for vcf_file in data/raw/*.vcf.gz; do
    sample_id=$(basename "$vcf_file" .vcf.gz)
    
    echo "Processing $sample_id..."
    
    # Validate
    vcf-agent validate "$vcf_file" || continue
    
    # Ingest
    vcf-agent ingest-vcf "$vcf_file" --sample-id "$sample_id"
    
    # Analyze
    vcf-agent ask "Analyze pathogenic variants in $sample_id" \
        --output "results/${sample_id}_analysis.json"
done
```

#### Monitoring Script
```bash
#!/bin/bash
# monitor_performance.sh

while true; do
    vcf-agent performance-report --format json > "logs/perf_$(date +%Y%m%d_%H%M%S).json"
    sleep 300  # Every 5 minutes
done
```

### Quality Assurance

#### Validation Workflow
```bash
# 1. File validation
vcf-agent validate input.vcf.gz --strict

# 2. Data quality check
vcf-agent quality-check input.vcf.gz --report quality_report.html

# 3. Comparison with reference
vcf-agent compare input.vcf.gz reference.vcf.gz --output comparison.json

# 4. Statistical validation
vcf-agent stats input.vcf.gz --detailed --output stats.json
```

#### Testing Pipeline
```bash
# Unit tests
pytest tests/unit/

# Integration tests
pytest tests/integration/

# End-to-end tests
vcf-agent test-e2e --test-data sample_test_data/
```

## Support and Resources

### Documentation
- [API Reference](api.md)
- [Data Stores Architecture](data_stores.md)
- [Security Policy](../SECURITY.md)
- [Contributing Guidelines](../CONTRIBUTING.md)

### Community
- **GitHub Issues**: Report bugs and request features
- **Discussions**: Community Q&A and best practices
- **Wiki**: Community-contributed guides and examples

### Professional Support
- **Enterprise Support**: Priority support and custom features
- **Training**: Workshops and certification programs
- **Consulting**: Custom implementation and optimization

### Getting Help

```bash
# Built-in help
vcf-agent --help
vcf-agent analyze --help

# Version and system info
vcf-agent --version
vcf-agent system-info

# Health check
vcf-agent health-check --comprehensive
```

---

**Next Steps**: After completing this user guide, explore the [API Reference](api.md) for detailed programming interfaces and the [Data Stores Documentation](data_stores.md) for advanced database operations. 