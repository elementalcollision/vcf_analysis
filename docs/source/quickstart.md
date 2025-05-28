# Quick Start Guide

Get up and running with the VCF Analysis Agent in under 10 minutes! This guide walks you through your first VCF analysis using AI-powered insights.

## Prerequisites

Before starting, ensure you have:
- ✅ Completed the [Installation Guide](installation.md)
- ✅ A VCF file to analyze (we provide sample data)
- ✅ At least one AI provider configured (OpenAI, Cerebras, or Ollama)

## Step 1: Verify Installation

First, let's make sure everything is working:

```bash
# Check the agent is installed
vcf-agent --version

# Expected output:
# VCF Analysis Agent v1.0.0
```

## Step 2: Download Sample Data

We'll use a sample VCF file for this quickstart:

```bash
# Download sample VCF file (if not already present)
curl -O https://raw.githubusercontent.com/your-org/vcf-analysis-agent/main/sample_data/example.vcf

# Or use the included sample data
ls sample_data/
```

## Step 3: Basic VCF Validation

Let's start with basic validation to ensure your VCF file is properly formatted:

```bash
# Validate the VCF file
vcf-agent validate sample_data/example.vcf
```

Expected output:
```
✅ VCF Validation Results
📁 File: sample_data/example.vcf
📊 Format: Valid VCF 4.2
🧬 Variants: 1,234 variants found
👥 Samples: 3 samples detected
✅ Status: PASSED - File is valid and ready for analysis
```

## Step 4: Your First AI Analysis

Now let's perform AI-powered analysis of the VCF file:

```bash
# Analyze with Ollama (local, no API key needed)
vcf-agent analyze sample_data/example.vcf --provider ollama

# Or with OpenAI (requires API key)
vcf-agent analyze sample_data/example.vcf --provider openai

# Or with Cerebras (requires API key)
vcf-agent analyze sample_data/example.vcf --provider cerebras
```

Expected output:
```
🧬 VCF Analysis Results
📁 File: sample_data/example.vcf
🤖 AI Provider: ollama (qwen3:4b)
⏱️  Analysis Time: 12.3 seconds

📊 Summary:
• Total Variants: 1,234
• Pathogenic: 23 variants
• Likely Pathogenic: 45 variants
• Benign: 892 variants
• Uncertain Significance: 274 variants

🎯 Key Findings:
• High-impact variants in BRCA1, TP53 genes
• 3 novel variants requiring further investigation
• Overall variant quality score: 8.7/10

💡 AI Insights:
The analysis reveals several clinically significant variants...
[Detailed AI-generated analysis follows]
```

## Step 5: Explore Advanced Features

### Vector Similarity Search

Find variants similar to a specific variant:

```bash
# Search for variants similar to a specific variant
vcf-agent search "chr1:123456:A>G" --top-k 10
```

### Graph Database Queries

Query relationships between variants and samples:

```bash
# Find all variants in a specific gene
vcf-agent query "MATCH (v:Variant) WHERE v.gene = 'BRCA1' RETURN v"

# Find variants shared between samples
vcf-agent query "MATCH (v:Variant)-[:ObservedIn]->(s:Sample) RETURN v, count(s) as sample_count ORDER BY sample_count DESC"
```

### Batch Processing

Process multiple VCF files at once:

```bash
# Analyze multiple files
vcf-agent batch-analyze sample_data/*.vcf --output results/
```

## Step 6: Interactive Analysis

Start an interactive session for exploratory analysis:

```bash
# Start interactive mode
vcf-agent interactive sample_data/example.vcf
```

In interactive mode, you can ask natural language questions:
```
> What are the most significant variants in this file?
> Show me all variants in the BRCA1 gene
> Compare variant quality scores across samples
> Find variants with uncertain clinical significance
```

## Step 7: Generate Reports

Create comprehensive analysis reports:

```bash
# Generate HTML report
vcf-agent report sample_data/example.vcf --format html --output analysis_report.html

# Generate PDF report
vcf-agent report sample_data/example.vcf --format pdf --output analysis_report.pdf

# Generate JSON for programmatic use
vcf-agent report sample_data/example.vcf --format json --output analysis_report.json
```

## Step 8: Monitor Performance

Check system performance and metrics:

```bash
# View performance metrics
vcf-agent metrics

# Check database status
vcf-agent status
```

Access the monitoring dashboard at: http://localhost:8000/metrics

## Common Use Cases

### Clinical Variant Analysis

```bash
# Focus on clinically significant variants
vcf-agent analyze patient.vcf --filter "clinical_significance=pathogenic,likely_pathogenic"
```

### Research Variant Discovery

```bash
# Find novel or rare variants
vcf-agent analyze research.vcf --filter "frequency<0.01" --include-novel
```

### Quality Control

```bash
# Comprehensive QC analysis
vcf-agent qc sample.vcf --include-metrics --generate-plots
```

### Comparative Analysis

```bash
# Compare two VCF files
vcf-agent compare file1.vcf file2.vcf --reference hg38.fa
```

## Next Steps

Now that you've completed the quickstart:

1. **Deep Dive**: Read the complete [User Guide](user_guide.md) for advanced features
2. **Configuration**: Customize settings in the [Configuration Guide](configuration.md)
3. **API Integration**: Explore programmatic access via the [API Documentation](api.md)
4. **Production Deployment**: Set up production systems with the [Deployment Guide](deployment.md)

## Troubleshooting

### Common Issues

**Analysis is slow**
```bash
# Use smaller batch sizes
vcf-agent analyze large.vcf --batch-size 1000

# Use local Ollama instead of API providers
vcf-agent analyze large.vcf --provider ollama
```

**Out of memory errors**
```bash
# Enable streaming mode for large files
vcf-agent analyze huge.vcf --streaming --memory-limit 8GB
```

**AI provider errors**
```bash
# Check API key configuration
echo $OPENAI_API_KEY

# Test with different provider
vcf-agent analyze sample.vcf --provider ollama
```

### Getting Help

- **Documentation**: Complete guides in `docs/`
- **Examples**: Sample code in `examples/`
- **Issues**: Report problems on GitHub
- **Community**: Join discussions and get support

## Performance Tips

1. **Use Docker**: Containerized deployment is optimized and recommended
2. **Local AI Models**: Ollama provides fast, private analysis without API costs
3. **Batch Processing**: Process multiple files together for efficiency
4. **Caching**: Enable caching to speed up repeated analyses
5. **Resource Allocation**: Allocate sufficient RAM for large VCF files

---

**Congratulations!** 🎉

You've successfully completed your first VCF analysis with AI-powered insights. The VCF Analysis Agent is now ready for your genomics research and clinical workflows.

**What's Next?**
- Explore advanced features in the [User Guide](user_guide.md)
- Set up production monitoring with [Prometheus](monitoring_with_prometheus.md)
- Integrate with your existing workflows using the [API](api.md) 