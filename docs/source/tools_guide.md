# VCF Agent Tools Guide

## Overview

The VCF Agent provides a comprehensive suite of AI-powered tools for genomic data analysis. This guide covers all available tools, usage patterns, and best practices for effective VCF file analysis.

## 🚀 Quick Start

### Basic Setup

```python
from src.vcf_agent.agent import get_agent_with_session
from src.vcf_agent.config import SessionConfig

# Create agent with natural conversation mode
config = SessionConfig(raw_mode=False)
agent = get_agent_with_session(config, "ollama")
```

### Simple Tool Usage

```python
# Natural language - agent automatically selects and executes tools
response = agent("Please validate the VCF file at sample_data/example.vcf")

# Direct tool calling
result = agent.validate_vcf("sample_data/example.vcf")
```

## 🛠️ Available Tools

### Core Validation Tools

#### echo
**Purpose**: Test agent connectivity and tool execution
**Usage**: Simple echo functionality for debugging

```python
# Direct call
result = agent.echo("Hello World!")
# Returns: "Echo: Hello World!"

# Natural language
response = agent("Echo the message 'System is working'")
```

#### validate_vcf
**Purpose**: Comprehensive VCF/BCF file validation
**Features**: 
- File existence check
- Format validation
- Index validation
- BCFtools compatibility

```python
# Direct call
result = agent.validate_vcf("sample_data/example.vcf")
# Returns: "✅ VCF file 'sample_data/example.vcf' is valid and passed all validation checks."

# Natural language
response = agent("Please validate the VCF file at sample_data/example.vcf")
response = agent("Check if my VCF file is properly formatted")
```

### BCFtools Integration Tools

The VCF Agent provides complete BCFtools integration with all major commands:

#### bcftools_view_tool
**Purpose**: View and subset VCF/BCF files
**Features**: Region filtering, sample selection, format conversion

```python
# Basic usage
result = agent.bcftools_view_tool("input.vcf", "output.vcf")

# Advanced filtering
result = agent.bcftools_view_tool(
    input_file="input.vcf",
    output_file="filtered.vcf",
    regions="chr1:1000000-2000000",
    samples="SAMPLE1,SAMPLE2",
    include_filter="QUAL>30"
)

# Natural language
response = agent("Extract variants from chromosome 1 between positions 1M and 2M")
response = agent("Show me only high-quality variants with QUAL > 30")
```

#### bcftools_query_tool
**Purpose**: Extract specific data from VCF files
**Features**: Custom format strings, flexible data extraction

```python
# Extract basic variant information
result = agent.bcftools_query_tool(
    input_file="input.vcf",
    format_string="%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n"
)

# Extract genotype information
result = agent.bcftools_query_tool(
    input_file="input.vcf",
    format_string="%CHROM\t%POS\t[%GT\t%DP]\n",
    samples="SAMPLE1"
)

# Natural language
response = agent("Extract chromosome, position, and quality scores from my VCF file")
response = agent("Get genotype information for sample SAMPLE1")
```

#### bcftools_filter_tool
**Purpose**: Filter variants based on criteria
**Features**: Include/exclude expressions, soft filtering

```python
# Quality filtering
result = agent.bcftools_filter_tool(
    input_file="input.vcf",
    output_file="filtered.vcf",
    include_expression="QUAL>30 && DP>10"
)

# Soft filtering (mark but don't remove)
result = agent.bcftools_filter_tool(
    input_file="input.vcf",
    output_file="filtered.vcf",
    soft_filter="LowQual",
    exclude_expression="QUAL<20"
)

# Natural language
response = agent("Filter variants to keep only high-quality ones with QUAL > 30")
response = agent("Remove variants with low depth coverage")
```

#### bcftools_norm_tool
**Purpose**: Normalize variants (left-align, split multiallelic)
**Features**: Reference-based normalization, multiallelic handling

```python
# Basic normalization
result = agent.bcftools_norm_tool(
    input_file="input.vcf",
    output_file="normalized.vcf",
    reference_fasta="reference.fa"
)

# Split multiallelic variants
result = agent.bcftools_norm_tool(
    input_file="input.vcf",
    output_file="normalized.vcf",
    multiallelics="-both"
)

# Natural language
response = agent("Normalize my VCF file and split multiallelic variants")
response = agent("Left-align indels using the reference genome")
```

#### bcftools_stats_tool
**Purpose**: Generate comprehensive VCF statistics
**Features**: Detailed metrics, sample-specific stats, quality distributions

```python
# Basic statistics
result = agent.bcftools_stats_tool("input.vcf")

# Statistics with reference
result = agent.bcftools_stats_tool(
    input_file="input.vcf",
    fasta_ref="reference.fa",
    output_file="stats.txt"
)

# Natural language
response = agent("Generate comprehensive statistics for my VCF file")
response = agent("Show me quality distributions and variant counts")
```

#### bcftools_annotate_tool
**Purpose**: Annotate variants with additional information
**Features**: Custom annotations, database integration

```python
# Annotate with database
result = agent.bcftools_annotate_tool(
    input_file="input.vcf",
    output_file="annotated.vcf",
    annotations_file="annotations.bed",
    columns="CHROM,FROM,TO,GENE"
)

# Natural language
response = agent("Annotate my variants with gene information")
response = agent("Add clinical significance annotations to my VCF file")
```

### AI-Powered Analysis Tools

#### ai_vcf_comparison_tool
**Purpose**: Intelligent comparison of VCF files with AI insights
**Features**: Multi-level comparison, clinical interpretation, statistical analysis

```python
# Basic comparison
result = agent.ai_vcf_comparison_tool("file1.vcf", "file2.vcf")

# Comprehensive comparison
result = agent.ai_vcf_comparison_tool(
    file1="patient1.vcf",
    file2="patient2.vcf",
    comparison_type="comprehensive"
)

# Clinical comparison
result = agent.ai_vcf_comparison_tool(
    file1="before_treatment.vcf",
    file2="after_treatment.vcf",
    comparison_type="clinical"
)

# Natural language
response = agent("Compare these two VCF files and highlight important differences")
response = agent("Analyze the differences between patient samples")
```

#### vcf_analysis_summary_tool
**Purpose**: AI-powered comprehensive VCF analysis
**Features**: Variant interpretation, clinical insights, quality assessment

```python
# Comprehensive analysis
result = agent.vcf_analysis_summary_tool("patient.vcf", "comprehensive")

# Clinical analysis
result = agent.vcf_analysis_summary_tool("patient.vcf", "clinical")

# Basic analysis
result = agent.vcf_analysis_summary_tool("patient.vcf", "basic")

# Natural language
response = agent("Provide a comprehensive analysis of my VCF file")
response = agent("Analyze this patient's variants for clinical significance")
```

#### vcf_summarization_tool
**Purpose**: Enhanced VCF summarization with LLM fallback
**Features**: Multi-model support, detailed summaries, error resilience

```python
# Comprehensive summary
result = agent.vcf_summarization_tool("patient.vcf", "comprehensive")

# Natural language
response = agent("Summarize the key findings in my VCF file")
response = agent("Give me an overview of the variants in this sample")
```

### Database Integration Tools

#### load_vcf_into_graph_db_tool
**Purpose**: Load VCF data into Kuzu graph database
**Features**: Batch processing, relationship modeling, performance optimization

```python
# Load VCF into graph database
result = agent.load_vcf_into_graph_db_tool(
    vcf_file="patient.vcf",
    sample_name="Patient_001",
    batch_size=1000
)

# Natural language
response = agent("Load my VCF file into the graph database as Patient_001")
response = agent("Import variants into the knowledge graph")
```

## 🔄 Usage Patterns

### 1. Natural Language Workflows

The agent automatically selects and executes appropriate tools based on natural language requests:

```python
# Single tool execution
response = agent("Validate my VCF file at sample_data/example.vcf")

# Multi-tool workflows
response = agent("""
Please perform a complete analysis of sample_data/patient.vcf:
1. First validate the file
2. Generate comprehensive statistics
3. Create a clinical analysis summary
4. Load the data into the graph database as 'Patient_123'
""")

# Comparative analysis
response = agent("Compare patient_before.vcf and patient_after.vcf to identify treatment effects")
```

### 2. Direct Tool Calling

For programmatic usage and precise control:

```python
# Sequential tool calls
validation = agent.validate_vcf("patient.vcf")
stats = agent.bcftools_stats_tool("patient.vcf")
analysis = agent.vcf_analysis_summary_tool("patient.vcf", "clinical")

# Conditional execution
if "valid" in validation.lower():
    summary = agent.vcf_summarization_tool("patient.vcf")
    graph_load = agent.load_vcf_into_graph_db_tool("patient.vcf", "Patient_001")
```

### 3. Tool Chaining

Combine multiple tools for complex workflows:

```python
# Quality control pipeline
response = agent("""
Quality control pipeline for sample_data/raw.vcf:
1. Validate the input file
2. Filter variants with QUAL > 30 and DP > 10
3. Normalize and left-align variants
4. Generate final statistics
5. Create analysis summary
""")

# Comparative genomics workflow
response = agent("""
Compare two patient samples:
1. Validate both files: patient1.vcf and patient2.vcf
2. Generate statistics for each
3. Perform AI-powered comparison
4. Summarize key differences
""")
```

## 📊 Performance and Monitoring

### Automatic Metrics Collection

All tools automatically collect performance metrics:

```python
# Metrics are automatically recorded for:
# - Execution time
# - Success/failure status
# - Error types and messages
# - Tool usage patterns

# Access metrics programmatically
from src.vcf_agent import metrics
# Metrics are available through the metrics module
```

### Error Handling

Tools provide comprehensive error handling:

```python
# Tools return user-friendly error messages
result = agent.validate_vcf("nonexistent.vcf")
# Returns: "❌ Error: File 'nonexistent.vcf' not found. Please check the file path."

# Natural language error handling
response = agent("Validate the file that_doesnt_exist.vcf")
# Agent provides helpful error message and suggestions
```

## 🎯 Best Practices

### 1. File Path Management

```python
# Use absolute paths for reliability
result = agent.validate_vcf("/full/path/to/file.vcf")

# Use relative paths from project root
result = agent.validate_vcf("sample_data/example.vcf")

# Check file existence first
response = agent("Check if sample_data/example.vcf exists and validate it")
```

### 2. Workflow Organization

```python
# Break complex workflows into steps
response = agent("""
Step 1: Validate input.vcf
Step 2: If valid, generate statistics
Step 3: Create analysis summary
Step 4: Load into database if analysis is successful
""")
```

### 3. Error Recovery

```python
# Use natural language for error recovery
response = agent("""
Try to validate my_file.vcf. If it fails, suggest what might be wrong
and how to fix it.
""")
```

### 4. Performance Optimization

```python
# Use appropriate batch sizes for large files
result = agent.load_vcf_into_graph_db_tool(
    "large_file.vcf", 
    "Sample_001", 
    batch_size=5000  # Larger batch for better performance
)

# Use specific analysis types for focused results
result = agent.vcf_analysis_summary_tool("file.vcf", "basic")  # Faster
result = agent.vcf_analysis_summary_tool("file.vcf", "comprehensive")  # More detailed
```

## 🔧 Configuration

### Agent Configuration

```python
# Configure for different use cases
config = SessionConfig(
    raw_mode=False,  # Enable chain-of-thought reasoning
    max_tokens=4000,  # Adjust for longer responses
    temperature=0.1   # Lower for more deterministic results
)

agent = get_agent_with_session(config, "ollama")
```

### Model Selection

```python
# Choose appropriate model for your needs
agent_ollama = get_agent_with_session(config, "ollama")      # Local, fast
agent_openai = get_agent_with_session(config, "openai")     # Cloud, powerful
agent_cerebras = get_agent_with_session(config, "cerebras") # Fast inference
```

## 🚨 Troubleshooting

### Common Issues

1. **Tool not executing automatically**
   ```python
   # Ensure raw_mode is False
   config = SessionConfig(raw_mode=False)
   ```

2. **File not found errors**
   ```python
   # Use absolute paths or verify working directory
   import os
   print(f"Current directory: {os.getcwd()}")
   ```

3. **BCFtools errors**
   ```python
   # Ensure BCFtools is installed and in PATH
   response = agent("Check if bcftools is available")
   ```

### Getting Help

```python
# Ask the agent for help
response = agent("What tools are available for VCF analysis?")
response = agent("How do I validate a VCF file?")
response = agent("What's the best way to compare two VCF files?")
```

## 📚 Examples

### Complete Analysis Pipeline

```python
# Comprehensive VCF analysis workflow
response = agent("""
Please perform a complete analysis of sample_data/patient.vcf:

1. Validate the VCF file format and structure
2. Generate comprehensive statistics including:
   - Variant counts by type
   - Quality score distributions
   - Sample-specific metrics
3. Create a clinical analysis summary focusing on:
   - Pathogenic variants
   - Gene associations
   - Clinical significance
4. Load the validated data into the graph database as 'Patient_001'
5. Provide a final summary with key findings and recommendations

Please provide detailed output for each step.
""")
```

### Quality Control Pipeline

```python
# Quality control and filtering workflow
response = agent("""
Quality control pipeline for sample_data/raw_variants.vcf:

1. Validate the input file
2. Generate initial statistics
3. Filter variants to keep only:
   - QUAL > 30
   - DP > 10
   - No missing genotypes
4. Normalize variants (left-align, split multiallelic)
5. Generate final statistics
6. Compare before/after filtering results
7. Create final analysis summary

Save filtered results to sample_data/qc_variants.vcf
""")
```

This comprehensive tools guide provides everything needed to effectively use the VCF Agent's tools for genomic data analysis. 