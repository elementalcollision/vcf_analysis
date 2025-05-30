# VCF Analysis Agent Usage Examples

**Last Updated**: May 29, 2025  
**Status**: Production Ready ✅  
**Interface Types**: Natural Language, Direct Tools, CLI, Data Operations

## Overview

The VCF Analysis Agent provides multiple interfaces for genomic analysis, from natural language conversations to direct tool usage and command-line operations. This guide provides comprehensive examples for all interaction methods.

## Natural Language Interface

### Basic Conversations

```python
from src.vcf_agent.agent import get_agent_with_session
from src.vcf_agent.config import SessionConfig

# Initialize agent
config = SessionConfig(raw_mode=False)
agent = get_agent_with_session(config, "ollama")

# Natural conversation
response = agent("Hello! What can you help me with?")
# → "I can help you analyze VCF files, validate formats, compare variants..."

# Get help with specific tasks
response = agent("How do I validate a VCF file?")
# → "I can validate VCF files using the validate_vcf tool..."
```

### Complex Analysis Workflows

```python
# Multi-step analysis request
response = agent("""
Analyze sample_data/patient.vcf:
1. Validate the file format
2. Find pathogenic variants  
3. Generate clinical summary
4. Load into graph database
""")
# → Executes multi-step workflow automatically

# Focused analysis with specific criteria
response = agent("""
Please analyze patient_sample.vcf for:
- High-quality variants (QUAL > 30)
- Variants in BRCA1 and BRCA2 genes
- Clinical significance assessment
- Generate a clinical report
""")

# Comparative analysis
response = agent("""
Compare before_treatment.vcf and after_treatment.vcf:
- Focus on quality differences
- Identify new variants
- Highlight clinical significance changes
""")
```

### Interactive Analysis Sessions

```python
# Start analysis session
agent = get_agent_with_session(config, "ollama")

# Initial request
response1 = agent("Validate sample_data/example.vcf")
print(response1)

# Follow-up analysis
response2 = agent("Now analyze it for pathogenic variants")
print(response2)

# Additional insights
response3 = agent("Show me similar variants in the database")
print(response3)

# Generate report
response4 = agent("Create a clinical summary report")
print(response4)
```

## Direct Tool Usage

### VCF Validation Tools

```python
# Basic validation
result = agent.validate_vcf("sample_data/example.vcf")
print(result)
# → "✅ VCF file 'sample_data/example.vcf' is valid and passed all validation checks."

# Validation with error details
try:
    result = agent.validate_vcf("invalid_file.vcf")
except Exception as e:
    print(f"Validation failed: {e}")

# Batch validation
vcf_files = ["file1.vcf", "file2.vcf", "file3.vcf"]
for vcf_file in vcf_files:
    try:
        result = agent.validate_vcf(vcf_file)
        print(f"✅ {vcf_file}: Valid")
    except Exception as e:
        print(f"❌ {vcf_file}: {e}")
```

### BCFtools Integration

```python
# Filter high-quality variants
result = agent.bcftools_filter_tool(
    input_file="input.vcf",
    output_file="filtered.vcf",
    include_expression="QUAL>30 && DP>10"
)
print(f"Filtering result: {result}")

# Extract specific variant information
result = agent.bcftools_query_tool(
    input_file="input.vcf",
    format_string="%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n"
)
print("Extracted variants:")
print(result)

# Generate comprehensive statistics
stats = agent.bcftools_stats_tool("input.vcf")
print("VCF Statistics:")
print(stats)

# Normalize variants
result = agent.bcftools_norm_tool(
    input_file="input.vcf",
    output_file="normalized.vcf",
    reference_fasta="reference.fa"
)

# View specific regions
result = agent.bcftools_view_tool(
    input_file="input.vcf",
    output_file="subset.vcf",
    regions="chr1:1000000-2000000"
)
```

### AI-Powered Analysis Tools

```python
# Comprehensive variant analysis
analysis = agent.vcf_analysis_summary_tool(
    vcf_file="patient.vcf",
    analysis_type="clinical"
)
print("Clinical Analysis:")
print(analysis)

# Compare two VCF files
comparison = agent.ai_vcf_comparison_tool(
    vcf_file1="before.vcf",
    vcf_file2="after.vcf",
    focus="quality_differences"
)
print("File Comparison:")
print(comparison)

# Focused pathogenicity analysis
pathogenicity = agent.vcf_analysis_summary_tool(
    vcf_file="variants.vcf",
    analysis_type="pathogenicity"
)
print("Pathogenicity Assessment:")
print(pathogenicity)

# Quality assessment
quality = agent.vcf_analysis_summary_tool(
    vcf_file="raw_variants.vcf", 
    analysis_type="quality"
)
print("Quality Assessment:")
print(quality)
```

### Graph Database Operations

```python
# Load VCF into graph database
result = agent.load_vcf_into_graph_db_tool(
    vcf_file="patient.vcf",
    sample_id="PATIENT_001"
)
print(f"Graph loading result: {result}")

# Load multiple samples
samples = [
    ("patient1.vcf", "PATIENT_001"),
    ("patient2.vcf", "PATIENT_002"),
    ("patient3.vcf", "PATIENT_003")
]

for vcf_file, sample_id in samples:
    result = agent.load_vcf_into_graph_db_tool(
        vcf_file=vcf_file,
        sample_id=sample_id
    )
    print(f"Loaded {sample_id}: {result}")
```

## Command Line Interface

### Quick Analysis Commands

```bash
# Basic analysis
vcf-agent analyze sample_data/example.vcf --output results/

# Validation only
vcf-agent validate sample_data/example.vcf

# Statistics generation
vcf-agent stats sample_data/example.vcf --format json

# Help and available commands
vcf-agent --help
vcf-agent analyze --help
```

### Comprehensive Workflows

```bash
# Full analysis workflow
vcf-agent workflow \
  --input patient.vcf \
  --validate \
  --ai-analysis \
  --graph-load \
  --output clinical_report.json

# Quality control pipeline
vcf-agent pipeline \
  --input raw_variants.vcf \
  --filter "QUAL>30 && DP>10" \
  --normalize \
  --annotate \
  --output high_quality_variants.vcf

# Comparative analysis
vcf-agent compare \
  --file1 before_treatment.vcf \
  --file2 after_treatment.vcf \
  --focus clinical_significance \
  --output comparison_report.html
```

### Batch Processing

```bash
# Process multiple files
vcf-agent batch process_list.txt --parallel 4

# Example process_list.txt:
# patient1.vcf,SAMPLE_001
# patient2.vcf,SAMPLE_002  
# patient3.vcf,SAMPLE_003

# Batch validation
find /data/vcf_files -name "*.vcf" | vcf-agent batch-validate --parallel 8

# Batch statistics
vcf-agent batch-stats \
  --input-dir /data/vcf_files \
  --output-dir /results/stats \
  --format json
```

### Search and Query Operations

```bash
# Search similar variants
vcf-agent search "pathogenic BRCA1 variant" --limit 10

# Search with filters
vcf-agent search \
  --query "high quality variant" \
  --filter "quality > 30" \
  --gene BRCA1,BRCA2 \
  --limit 20

# Graph queries
vcf-agent graph-query \
  --sample PATIENT_001 \
  --relationship "similar_to" \
  --limit 15
```

## Data Store Operations

### Database Initialization

```python
from src.vcf_agent.data_store_manager import create_data_store_manager

# Initialize data store manager
manager = create_data_store_manager(
    lancedb_path="./data/lancedb",
    kuzu_path="./data/kuzu_db"
)

# Verify database connections
try:
    status = manager.get_status()
    print(f"Database status: {status}")
except Exception as e:
    print(f"Database connection failed: {e}")
```

### Sample and Variant Management

```python
# Add single sample with variants
sample_data = {"id": "SAMPLE_001", "name": "Patient 1"}
variants_data = [{
    "id": "chr1-123456-A-G",
    "chr": "1", 
    "pos": 123456,
    "ref": "A", 
    "alt": "G",
    "clinical_significance": "Pathogenic"
}]

result = manager.add_sample_with_variants(
    sample_data=sample_data,
    variants_data=variants_data
)
print(f"Sample addition result: {result}")

# Batch add multiple samples
samples_batch = [
    {
        "sample_data": {"id": "SAMPLE_002", "name": "Patient 2"},
        "variants_data": [
            {"id": "chr2-234567-C-T", "chr": "2", "pos": 234567, "ref": "C", "alt": "T"}
        ]
    },
    {
        "sample_data": {"id": "SAMPLE_003", "name": "Patient 3"}, 
        "variants_data": [
            {"id": "chr3-345678-G-A", "chr": "3", "pos": 345678, "ref": "G", "alt": "A"}
        ]
    }
]

for batch_item in samples_batch:
    result = manager.add_sample_with_variants(**batch_item)
    print(f"Batch result: {result}")
```

### Search and Analysis Operations

```python
# Search for similar variants
results = manager.search_variants(
    query="pathogenic variant in BRCA1",
    search_type="hybrid",
    limit=10
)
print("Similar variants found:")
for result in results:
    print(f"- {result}")

# Metadata filtering
filtered_results = manager.search_variants(
    query="high quality variant",
    search_type="metadata",
    filters={"quality_score": ">30", "chromosome": "1"},
    limit=15
)

# Vector similarity search
vector_results = manager.search_variants(
    query="missense mutation pathogenic",
    search_type="vector",
    limit=5
)

# Get sample analysis
analysis = manager.get_sample_analysis("SAMPLE_001")
print("Sample analysis:")
print(analysis)
```

### Advanced Data Operations

```python
# Get database statistics
stats = manager.get_database_stats()
print("Database statistics:")
print(f"Total samples: {stats.get('total_samples', 0)}")
print(f"Total variants: {stats.get('total_variants', 0)}")
print(f"Database size: {stats.get('database_size_mb', 0)} MB")

# Export data
export_result = manager.export_sample_data(
    sample_id="SAMPLE_001",
    output_format="json",
    output_file="sample_001_export.json"
)

# Import data
import_result = manager.import_sample_data(
    input_file="sample_data.json",
    data_format="json"
)

# Cleanup operations
cleanup_result = manager.cleanup_old_data(days_old=30)
print(f"Cleanup result: {cleanup_result}")
```

## Error Handling Patterns

### Graceful Error Handling

```python
def robust_vcf_analysis(vcf_file):
    try:
        # Primary analysis path
        agent = get_agent_with_session(config, "ollama")
        
        # Validate first
        validation_result = agent.validate_vcf(vcf_file)
        print(f"Validation: {validation_result}")
        
        # Perform analysis
        analysis_result = agent.vcf_analysis_summary_tool(vcf_file)
        return analysis_result
        
    except FileNotFoundError:
        print(f"Error: VCF file '{vcf_file}' not found")
        return None
        
    except ValidationError as e:
        print(f"VCF validation failed: {e}")
        return None
        
    except Exception as e:
        print(f"Analysis failed: {e}")
        # Fallback to basic stats
        try:
            return agent.bcftools_stats_tool(vcf_file)
        except:
            return None

# Usage
result = robust_vcf_analysis("patient.vcf")
if result:
    print("Analysis successful:")
    print(result)
else:
    print("Analysis failed - check file and try again")
```

### Retry Mechanisms

```python
import time
from typing import Optional

def retry_analysis(vcf_file: str, max_retries: int = 3) -> Optional[str]:
    for attempt in range(max_retries):
        try:
            agent = get_agent_with_session(config, "ollama")
            result = agent.vcf_analysis_summary_tool(vcf_file)
            return result
            
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt  # Exponential backoff
                print(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print("All retry attempts exhausted")
                return None

# Usage
result = retry_analysis("large_file.vcf")
```

## Performance Optimization Tips

### Batch Processing Optimization

```python
# Efficient batch processing
def process_vcf_batch(vcf_files: list, batch_size: int = 10):
    agent = get_agent_with_session(config, "ollama")
    results = []
    
    for i in range(0, len(vcf_files), batch_size):
        batch = vcf_files[i:i + batch_size]
        batch_results = []
        
        for vcf_file in batch:
            try:
                result = agent.validate_vcf(vcf_file)
                batch_results.append({"file": vcf_file, "status": "valid"})
            except Exception as e:
                batch_results.append({"file": vcf_file, "status": "invalid", "error": str(e)})
        
        results.extend(batch_results)
        print(f"Processed batch {i//batch_size + 1}/{(len(vcf_files)-1)//batch_size + 1}")
    
    return results

# Usage
vcf_files = ["file1.vcf", "file2.vcf", "file3.vcf"]  # ... many files
results = process_vcf_batch(vcf_files, batch_size=5)
```

### Memory-Conscious Operations

```python
# Memory-optimized configuration
from src.vcf_agent.config import SessionConfig, MemoryOptimizationConfig

memory_config = MemoryOptimizationConfig(
    optimization_level="aggressive",
    memory_cleanup_threshold_mb=100,
    streaming_batch_size=25
)

config = SessionConfig(
    raw_mode=False,
    memory_optimization=memory_config
)

agent = get_agent_with_session(config, "ollama")

# Use streaming for large files
def analyze_large_vcf(vcf_file: str):
    try:
        # Enable streaming mode for large files
        result = agent.vcf_analysis_summary_tool(
            vcf_file=vcf_file,
            analysis_type="streaming"
        )
        return result
    except MemoryError:
        print("Memory limit reached - consider splitting the file")
        return None
```

## Integration Examples

### Workflow Integration

```python
def genomic_analysis_pipeline(vcf_file: str, sample_id: str):
    """Complete genomic analysis pipeline"""
    
    agent = get_agent_with_session(config, "ollama")
    manager = create_data_store_manager()
    
    pipeline_results = {}
    
    # Step 1: Validation
    print("Step 1: Validating VCF file...")
    try:
        validation = agent.validate_vcf(vcf_file)
        pipeline_results["validation"] = validation
        print("✅ Validation passed")
    except Exception as e:
        print(f"❌ Validation failed: {e}")
        return pipeline_results
    
    # Step 2: Quality filtering
    print("Step 2: Quality filtering...")
    try:
        filtered_file = f"filtered_{vcf_file}"
        filter_result = agent.bcftools_filter_tool(
            input_file=vcf_file,
            output_file=filtered_file,
            include_expression="QUAL>30 && DP>10"
        )
        pipeline_results["filtering"] = filter_result
        print("✅ Filtering completed")
    except Exception as e:
        print(f"❌ Filtering failed: {e}")
        filtered_file = vcf_file  # Use original file
    
    # Step 3: AI analysis
    print("Step 3: AI analysis...")
    try:
        analysis = agent.vcf_analysis_summary_tool(filtered_file)
        pipeline_results["analysis"] = analysis
        print("✅ Analysis completed")
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
    
    # Step 4: Database loading
    print("Step 4: Loading into database...")
    try:
        db_result = agent.load_vcf_into_graph_db_tool(filtered_file, sample_id)
        pipeline_results["database"] = db_result
        print("✅ Database loading completed")
    except Exception as e:
        print(f"❌ Database loading failed: {e}")
    
    return pipeline_results

# Usage
results = genomic_analysis_pipeline("patient.vcf", "PATIENT_001")
print("Pipeline completed:")
for step, result in results.items():
    print(f"{step}: {result}")
```

## Links and References

- [Architecture Guide](ARCHITECTURE_GUIDE.md)
- [Memory Optimization Features](MEMORY_OPTIMIZATION_FEATURES.md)
- [Production Monitoring](PRODUCTION_MONITORING.md)
- [Available Tools Guide](TOOLS_GUIDE.md)
- [Developer Guide](DEVELOPER_GUIDE.md)

---

**Next Steps**: Explore specific tool documentation for advanced usage patterns and optimization techniques. 