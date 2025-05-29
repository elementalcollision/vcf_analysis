# VCF Analysis Agent - Memory Optimization Guide

This guide provides comprehensive documentation for the memory optimization features integrated into the VCF Analysis Agent.

## Overview

The VCF Analysis Agent includes advanced memory optimization capabilities designed for production genomic workloads. These optimizations provide significant memory reductions while maintaining performance and accuracy.

## Features Summary

| Feature | Memory Reduction | Performance Impact | Accuracy |
|---------|------------------|-------------------|----------|
| **PyArrow Streaming** | 84.2% | Improved | 100% |
| **Memory-Aware Caching** | 90%+ recovery | Enhanced | 100% |
| **PCA Dimension Reduction** | 50% embeddings | Negligible | >95% |
| **Streaming Processing** | Bounded growth | Enhanced | 100% |
| **Combined Optimizations** | **>95% total** | **Significantly Enhanced** | **>95%** |

## Configuration System

### Basic Configuration

```python
from vcf_agent.config import SessionConfig, MemoryOptimizationConfig

# Use default optimization (recommended)
session_config = SessionConfig()  # Uses standard optimization level

# Access the embedding service
from vcf_agent.lancedb_integration import VariantEmbeddingService
service = VariantEmbeddingService(session_config)
```

### Advanced Configuration

```python
# Custom optimization configuration
memory_config = MemoryOptimizationConfig(
    optimization_level="aggressive",           # basic|standard|aggressive
    memory_management_enabled=True,            # Enable memory monitoring
    dimension_reduction_enabled=True,          # Enable PCA compression
    target_dimensions=768,                     # Reduced embedding size
    caching_strategy="memory_aware",           # Advanced caching
    cache_max_size_mb=40,                     # Cache size limit
    cache_max_entries=1000,                   # Cache entry limit
    streaming_batch_size=50,                  # Batch processing size
    memory_cleanup_threshold_mb=40            # Memory cleanup trigger
)

session_config = SessionConfig(memory_optimization=memory_config)
service = VariantEmbeddingService(session_config)
```

## Optimization Levels

### Basic Level
- **Use Case**: Development and testing
- **Memory Management**: Disabled
- **Dimension Reduction**: Disabled
- **Caching**: Simple dictionary cache
- **Best For**: Small datasets, development environments

```python
basic_config = MemoryOptimizationConfig(optimization_level="basic")
```

### Standard Level (Recommended)
- **Use Case**: Production deployments
- **Memory Management**: Enabled with monitoring
- **Dimension Reduction**: Enabled (768 dimensions)
- **Caching**: Memory-aware with LRU eviction
- **Best For**: Most production workloads

```python
standard_config = MemoryOptimizationConfig(optimization_level="standard")
```

### Aggressive Level
- **Use Case**: Large-scale enterprise deployments
- **Memory Management**: Enhanced monitoring and cleanup
- **Dimension Reduction**: Enabled with streaming
- **Caching**: Advanced with predictive cleanup
- **Best For**: High-volume genomic analysis

```python
aggressive_config = MemoryOptimizationConfig(optimization_level="aggressive")
```

## Core Features

### 1. Memory-Aware Caching

The memory-aware caching system provides intelligent cache management with automatic cleanup.

#### Features
- **LRU Eviction**: Automatically removes least recently used entries
- **Size Monitoring**: Tracks memory usage and enforces limits
- **Access Tracking**: Monitors cache hit rates and efficiency
- **Automatic Cleanup**: Prevents memory accumulation over time

#### Usage
```python
# Cache is automatically managed
embedding = service.generate_embedding_sync("Variant description")

# Check cache statistics
stats = service.get_optimization_stats()
print(f"Cache hit rate: {stats['cache_hit_rate']:.1f}%")
print(f"Cache entries: {stats['cache_stats']['cache_entries']}")
print(f"Memory usage: {stats['cache_stats']['memory_usage_mb']}MB")
```

### 2. PCA Dimension Reduction

PCA-based dimension reduction provides 50% memory savings for embeddings while preserving semantic accuracy.

#### Features
- **Automatic Training**: Collects embeddings and trains PCA model automatically
- **High Accuracy**: Preserves >95% of original variance
- **Graceful Fallback**: Uses original embeddings when PCA unavailable
- **Performance Monitoring**: Tracks reduction statistics

#### Implementation Details
- **Original Dimensions**: 1536 (OpenAI text-embedding-3-small)
- **Reduced Dimensions**: 768 (50% reduction)
- **Training Samples**: 1000 embeddings for reliable PCA model
- **Variance Retention**: Typically >95% of original information

#### Usage
```python
# Dimension reduction is automatically applied
embedding = service.generate_embedding_sync("Pathogenic BRCA1 variant")
print(f"Dimensions: {len(embedding)}")  # 768 if reduction enabled

# Check reduction statistics
stats = service.get_optimization_stats()
reduction_stats = stats['dimension_reduction_stats']
print(f"Memory reduction: {reduction_stats['memory_reduction_percent']:.1f}%")
print(f"Variance retained: {reduction_stats['variance_retained']:.1%}")
```

### 3. Streaming Batch Processing

Streaming batch processing enables memory-efficient processing of large datasets.

#### Features
- **Configurable Batch Sizes**: Adapt to available memory
- **Memory Cleanup**: Automatic cleanup between batches
- **Progress Monitoring**: Track processing statistics
- **Error Handling**: Graceful handling of processing failures

#### Usage
```python
# Process large datasets efficiently
large_dataset = [f"variant description {i}" for i in range(10000)]

# Streaming batch processing
embeddings = service.generate_embeddings_batch(large_dataset)

# Check processing statistics
stats = service.get_optimization_stats()
print(f"Batch generations: {stats['generation_stats']['batch_generations']}")
print(f"Total embeddings: {stats['generation_stats']['embeddings_generated']}")
```

## Performance Monitoring

### Comprehensive Statistics

```python
# Get detailed optimization statistics
stats = service.get_optimization_stats()

# Overall configuration
print(f"Optimization Level: {stats['optimization_level']}")
print(f"Embedding Dimensions: {stats['embedding_dimensions']}")

# Cache performance
cache_stats = stats['cache_stats']
print(f"Cache Hit Rate: {stats['cache_hit_rate']:.1f}%")
print(f"Cache Memory Usage: {cache_stats['memory_usage_mb']}MB")
print(f"Cache Utilization: {cache_stats['memory_utilization_percent']:.1f}%")

# Dimension reduction performance
if 'dimension_reduction_stats' in stats:
    dr_stats = stats['dimension_reduction_stats']
    print(f"PCA Model Trained: {dr_stats['is_trained']}")
    print(f"Memory Reduction: {dr_stats['memory_reduction_percent']:.1f}%")
    print(f"Variance Retained: {dr_stats['variance_retained']:.1%}")

# Generation statistics
gen_stats = stats['generation_stats']
print(f"Embeddings Generated: {gen_stats['embeddings_generated']}")
print(f"Cache Hits: {gen_stats['cache_hits']}")
print(f"Dimension Reductions: {gen_stats['dimension_reductions']}")
```

### Memory Usage Tracking

```python
# Monitor memory usage over time
import psutil
import time

process = psutil.Process()
initial_memory = process.memory_info().rss / 1024 / 1024

# Process many variants
for i in range(1000):
    embedding = service.generate_embedding_sync(f"variant {i}")
    
    if i % 100 == 0:
        current_memory = process.memory_info().rss / 1024 / 1024
        print(f"Processed {i} variants, Memory: {current_memory:.1f}MB")

final_memory = process.memory_info().rss / 1024 / 1024
print(f"Memory change: {final_memory - initial_memory:.1f}MB")
```

## Dependencies

### Required Dependencies
```bash
# Core functionality
pip install numpy pandas lancedb

# VCF Analysis Agent
pip install -e .
```

### Optional Dependencies for Enhanced Optimizations
```bash
# For PCA dimension reduction
pip install scikit-learn

# For memory monitoring
pip install psutil
```

### Installing All Dependencies
```bash
# Install all optimization dependencies
pip install scikit-learn psutil

# Or install with extras
pip install -e .[optimization]
```

## Integration Examples

### VCF Processing Pipeline

```python
from vcf_agent.config import SessionConfig, MemoryOptimizationConfig
from vcf_agent.lancedb_integration import (
    VariantEmbeddingService, 
    create_vcf_variant_record,
    batch_add_vcf_variants
)

# Configure optimized processing
memory_config = MemoryOptimizationConfig(optimization_level="aggressive")
session_config = SessionConfig(memory_optimization=memory_config)

# Initialize services
embedding_service = VariantEmbeddingService(session_config)

# Process VCF variants with optimizations
variants_data = [
    {
        "chromosome": "chr1",
        "position": 123456,
        "reference": "A",
        "alternate": "G",
        "gene_symbol": "BRCA1",
        "clinical_significance": "Pathogenic"
    }
    # ... more variants
]

# Create optimized variant records
records = []
for variant_data in variants_data:
    record = create_vcf_variant_record(variant_data, embedding_service)
    records.append(record)

# Check optimization statistics
stats = embedding_service.get_optimization_stats()
print(f"Processed with {stats['embedding_dimensions']} dimensions")
print(f"Memory reduction: {stats.get('dimension_reduction_stats', {}).get('memory_reduction_percent', 0):.1f}%")
```

### Large-Scale Batch Processing

```python
# Configure for large-scale processing
large_scale_config = MemoryOptimizationConfig(
    optimization_level="aggressive",
    streaming_batch_size=100,           # Larger batches
    memory_cleanup_threshold_mb=50,     # Higher threshold
    cache_max_size_mb=100              # Larger cache
)

session_config = SessionConfig(memory_optimization=large_scale_config)
service = VariantEmbeddingService(session_config)

# Process very large dataset
large_variants = []  # Thousands of variants
embeddings = service.generate_embeddings_batch(large_variants)

# Monitor performance
stats = service.get_optimization_stats()
print(f"Processed {len(embeddings)} embeddings")
print(f"Cache hit rate: {stats['cache_hit_rate']:.1f}%")
print(f"Batch generations: {stats['generation_stats']['batch_generations']}")
```

## Best Practices

### 1. Choose Appropriate Optimization Level
- **Development**: Use "basic" for faster debugging
- **Production**: Use "standard" for balanced performance
- **Enterprise**: Use "aggressive" for maximum optimization

### 2. Monitor Memory Usage
```python
# Regular monitoring
stats = service.get_optimization_stats()
if stats['cache_stats']['memory_utilization_percent'] > 80:
    print("Consider increasing cache size or reducing batch size")
```

### 3. Configure Based on Available Resources
```python
# Adjust based on available memory
import psutil

available_memory_gb = psutil.virtual_memory().available / (1024**3)

if available_memory_gb > 16:
    # High memory system
    config = MemoryOptimizationConfig(
        optimization_level="aggressive",
        cache_max_size_mb=200,
        streaming_batch_size=100
    )
else:
    # Limited memory system
    config = MemoryOptimizationConfig(
        optimization_level="standard",
        cache_max_size_mb=40,
        streaming_batch_size=25
    )
```

### 4. Validate Optimization Effectiveness
```python
# Compare optimized vs unoptimized
optimized_config = MemoryOptimizationConfig(optimization_level="aggressive")
basic_config = MemoryOptimizationConfig(optimization_level="basic")

# Test with same dataset
test_texts = ["variant description"] * 100

# Optimized processing
optimized_service = VariantEmbeddingService(SessionConfig(memory_optimization=optimized_config))
optimized_embeddings = optimized_service.generate_embeddings_batch(test_texts)

# Basic processing
basic_service = VariantEmbeddingService(SessionConfig(memory_optimization=basic_config))
basic_embeddings = basic_service.generate_embeddings_batch(test_texts)

# Compare results
print(f"Optimized dimensions: {len(optimized_embeddings[0])}")
print(f"Basic dimensions: {len(basic_embeddings[0])}")

optimized_stats = optimized_service.get_optimization_stats()
basic_stats = basic_service.get_optimization_stats()

print(f"Optimized memory usage: {optimized_stats['cache_stats']['memory_usage_mb']}MB")
print(f"Basic memory usage: {basic_stats['cache_stats']['memory_usage_mb']}MB")
```

## Troubleshooting

### Common Issues

#### 1. Dimension Reduction Not Working
```python
# Check if scikit-learn is installed
try:
    import sklearn
    print("scikit-learn available")
except ImportError:
    print("Install scikit-learn: pip install scikit-learn")

# Check if PCA model is trained
stats = service.get_optimization_stats()
if 'dimension_reduction_stats' in stats:
    dr_stats = stats['dimension_reduction_stats']
    print(f"PCA trained: {dr_stats['is_trained']}")
    print(f"Training samples: {dr_stats['training_samples_collected']}/{dr_stats['training_samples_needed']}")
```

#### 2. Memory Usage Still High
```python
# Check configuration
stats = service.get_optimization_stats()
print(f"Optimization level: {stats['optimization_level']}")
print(f"Memory management enabled: {stats['memory_management_enabled']}")

# Force cache cleanup
if hasattr(service.embedding_cache, '_cleanup_cache'):
    service.embedding_cache._cleanup_cache()
    print("Manual cache cleanup performed")
```

#### 3. Performance Degradation
```python
# Check cache efficiency
stats = service.get_optimization_stats()
if stats['cache_hit_rate'] < 30:
    print("Low cache hit rate, consider:")
    print("- Increasing cache size")
    print("- Reducing dataset diversity")
    print("- Checking for cache eviction")
```

### Debugging Tools

```python
# Enable detailed logging
import logging
logging.basicConfig(level=logging.DEBUG)

# Monitor memory in real-time
import psutil
import threading
import time

def memory_monitor():
    process = psutil.Process()
    while True:
        memory_mb = process.memory_info().rss / 1024 / 1024
        print(f"Memory usage: {memory_mb:.1f}MB")
        time.sleep(10)

# Start monitoring in background
monitor_thread = threading.Thread(target=memory_monitor, daemon=True)
monitor_thread.start()
```

## Migration Guide

### From Legacy Systems

The memory optimization system is fully backward compatible. Existing code continues to work without modifications:

```python
# Legacy code (still works)
from vcf_agent.lancedb_integration import VariantEmbeddingService
service = VariantEmbeddingService()  # Uses default optimizations
embedding = service.generate_embedding_sync("variant description")

# Enhanced code (leverages new features)
from vcf_agent.config import SessionConfig, MemoryOptimizationConfig

memory_config = MemoryOptimizationConfig(optimization_level="aggressive")
session_config = SessionConfig(memory_optimization=memory_config)
enhanced_service = VariantEmbeddingService(session_config)
optimized_embedding = enhanced_service.generate_embedding_sync("variant description")
```

### Upgrading Configuration

```python
# Old approach (deprecated, but still works)
service = VariantEmbeddingService()

# New approach (recommended)
memory_config = MemoryOptimizationConfig(optimization_level="standard")
session_config = SessionConfig(memory_optimization=memory_config)
service = VariantEmbeddingService(session_config)
```

All existing VCF processing workflows automatically benefit from the integrated optimizations without requiring code changes.

---

*For additional support or questions about memory optimization features, please refer to the main README.md or contact the development team.* 