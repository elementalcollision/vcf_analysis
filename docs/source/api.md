# VCF Agent API Reference

## Overview

This document provides comprehensive API reference for the VCF Agent data stores implementation, including the UnifiedDataStoreManager, LanceDB integration, Kuzu graph database, and embedding services.

## Core Classes

### UnifiedDataStoreManager

The central orchestrator for all data operations, providing a unified interface for both LanceDB and Kuzu databases.

```python
from src.vcf_agent.data_store_manager import UnifiedDataStoreManager, DataStoreConfig
from src.vcf_agent.config import SessionConfig
```

#### Constructor

```python
def __init__(
    self, 
    config: Optional[DataStoreConfig] = None, 
    session_config: Optional[SessionConfig] = None
)
```

**Parameters:**
- `config`: Data store configuration (optional, uses defaults if not provided)
- `session_config`: Session configuration for AI models (optional)

**Example:**
```python
config = DataStoreConfig(
    lancedb_path="./data/lancedb",
    kuzu_path="./data/kuzu_db",
    sync_batch_size=2000,
    max_workers=8
)

manager = UnifiedDataStoreManager(config=config)
```

#### Methods

##### add_sample_with_variants

```python
def add_sample_with_variants(
    self,
    sample_data: Dict[str, Any],
    variants_data: List[Dict[str, Any]],
    genes_data: Optional[List[Dict[str, Any]]] = None,
    analysis_data: Optional[List[Dict[str, Any]]] = None
) -> Dict[str, Any]
```

Add a sample with its variants to both databases with full synchronization.

**Parameters:**
- `sample_data`: Sample information dictionary
- `variants_data`: List of variant dictionaries
- `genes_data`: Optional list of gene dictionaries
- `analysis_data`: Optional list of analysis dictionaries

**Returns:**
Dictionary with operation results and performance metrics.

**Example:**
```python
sample_data = {
    "id": "SAMPLE_001",
    "name": "Patient Sample 1",
    "type": "germline",
    "metadata": '{"patient_id": "P001"}'
}

variants_data = [
    {
        "id": "chr1-123456-A-G",
        "chr": "1",
        "pos": 123456,
        "ref": "A",
        "alt": "G",
        "genotype": "0/1",
        "quality_score": 99.5,
        "clinical_significance": "Pathogenic"
    }
]

result = manager.add_sample_with_variants(
    sample_data=sample_data,
    variants_data=variants_data
)
```

##### search_variants

```python
def search_variants(
    self,
    query: str,
    search_type: str = "hybrid",
    metadata_filter: Optional[str] = None,
    limit: int = 10,
    similarity_threshold: float = 0.7
) -> Dict[str, Any]
```

Search variants using multiple search strategies.

**Parameters:**
- `query`: Search query text
- `search_type`: Type of search ("hybrid", "vector", "graph")
- `metadata_filter`: Optional metadata filter for hybrid search
- `limit`: Maximum number of results
- `similarity_threshold`: Minimum similarity score

**Returns:**
Dictionary with search results and metadata.

**Example:**
```python
results = manager.search_variants(
    query="pathogenic variant in BRCA1 gene",
    search_type="hybrid",
    metadata_filter="chromosome = '17' AND gene_symbol = 'BRCA1'",
    limit=10,
    similarity_threshold=0.8
)
```

##### get_sample_analysis

```python
def get_sample_analysis(self, sample_id: str) -> Dict[str, Any]
```

Get comprehensive analysis for a sample using both databases.

**Parameters:**
- `sample_id`: ID of the sample to analyze

**Returns:**
Dictionary with comprehensive sample analysis.

**Example:**
```python
analysis = manager.get_sample_analysis("SAMPLE_001")
print(f"Total variants: {analysis['statistics']['total_variants']}")
```

##### get_comprehensive_statistics

```python
def get_comprehensive_statistics(self) -> Dict[str, Any]
```

Get comprehensive statistics from both databases.

**Returns:**
Dictionary with statistics from both LanceDB and Kuzu.

**Example:**
```python
stats = manager.get_comprehensive_statistics()
print(f"LanceDB variants: {stats['lancedb_stats']['total_variants']}")
print(f"Kuzu nodes: {stats['kuzu_stats']['node_counts']}")
```

### DataStoreConfig

Configuration class for the unified data store manager.

```python
@dataclass
class DataStoreConfig:
    lancedb_path: str = "./lancedb"
    kuzu_path: str = "./kuzu_db"
    lancedb_table_name: str = "vcf_variants"
    enable_sync: bool = True
    sync_batch_size: int = 1000
    performance_monitoring: bool = True
    max_workers: int = 4
```

**Example:**
```python
config = DataStoreConfig(
    lancedb_path="/data/lancedb",
    kuzu_path="/data/kuzu",
    sync_batch_size=2000,
    max_workers=8,
    performance_monitoring=True
)
```

### PerformanceMetrics

Performance tracking for data store operations.

```python
@dataclass
class PerformanceMetrics:
    operation: str
    start_time: float
    end_time: float
    duration: float
    records_processed: int
    success: bool
    error_message: Optional[str] = None
    
    @property
    def records_per_second(self) -> float
```

## LanceDB Integration

### VCFVariant Model

Enhanced Pydantic model for VCF variant records.

```python
from src.vcf_agent.lancedb_integration import VCFVariant
```

```python
class VCFVariant(LanceModel):
    variant_id: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    variant_description: str
    variant_vector: Vector(1536)  # OpenAI text-embedding-3-small dimension
    analysis_summary: str
    sample_id: str
    quality_score: Optional[float] = None
    filter_status: Optional[str] = None
    genotype: Optional[str] = None
    allele_frequency: Optional[float] = None
    clinical_significance: Optional[str] = None
    gene_symbol: Optional[str] = None
    consequence: Optional[str] = None
    created_at: datetime
    updated_at: datetime
```

### VariantEmbeddingService

Service for generating embeddings for VCF variants.

```python
from src.vcf_agent.lancedb_integration import VariantEmbeddingService
```

#### Constructor

```python
def __init__(self, session_config: Optional[SessionConfig] = None)
```

#### Methods

##### generate_variant_description

```python
def generate_variant_description(self, variant_data: Dict) -> str
```

Generate a human-readable description for a variant.

**Example:**
```python
service = VariantEmbeddingService()
variant_data = {
    "chromosome": "1",
    "position": 123456,
    "reference": "A",
    "alternate": "G",
    "gene_symbol": "BRCA1"
}
description = service.generate_variant_description(variant_data)
# Returns: "Variant at chromosome 1 position 123456 reference allele A to alternate allele G in gene BRCA1"
```

##### generate_embedding_sync

```python
def generate_embedding_sync(self, text: str) -> List[float]
```

Synchronous wrapper for embedding generation.

**Example:**
```python
embedding = service.generate_embedding_sync("pathogenic variant in BRCA1")
# Returns: [0.1, 0.2, 0.3, ...] (1536-dimensional vector)
```

### Core Functions

#### get_or_create_vcf_table

```python
def get_or_create_vcf_table(db, table_name: str = "vcf_variants") -> LanceTable
```

Opens an existing VCF variants LanceDB table or creates a new one.

**Example:**
```python
from src.vcf_agent.lancedb_integration import get_db, get_or_create_vcf_table

db = get_db("./lancedb")
table = get_or_create_vcf_table(db, "vcf_variants")
```

#### batch_add_vcf_variants

```python
def batch_add_vcf_variants(
    table: LanceTable,
    variants_data: List[Dict],
    embedding_service: Optional[VariantEmbeddingService] = None,
    batch_size: int = 1000,
    max_workers: int = 4
) -> int
```

Add a batch of VCF variant records with optimized performance.

**Example:**
```python
variants_data = [
    {
        "chromosome": "1",
        "position": 123456,
        "reference": "A",
        "alternate": "G",
        "sample_id": "SAMPLE_001"
    }
]

count = batch_add_vcf_variants(
    table=table,
    variants_data=variants_data,
    batch_size=500,
    max_workers=8
)
```

#### hybrid_search_variants

```python
def hybrid_search_variants(
    table: LanceTable,
    query_text: str,
    embedding_service: Optional[VariantEmbeddingService] = None,
    metadata_filter: Optional[str] = None,
    limit: int = 10,
    similarity_threshold: float = 0.7
) -> pd.DataFrame
```

Perform hybrid search combining vector similarity and metadata filtering.

**Example:**
```python
results = hybrid_search_variants(
    table=table,
    query_text="pathogenic BRCA1 variant",
    metadata_filter="chromosome = '17'",
    limit=5,
    similarity_threshold=0.8
)
```

#### search_similar_variants

```python
def search_similar_variants(
    table: LanceTable,
    reference_variant_id: str,
    limit: int = 10,
    include_metadata: bool = True
) -> pd.DataFrame
```

Find variants similar to a reference variant using vector similarity.

**Example:**
```python
similar = search_similar_variants(
    table=table,
    reference_variant_id="chr1-123456-A-G",
    limit=10
)
```

## Kuzu Graph Integration

### Core Functions

#### get_kuzu_db_connection

```python
def get_kuzu_db_connection(
    db_path: str = "./kuzu_db", 
    read_only: bool = False
) -> kuzu.Connection
```

Initialize a Kuzu database connection.

**Example:**
```python
from src.vcf_agent.graph_integration import get_kuzu_db_connection

conn = get_kuzu_db_connection("./kuzu_db")
```

#### create_enhanced_schema

```python
def create_enhanced_schema(conn: kuzu.Connection) -> None
```

Creates the enhanced schema for genomic data.

**Example:**
```python
from src.vcf_agent.graph_integration import create_enhanced_schema

create_enhanced_schema(conn)
```

#### batch_add_genomic_data

```python
def batch_add_genomic_data(
    conn: kuzu.Connection,
    samples: List[Dict[str, Any]],
    variants: List[Dict[str, Any]],
    genes: List[Dict[str, Any]],
    relationships: List[Dict[str, Any]],
    batch_size: int = 1000
) -> Dict[str, int]
```

Batch add genomic data with optimized performance.

**Example:**
```python
samples = [{"id": "SAMPLE_001", "name": "Patient 1", "type": "germline"}]
variants = [{"id": "chr1-123456-A-G", "chr": "1", "pos": 123456, "ref": "A", "alt": "G"}]
genes = [{"id": "ENSG00000012048", "symbol": "BRCA1", "chromosome": "17"}]
relationships = [
    {
        "type": "HasVariant",
        "sample_id": "SAMPLE_001",
        "variant_id": "chr1-123456-A-G",
        "properties": {"genotype": "0/1", "quality": 99.5}
    }
]

counts = batch_add_genomic_data(conn, samples, variants, genes, relationships)
```

#### find_sample_variants

```python
def find_sample_variants(
    conn: kuzu.Connection, 
    sample_id: str, 
    limit: int = 100
) -> List[Dict[str, Any]]
```

Find all variants for a specific sample.

**Example:**
```python
variants = find_sample_variants(conn, "SAMPLE_001", limit=50)
```

#### find_variant_genes

```python
def find_variant_genes(
    conn: kuzu.Connection, 
    variant_id: str
) -> List[Dict[str, Any]]
```

Find all genes associated with a specific variant.

**Example:**
```python
genes = find_variant_genes(conn, "chr1-123456-A-G")
```

#### find_similar_samples

```python
def find_similar_samples(
    conn: kuzu.Connection, 
    sample_id: str, 
    min_shared_variants: int = 5
) -> List[Dict[str, Any]]
```

Find samples that share variants with the given sample.

**Example:**
```python
similar_samples = find_similar_samples(conn, "SAMPLE_001", min_shared_variants=10)
```

#### get_genomic_statistics

```python
def get_genomic_statistics(conn: kuzu.Connection) -> Dict[str, Any]
```

Get comprehensive statistics about the genomic data in the graph.

**Example:**
```python
stats = get_genomic_statistics(conn)
print(f"Total samples: {stats['node_counts']['sample']}")
print(f"Total variants: {stats['node_counts']['variant']}")
```

## Convenience Functions

### create_data_store_manager

```python
def create_data_store_manager(
    lancedb_path: str = "./lancedb",
    kuzu_path: str = "./kuzu_db",
    session_config: Optional[SessionConfig] = None
) -> UnifiedDataStoreManager
```

Create a unified data store manager with the specified configuration.

**Example:**
```python
from src.vcf_agent.data_store_manager import create_data_store_manager

manager = create_data_store_manager(
    lancedb_path="/data/lancedb",
    kuzu_path="/data/kuzu_db"
)
```

## Error Handling

### Common Exceptions

#### DatabaseConnectionError
Raised when database connections fail.

```python
try:
    manager = create_data_store_manager()
except DatabaseConnectionError as e:
    logger.error(f"Failed to connect to databases: {e}")
```

#### ValidationError
Raised when data validation fails.

```python
try:
    result = manager.add_sample_with_variants(invalid_data)
except ValidationError as e:
    logger.error(f"Data validation failed: {e}")
```

#### PerformanceError
Raised when performance targets are not met.

```python
try:
    result = manager.search_variants(query, timeout=100)
except PerformanceError as e:
    logger.warning(f"Performance target missed: {e}")
```

## Performance Optimization

### Batch Size Tuning

```python
# For high-throughput ingestion
config = DataStoreConfig(
    sync_batch_size=5000,  # Larger batches for better throughput
    max_workers=16         # More workers for parallel processing
)

# For memory-constrained environments
config = DataStoreConfig(
    sync_batch_size=500,   # Smaller batches to reduce memory usage
    max_workers=2          # Fewer workers to reduce resource usage
)
```

### Connection Pooling

```python
# Reuse manager instances for better performance
class DataStorePool:
    def __init__(self):
        self.manager = create_data_store_manager()
    
    def get_manager(self):
        return self.manager

# Global pool instance
pool = DataStorePool()
manager = pool.get_manager()
```

### Monitoring Performance

```python
# Enable detailed performance monitoring
config = DataStoreConfig(performance_monitoring=True)
manager = UnifiedDataStoreManager(config)

# Get performance statistics
stats = manager._get_performance_statistics()
for operation, metrics in stats.items():
    print(f"{operation}: {metrics['avg_duration']:.3f}s avg")
```

## Integration Examples

### Complete Workflow Example

```python
from src.vcf_agent.data_store_manager import create_data_store_manager
from src.vcf_agent.config import SessionConfig

# 1. Initialize manager
session_config = SessionConfig(model_provider="openai")
manager = create_data_store_manager(session_config=session_config)

# 2. Add sample data
sample_data = {
    "id": "SAMPLE_001",
    "name": "Patient Sample",
    "type": "germline"
}

variants_data = [
    {
        "id": "chr17-43044295-G-A",
        "chr": "17",
        "pos": 43044295,
        "ref": "G",
        "alt": "A",
        "gene_symbol": "BRCA1",
        "clinical_significance": "Pathogenic",
        "genotype": "0/1"
    }
]

# 3. Ingest data
result = manager.add_sample_with_variants(
    sample_data=sample_data,
    variants_data=variants_data
)

print(f"Ingested {result['performance']['total_variants']} variants")

# 4. Search for similar variants
search_results = manager.search_variants(
    query="pathogenic BRCA1 mutation",
    search_type="hybrid",
    limit=10
)

print(f"Found {len(search_results['results'])} similar variants")

# 5. Analyze sample
analysis = manager.get_sample_analysis("SAMPLE_001")
print(f"Sample has {analysis['statistics']['total_variants']} variants")

# 6. Get comprehensive statistics
stats = manager.get_comprehensive_statistics()
print(f"Database contains {stats['lancedb_stats']['total_variants']} total variants")
```

### Async Operations Example

```python
import asyncio
from concurrent.futures import ThreadPoolExecutor

async def process_multiple_samples(manager, samples_data):
    """Process multiple samples concurrently."""
    
    async def process_sample(sample_data):
        loop = asyncio.get_event_loop()
        with ThreadPoolExecutor() as executor:
            result = await loop.run_in_executor(
                executor,
                manager.add_sample_with_variants,
                sample_data['sample'],
                sample_data['variants']
            )
        return result
    
    # Process all samples concurrently
    tasks = [process_sample(sample) for sample in samples_data]
    results = await asyncio.gather(*tasks)
    
    return results

# Usage
samples_data = [
    {
        "sample": {"id": "SAMPLE_001", "name": "Patient 1"},
        "variants": [{"id": "var1", "chr": "1", "pos": 123, "ref": "A", "alt": "G"}]
    },
    {
        "sample": {"id": "SAMPLE_002", "name": "Patient 2"},
        "variants": [{"id": "var2", "chr": "2", "pos": 456, "ref": "T", "alt": "C"}]
    }
]

results = asyncio.run(process_multiple_samples(manager, samples_data))
```

## Testing

### Unit Testing Example

```python
import pytest
from unittest.mock import Mock, patch
from src.vcf_agent.data_store_manager import UnifiedDataStoreManager

class TestUnifiedDataStoreManager:
    
    @pytest.fixture
    def manager(self):
        with patch('src.vcf_agent.data_store_manager.get_lancedb'), \
             patch('src.vcf_agent.data_store_manager.get_kuzu_db_connection'):
            return UnifiedDataStoreManager()
    
    def test_add_sample_with_variants(self, manager):
        sample_data = {"id": "TEST_SAMPLE", "name": "Test"}
        variants_data = [{"id": "test_var", "chr": "1", "pos": 123}]
        
        result = manager.add_sample_with_variants(sample_data, variants_data)
        
        assert result["sample_id"] == "TEST_SAMPLE"
        assert result["sync_status"] == "success"
```

### Integration Testing Example

```python
import tempfile
import shutil
from src.vcf_agent.data_store_manager import create_data_store_manager

class TestIntegration:
    
    @pytest.fixture
    def temp_manager(self):
        temp_dir = tempfile.mkdtemp()
        manager = create_data_store_manager(
            lancedb_path=f"{temp_dir}/lancedb",
            kuzu_path=f"{temp_dir}/kuzu"
        )
        yield manager
        shutil.rmtree(temp_dir)
    
    def test_end_to_end_workflow(self, temp_manager):
        # Test complete workflow
        sample_data = {"id": "INTEGRATION_TEST", "name": "Test Sample"}
        variants_data = [{"id": "test_variant", "chr": "1", "pos": 12345}]
        
        # Add data
        result = temp_manager.add_sample_with_variants(sample_data, variants_data)
        assert result["sync_status"] == "success"
        
        # Search data
        search_results = temp_manager.search_variants("test variant")
        assert len(search_results["results"]) >= 0
        
        # Analyze sample
        analysis = temp_manager.get_sample_analysis("INTEGRATION_TEST")
        assert analysis["sample_id"] == "INTEGRATION_TEST"
```

## Best Practices

### 1. Resource Management

```python
# Always close managers when done
try:
    manager = create_data_store_manager()
    # Use manager...
finally:
    manager.close()

# Or use context manager pattern
class DataStoreContext:
    def __enter__(self):
        self.manager = create_data_store_manager()
        return self.manager
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.manager.close()

# Usage
with DataStoreContext() as manager:
    result = manager.add_sample_with_variants(sample_data, variants_data)
```

### 2. Error Handling

```python
import logging

logger = logging.getLogger(__name__)

def safe_data_operation(manager, operation_func, *args, **kwargs):
    """Safely execute data operations with proper error handling."""
    try:
        return operation_func(*args, **kwargs)
    except ValidationError as e:
        logger.error(f"Data validation failed: {e}")
        return {"error": "validation_failed", "message": str(e)}
    except PerformanceError as e:
        logger.warning(f"Performance target missed: {e}")
        return {"error": "performance_warning", "message": str(e)}
    except Exception as e:
        logger.exception(f"Unexpected error in data operation: {e}")
        return {"error": "unexpected_error", "message": str(e)}

# Usage
result = safe_data_operation(
    manager,
    manager.add_sample_with_variants,
    sample_data,
    variants_data
)
```

### 3. Performance Monitoring

```python
import time
from contextlib import contextmanager

@contextmanager
def performance_monitor(operation_name):
    """Monitor performance of operations."""
    start_time = time.time()
    try:
        yield
    finally:
        duration = time.time() - start_time
        logger.info(f"{operation_name} completed in {duration:.2f}s")

# Usage
with performance_monitor("sample_ingestion"):
    result = manager.add_sample_with_variants(sample_data, variants_data)
```

### 4. Configuration Management

```python
import os
from src.vcf_agent.data_store_manager import DataStoreConfig

def create_production_config():
    """Create production-optimized configuration."""
    return DataStoreConfig(
        lancedb_path=os.getenv("LANCEDB_PATH", "/data/lancedb"),
        kuzu_path=os.getenv("KUZU_PATH", "/data/kuzu"),
        sync_batch_size=int(os.getenv("BATCH_SIZE", "2000")),
        max_workers=int(os.getenv("MAX_WORKERS", "8")),
        performance_monitoring=True
    )

def create_development_config():
    """Create development-optimized configuration."""
    return DataStoreConfig(
        lancedb_path="./dev_data/lancedb",
        kuzu_path="./dev_data/kuzu",
        sync_batch_size=100,
        max_workers=2,
        performance_monitoring=True
    )

# Usage
config = create_production_config() if os.getenv("ENV") == "production" else create_development_config()
manager = UnifiedDataStoreManager(config)
```

## Migration and Versioning

### Schema Migration

```python
def migrate_schema_v1_to_v2(manager):
    """Migrate from schema v1 to v2."""
    # Add new fields to existing records
    # Update relationship structures
    # Reindex data if necessary
    pass

def check_schema_version(manager):
    """Check current schema version."""
    stats = manager.get_comprehensive_statistics()
    return stats.get("schema_version", "1.0")

# Usage
current_version = check_schema_version(manager)
if current_version < "2.0":
    migrate_schema_v1_to_v2(manager)
```

### Data Export/Import

```python
def export_data(manager, output_path):
    """Export data for backup or migration."""
    stats = manager.get_comprehensive_statistics()
    
    # Export LanceDB data
    # Export Kuzu data
    # Save metadata
    
    return {"exported_records": stats["total_records"]}

def import_data(manager, input_path):
    """Import data from backup or migration."""
    # Load and validate data
    # Import to both databases
    # Verify consistency
    
    return {"imported_records": 0}
```

This comprehensive API reference provides detailed documentation for all components of the VCF Agent data stores implementation, including usage examples, best practices, and integration patterns. 