# Configuration Guide

This guide covers all configuration options for the VCF Analysis Agent, from basic setup to advanced production configurations.

## Overview

The VCF Analysis Agent can be configured through multiple methods:

1. **Environment Variables** - For system-wide settings
2. **Configuration Files** - For structured configuration
3. **Command Line Arguments** - For runtime overrides
4. **Session Configuration** - For programmatic control

## Environment Variables

### Core Configuration

```bash
# Database Paths
LANCEDB_PATH=./lancedb                    # LanceDB storage location
KUZU_DB_PATH=./kuzu_db                    # Kuzu database location

# Logging
LOG_LEVEL=INFO                            # DEBUG, INFO, WARNING, ERROR
VCF_AGENT_RAW_MODE=false                  # Enable/disable chain-of-thought

# Performance
VCF_AGENT_BATCH_SIZE=5000                 # Default batch size for processing
VCF_AGENT_MAX_WORKERS=4                   # Maximum concurrent workers
```

### AI Provider Configuration

```bash
# OpenAI
OPENAI_API_KEY=sk-your-key-here
OPENAI_MODEL=gpt-4                        # Default model
OPENAI_BASE_URL=https://api.openai.com/v1 # Custom endpoint (optional)

# Cerebras
CEREBRAS_API_KEY=your-cerebras-key
CEREBRAS_MODEL=cerebras-gpt               # Default model

# Ollama
OLLAMA_BASE_URL=http://localhost:11434    # Ollama server URL
OLLAMA_MODEL=qwen3:4b                     # Default local model
```

### Monitoring & Observability

```bash
# Prometheus Metrics
VCF_AGENT_METRICS_PORT=8000               # Metrics server port
VCF_AGENT_PUSHGATEWAY_URL=http://localhost:9091  # Pushgateway for CLI metrics

# OpenTelemetry Tracing
OTEL_EXPORTER_OTLP_TRACES_ENDPOINT=http://localhost:4317
OTEL_SERVICE_NAME=vcf-agent-core
OTEL_RESOURCE_ATTRIBUTES=service.version=1.0.0

# Jaeger (alternative to OTLP)
JAEGER_AGENT_HOST=localhost
JAEGER_AGENT_PORT=6831
```

### Security & Access Control

```bash
# API Security
VCF_AGENT_API_KEY=your-secure-api-key     # For API access
VCF_AGENT_ALLOWED_ORIGINS=localhost,*.example.com  # CORS origins

# Data Protection
VCF_AGENT_ENCRYPT_AT_REST=true            # Enable database encryption
VCF_AGENT_AUDIT_LOG=true                  # Enable audit logging
```

## Configuration Files

### Main Configuration File

Create `config/vcf_agent.yaml`:

```yaml
# VCF Analysis Agent Configuration
version: "1.0"

# Database Configuration
databases:
  lancedb:
    path: "./lancedb"
    cache_size: "1GB"
    enable_compression: true
  
  kuzu:
    path: "./kuzu_db"
    buffer_pool_size: "512MB"
    max_num_threads: 4

# AI Providers
ai_providers:
  default: "ollama"
  
  openai:
    model: "gpt-4"
    max_tokens: 4096
    temperature: 0.1
    timeout: 30
  
  cerebras:
    model: "cerebras-gpt"
    max_tokens: 2048
    temperature: 0.1
  
  ollama:
    base_url: "http://localhost:11434"
    model: "qwen3:4b"
    timeout: 60

# Performance Settings
performance:
  batch_size: 5000
  max_workers: 4
  enable_caching: true
  cache_ttl: 3600  # seconds
  
  # Memory management
  max_memory_usage: "8GB"
  enable_streaming: true
  
  # Optimization features
  enable_embedding_cache: true
  embedding_cache_size: 10000
  enable_query_batching: true
  batch_query_size: 50

# Monitoring
monitoring:
  enable_metrics: true
  metrics_port: 8000
  enable_tracing: true
  enable_profiling: false
  
  # Health checks
  health_check_interval: 30
  database_health_check: true

# Security
security:
  enable_audit_log: true
  audit_log_path: "./logs/audit.log"
  encrypt_at_rest: false
  require_api_key: false
  
  # Rate limiting
  enable_rate_limiting: true
  requests_per_minute: 100

# Logging
logging:
  level: "INFO"
  format: "json"
  output: "stdout"
  
  # File logging
  enable_file_logging: true
  log_file: "./logs/vcf_agent.log"
  max_file_size: "100MB"
  backup_count: 5
```

### Docker Configuration

For Docker deployments, use `docker-compose.yml`:

```yaml
version: '3.8'

services:
  vcf-agent:
    build: .
    ports:
      - "8000:8000"  # Metrics
      - "8080:8080"  # API
    environment:
      - LOG_LEVEL=INFO
      - VCF_AGENT_METRICS_PORT=8000
      - LANCEDB_PATH=/data/lancedb
      - KUZU_DB_PATH=/data/kuzu_db
    volumes:
      - ./data:/data
      - ./config:/app/config
      - ./logs:/app/logs
    depends_on:
      - prometheus
      - jaeger

  prometheus:
    image: prom/prometheus:latest
    ports:
      - "9090:9090"
    volumes:
      - ./config/prometheus.yml:/etc/prometheus/prometheus.yml

  jaeger:
    image: jaegertracing/all-in-one:latest
    ports:
      - "16686:16686"
      - "14268:14268"
    environment:
      - COLLECTOR_OTLP_ENABLED=true
```

## Session Configuration

For programmatic control, use the `SessionConfig` class:

```python
from vcf_agent.config import SessionConfig

# Basic configuration
config = SessionConfig(
    model_provider="ollama",
    ollama_model_name="qwen3:4b",
    raw_mode=False
)

# Advanced configuration
config = SessionConfig(
    model_provider="openai",
    credentials_file="./credentials.json",
    reference_fasta="./reference/hg38.fa",
    raw_mode=True
)

# Use with agent
from vcf_agent.agent import get_agent_with_session
agent = get_agent_with_session(session_config=config)
```

## Command Line Configuration

Override settings via command line arguments:

```bash
# Basic analysis with custom settings
vcf-agent analyze sample.vcf \
  --provider openai \
  --model gpt-4 \
  --batch-size 1000 \
  --output results.json

# Advanced configuration
vcf-agent analyze large.vcf \
  --provider ollama \
  --streaming \
  --memory-limit 16GB \
  --workers 8 \
  --enable-cache \
  --verbose
```

## Provider-Specific Configuration

### OpenAI Configuration

```bash
# Environment variables
export OPENAI_API_KEY="sk-your-key-here"
export OPENAI_MODEL="gpt-4"
export OPENAI_MAX_TOKENS=4096
export OPENAI_TEMPERATURE=0.1

# Or in credentials.json
{
  "openai": {
    "api_key": "sk-your-key-here",
    "model": "gpt-4",
    "max_tokens": 4096,
    "temperature": 0.1,
    "organization": "org-your-org-id"
  }
}
```

### Cerebras Configuration

```bash
# Environment variables
export CEREBRAS_API_KEY="your-cerebras-key"
export CEREBRAS_MODEL="cerebras-gpt"

# Or in credentials.json
{
  "cerebras": {
    "api_key": "your-cerebras-key",
    "model": "cerebras-gpt",
    "max_tokens": 2048
  }
}
```

### Ollama Configuration

```bash
# Environment variables
export OLLAMA_BASE_URL="http://localhost:11434"
export OLLAMA_MODEL="qwen3:4b"

# Pull models
ollama pull qwen3:4b
ollama pull llama2:7b
ollama pull codellama:13b

# List available models
ollama list
```

## Database Configuration

### LanceDB Settings

```python
# Advanced LanceDB configuration
lancedb_config = {
    "path": "./lancedb",
    "cache_size": "2GB",
    "enable_compression": True,
    "compression_level": 6,
    "index_type": "IVF_PQ",
    "num_partitions": 256,
    "num_sub_quantizers": 96
}
```

### Kuzu Settings

```python
# Advanced Kuzu configuration
kuzu_config = {
    "path": "./kuzu_db",
    "buffer_pool_size": "1GB",
    "max_num_threads": 8,
    "enable_compression": True,
    "checkpoint_threshold": 1000000
}
```

## Performance Tuning

### Memory Optimization

```yaml
performance:
  # Batch processing
  batch_size: 10000              # Larger batches for better throughput
  max_workers: 8                 # Match CPU cores
  
  # Memory limits
  max_memory_usage: "16GB"       # Set based on available RAM
  enable_streaming: true         # For large files
  
  # Caching
  enable_embedding_cache: true
  embedding_cache_size: 50000    # Increase for better hit rates
  cache_ttl: 7200               # 2 hours
```

### Network Optimization

```yaml
ai_providers:
  openai:
    timeout: 60                  # Increase for large requests
    max_retries: 3
    retry_delay: 1
    
  ollama:
    timeout: 120                 # Local models may need more time
    keep_alive: "5m"            # Keep model loaded
```

## Security Configuration

### API Security

```yaml
security:
  # Authentication
  require_api_key: true
  api_key_header: "X-VCF-Agent-Key"
  
  # Authorization
  enable_rbac: true
  roles_file: "./config/roles.yaml"
  
  # Rate limiting
  enable_rate_limiting: true
  requests_per_minute: 100
  burst_limit: 20
```

### Data Protection

```yaml
security:
  # Encryption
  encrypt_at_rest: true
  encryption_key_file: "./config/encryption.key"
  
  # Audit logging
  enable_audit_log: true
  audit_log_path: "./logs/audit.log"
  audit_events: ["access", "modify", "delete"]
  
  # Data retention
  data_retention_days: 90
  auto_cleanup: true
```

## Monitoring Configuration

### Prometheus Metrics

```yaml
monitoring:
  prometheus:
    enabled: true
    port: 8000
    path: "/metrics"
    
    # Custom metrics
    custom_metrics:
      - name: "vcf_processing_time"
        type: "histogram"
        help: "Time to process VCF files"
```

### OpenTelemetry Tracing

```yaml
monitoring:
  tracing:
    enabled: true
    exporter: "otlp"
    endpoint: "http://localhost:4317"
    
    # Sampling
    sampling_rate: 0.1
    max_spans_per_trace: 1000
```

## Troubleshooting Configuration

### Debug Mode

```bash
# Enable debug logging
export LOG_LEVEL=DEBUG
export VCF_AGENT_DEBUG=true

# Run with verbose output
vcf-agent analyze sample.vcf --verbose --debug
```

### Configuration Validation

```bash
# Validate configuration
vcf-agent config validate

# Show current configuration
vcf-agent config show

# Test database connections
vcf-agent config test-db
```

### Common Issues

**Configuration not loading**
```bash
# Check configuration file path
vcf-agent config show --config-file ./config/vcf_agent.yaml

# Validate YAML syntax
python -c "import yaml; yaml.safe_load(open('config/vcf_agent.yaml'))"
```

**Environment variables not working**
```bash
# Check environment
env | grep VCF_AGENT

# Test specific variable
echo $OPENAI_API_KEY
```

## Best Practices

1. **Use Configuration Files**: For complex setups, prefer YAML configuration over environment variables
2. **Secure Credentials**: Store API keys in separate credential files, not in main configuration
3. **Environment-Specific Configs**: Use different configurations for development, staging, and production
4. **Monitor Resource Usage**: Configure appropriate memory and CPU limits
5. **Enable Logging**: Use structured logging for better observability
6. **Regular Backups**: Backup configuration files and database paths
7. **Version Control**: Track configuration changes in version control (excluding secrets)

---

**Configuration Complete!** 🎉

Your VCF Analysis Agent is now configured for optimal performance and security. For production deployments, review the [Deployment Guide](deployment.md) for additional considerations. 