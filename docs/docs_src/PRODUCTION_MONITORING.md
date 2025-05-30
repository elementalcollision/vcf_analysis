# Production Monitoring & Observability

**Last Updated**: May 29, 2025  
**Status**: Production Ready ✅  
**Observability Coverage**: 100%

## Overview

The VCF Analysis Agent includes comprehensive production monitoring designed for enterprise genomic workloads. The observability stack provides real-time metrics, distributed tracing, and proactive alerting across all system components.

## Monitoring Stack

### Core Components
- **Prometheus**: Metrics collection and alerting
- **Grafana**: Visualization and dashboards
- **Jaeger**: Distributed tracing
- **OpenTelemetry**: Instrumentation and data collection

### Access Points
```yaml
Service URLs (Production):
  Grafana Dashboard: http://localhost:3000
  Prometheus Metrics: http://localhost:9090
  Jaeger Tracing: http://localhost:16686
  VCF Agent API: http://localhost:8080
```

## VCF Agent Dashboard

### Key Metrics

Access real-time monitoring at `http://localhost:3000` after production deployment:

```yaml
VCF Agent Dashboard Metrics:
  Request Processing:
    - Request rate (requests/second)
    - Error percentage and status codes
    - Response time percentiles (P50, P95, P99)
    
  VCF Processing Performance:
    - Variants processed per second
    - File processing latency
    - Batch processing efficiency
    
  AI & Embedding Services:
    - AI embedding generation latency
    - AI provider request distribution
    - Embedding cache hit rate
    
  Memory Optimization:
    - Memory optimization effectiveness (%)
    - Memory usage per 100 variants
    - Cache performance metrics
    
  System Resources:
    - CPU utilization percentage
    - Memory usage and available memory
    - Disk I/O and network traffic
```

### Dashboard Panels

#### 1. Request Overview Panel
- Total requests per minute
- Success rate percentage
- Error rate trends
- Response time distribution

#### 2. VCF Processing Panel
- Variants processed per second
- File upload and processing times
- Processing queue length
- Throughput trends

#### 3. AI Services Panel
- Embedding generation performance
- AI provider usage distribution
- AI request latency breakdown
- Model performance metrics

#### 4. Memory Optimization Panel
- Memory reduction percentage
- Cache hit rate trends
- Memory usage efficiency
- Optimization impact metrics

#### 5. System Health Panel
- CPU and memory utilization
- Container health status
- Service availability
- Resource alerts

## Alert Rules (Prometheus)

### Critical Alerts
```yaml
Critical Alerts:
  Service Down:
    condition: up == 0
    duration: 30s
    severity: critical
    
  High Error Rate:
    condition: rate(http_requests_total{status=~"5.."}[5m]) > 0.10
    duration: 2m
    severity: critical
    
  Memory Optimization Degraded:
    condition: vcf_agent_memory_optimization_reduction_percent < 40
    duration: 5m
    severity: critical
    
  VCF Processing Stopped:
    condition: rate(vcf_variants_processed_total[5m]) == 0
    duration: 3m
    severity: critical
```

### Warning Alerts
```yaml
Warning Alerts:
  High Latency:
    condition: histogram_quantile(0.95, rate(http_request_duration_seconds_bucket[5m])) > 2.0
    duration: 5m
    severity: warning
    
  High CPU Usage:
    condition: avg(cpu_usage_percent) > 80
    duration: 10m
    severity: warning
    
  High Memory Usage:
    condition: memory_usage_percent > 85
    duration: 5m
    severity: warning
    
  Low Cache Hit Rate:
    condition: vcf_agent_cache_hit_rate < 70
    duration: 10m
    severity: warning
```

### Info Alerts
```yaml
Info Alerts:
  High Request Rate:
    condition: rate(http_requests_total[5m]) > 100
    duration: 5m
    severity: info
    
  Memory Optimization Active:
    condition: vcf_agent_memory_optimization_enabled == 1
    duration: 1m
    severity: info
```

## Distributed Tracing

### Jaeger Configuration

The VCF Agent includes comprehensive distributed tracing for debugging and performance analysis:

```yaml
Tracing Coverage:
  HTTP Requests: All API endpoints traced
  VCF Processing: File upload to completion
  AI Services: Embedding generation and caching
  Database Operations: LanceDB and Kuzu queries
  Memory Operations: Optimization and cleanup
```

### Trace Spans
- **HTTP Request Span**: Complete request lifecycle
- **VCF Processing Span**: File parsing and validation
- **Embedding Generation Span**: AI service interactions
- **Database Query Span**: Vector and graph operations
- **Memory Optimization Span**: Cache and cleanup operations

### Performance Analysis
```yaml
Key Trace Metrics:
  Total Request Time: End-to-end latency
  Service Dependencies: Inter-service call patterns
  Bottleneck Identification: Slow operations detection
  Error Propagation: Failure analysis across services
```

## Production Security Features

### Container Security
```yaml
Security Hardening:
  Container Execution:
    user: appuser (UID: 1000, non-root)
    capabilities: minimal (dropped ALL, added NET_BIND_SERVICE)
    
  Filesystem Security:
    root_filesystem: read-only
    temp_mounts: /tmp, /var/tmp (writable)
    secrets_mount: /app/secrets (read-only, 600 permissions)
    
  Network Isolation:
    app_network: application services
    observability_network: monitoring stack
    external_access: controlled ports only
    
  TLS Configuration:
    certificates: production-ready setup
    encryption: end-to-end encryption capable
    protocols: TLS 1.2+ only
```

### Secrets Management
```yaml
Secrets Configuration:
  Storage: External files (not in container images)
  Permissions: 600 (owner read/write only)
  Rotation: Automated rotation capability
  Access: Principle of least privilege
```

## Metrics Collection

### Prometheus Metrics

#### VCF Agent Specific Metrics
```yaml
Application Metrics:
  vcf_agent_requests_total: Total HTTP requests
  vcf_agent_request_duration_seconds: Request latency histogram
  vcf_agent_variants_processed_total: Variants processed counter
  vcf_agent_memory_optimization_reduction_percent: Memory reduction percentage
  vcf_agent_cache_hit_rate: Cache hit rate percentage
  vcf_agent_ai_requests_total: AI service requests by provider
  vcf_agent_embedding_generation_duration_seconds: Embedding latency
```

#### System Metrics
```yaml
System Metrics:
  process_cpu_usage_percent: CPU utilization
  process_memory_usage_bytes: Memory usage
  process_memory_available_bytes: Available memory
  disk_io_operations_total: Disk I/O operations
  network_bytes_transmitted_total: Network traffic
```

### Metric Labels
```yaml
Common Labels:
  service: vcf-agent
  environment: production|staging|development
  version: application version
  instance: container instance ID
  provider: ai_provider_name (for AI metrics)
  status_code: HTTP status codes
```

## Environment Configuration

### Production vs Development

#### Production Monitoring (10% Sampling)
```yaml
production:
  opentelemetry:
    sampling_rate: 0.1
    batch_timeout: 5s
    max_batch_size: 512
  prometheus:
    scrape_interval: 15s
    retention: 30d
  grafana:
    refresh_interval: 30s
```

#### Development Monitoring (100% Sampling)
```yaml
development:
  opentelemetry:
    sampling_rate: 1.0
    batch_timeout: 1s
    max_batch_size: 128
  prometheus:
    scrape_interval: 5s
    retention: 7d
  grafana:
    refresh_interval: 5s
```

## Health Checks

### Service Health Endpoints
```yaml
Health Check Endpoints:
  VCF Agent: GET /health
  Response: {"status": "healthy", "checks": {...}}
  Timeout: 2 seconds
  
  Prometheus: GET /api/v1/status/config
  Grafana: GET /api/health
  Jaeger: GET /
```

### Health Check Content
```json
{
  "status": "healthy",
  "timestamp": "2025-05-29T12:36:00Z",
  "checks": {
    "database": "healthy",
    "memory_optimization": "active",
    "cache": "operational",
    "ai_services": "available"
  },
  "metrics": {
    "memory_reduction_percent": 95.2,
    "cache_hit_rate": 87.3,
    "variants_per_second": 28.1
  }
}
```

## Monitoring Best Practices

### Dashboard Organization
1. **Executive Summary**: High-level KPIs and status
2. **Operational View**: Day-to-day monitoring metrics
3. **Technical Deep Dive**: Detailed performance analysis
4. **Troubleshooting**: Error analysis and debugging

### Alert Management
1. **Alert Fatigue Prevention**: Tuned thresholds to reduce noise
2. **Escalation Paths**: Clear escalation procedures
3. **Runbook Integration**: Links to troubleshooting procedures
4. **Alert Grouping**: Related alerts grouped together

### Performance Optimization
1. **Metric Retention**: Balanced retention vs. storage costs
2. **Sampling Strategies**: Environment-specific sampling rates
3. **Dashboard Performance**: Optimized queries and refresh rates
4. **Resource Allocation**: Right-sized monitoring infrastructure

## Troubleshooting Guide

### Common Issues

#### 1. High Memory Usage Alert
```bash
# Check memory optimization status
curl http://localhost:8080/health | jq '.metrics.memory_reduction_percent'

# Force memory cleanup
kubectl exec -it vcf-agent-pod -- curl -X POST http://localhost:8080/admin/cleanup
```

#### 2. Low Cache Hit Rate
```bash
# Check cache statistics
curl http://localhost:8080/metrics | grep cache_hit_rate

# Clear and rebuild cache
kubectl exec -it vcf-agent-pod -- curl -X POST http://localhost:8080/admin/cache/reset
```

#### 3. High Error Rate
```bash
# Check error logs
kubectl logs vcf-agent-pod --tail=100 | grep ERROR

# Check error distribution
curl http://localhost:9090/api/v1/query?query='rate(http_requests_total{status=~"5.."}[5m])'
```

#### 4. Service Discovery Issues
```bash
# Verify Prometheus targets
curl http://localhost:9090/api/v1/targets

# Check service health
curl http://localhost:8080/health
curl http://localhost:9090/-/healthy
curl http://localhost:3000/api/health
```

## Integration with CI/CD

### Automated Monitoring Setup
```yaml
Deployment Pipeline Integration:
  1. Deploy monitoring stack
  2. Validate service discovery
  3. Execute health checks
  4. Verify metric collection
  5. Test alert firing
  6. Validate dashboard rendering
```

### Monitoring Validation
```bash
# Automated monitoring validation script
scripts/validate-monitoring.sh
```

## Links and References

- [Production Deployment Guide](deployment/production-deployment-runbook.md)
- [Grafana Dashboard Configuration](../config/grafana/dashboards/)
- [Prometheus Alert Rules](../config/prometheus/rules/)
- [OpenTelemetry Configuration](../config/otel/)

---

**Next Steps**: After monitoring setup, configure alert notifications and establish operational procedures. 