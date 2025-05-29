# VCF Agent Production Deployment Runbook

## Table of Contents
1. [Pre-Deployment Checklist](#pre-deployment-checklist)
2. [Deployment Procedures](#deployment-procedures)
3. [Post-Deployment Validation](#post-deployment-validation)
4. [Troubleshooting Guide](#troubleshooting-guide)
5. [Scaling Procedures](#scaling-procedures)
6. [Rollback Procedures](#rollback-procedures)
7. [Monitoring and Alerting](#monitoring-and-alerting)

## Pre-Deployment Checklist

### Infrastructure Requirements
- [ ] Docker Engine 20.10+ installed
- [ ] Docker Compose 2.0+ installed
- [ ] Minimum 4GB RAM available
- [ ] Minimum 2 CPU cores available
- [ ] 50GB+ disk space available
- [ ] Network connectivity to AI providers (OpenAI, Anthropic)
- [ ] SSL certificates available (production only)

### Security Requirements
- [ ] API keys securely stored in secrets management
- [ ] TLS certificates valid and properly configured
- [ ] Network security groups configured
- [ ] Firewall rules applied
- [ ] Secret rotation schedule established

### Environment Configuration
- [ ] Environment variables configured
- [ ] Configuration files validated
- [ ] Secrets files created and secured
- [ ] Backup storage configured
- [ ] Monitoring endpoints accessible

## Deployment Procedures

### 1. Initial Production Deployment

#### Step 1: Prepare Environment
```bash
# Create deployment directory
mkdir -p /opt/vcf-agent
cd /opt/vcf-agent

# Clone repository
git clone https://github.com/your-org/vcf-agent.git .
git checkout v1.0.0  # Use specific version tag

# Create secrets directory
mkdir -p secrets
chmod 700 secrets
```

#### Step 2: Configure Secrets
```bash
# Create API key files
echo "your-openai-api-key" > secrets/openai_api_key.txt
echo "your-anthropic-api-key" > secrets/anthropic_api_key.txt
echo "your-jwt-secret" > secrets/jwt_secret.txt

# Secure permissions
chmod 600 secrets/*.txt
chown root:root secrets/*.txt
```

#### Step 3: Configure Environment
```bash
# Create production environment file
cat > .env.production << EOF
ENVIRONMENT=production
VERSION=v1.0.0
LOG_LEVEL=INFO
OTEL_SAMPLING_RATE=0.1
GRAFANA_ADMIN_PASSWORD=secure-password-here
OPENAI_API_KEY_FILE=./secrets/openai_api_key.txt
ANTHROPIC_API_KEY_FILE=./secrets/anthropic_api_key.txt
OTEL_CREDENTIALS_FILE=./secrets/otel_credentials.txt
EOF
```

#### Step 4: Deploy Services
```bash
# Build and start services
docker-compose -f docker-compose.production.yml --env-file .env.production up -d

# Verify all services are running
docker-compose -f docker-compose.production.yml ps
```

#### Step 5: Validate Deployment
```bash
# Check service health
curl -f http://localhost:8080/health
curl -f http://localhost:9090/-/healthy  # Prometheus
curl -f http://localhost:3000/api/health  # Grafana
curl -f http://localhost:16686/  # Jaeger

# Check logs for errors
docker-compose -f docker-compose.production.yml logs vcf-agent
```

### 2. Rolling Updates

#### Step 1: Prepare Update
```bash
# Pull latest changes
git fetch origin
git checkout v1.1.0  # New version

# Build new image
docker-compose -f docker-compose.production.yml build vcf-agent
```

#### Step 2: Rolling Update
```bash
# Update with zero downtime
docker-compose -f docker-compose.production.yml up -d --no-deps vcf-agent

# Verify new version is running
docker-compose -f docker-compose.production.yml exec vcf-agent python -c "import vcf_agent; print(vcf_agent.__version__)"
```

#### Step 3: Validate Update
```bash
# Run health checks
./scripts/health-check.sh

# Monitor metrics for 10 minutes
./scripts/monitor-deployment.sh
```

## Post-Deployment Validation

### Health Check Script
```bash
#!/bin/bash
# scripts/health-check.sh

set -e

echo "Running post-deployment health checks..."

# Service health checks
services=("vcf-agent:8080" "prometheus:9090" "grafana:3000" "jaeger:16686")
for service in "${services[@]}"; do
    IFS=':' read -r name port <<< "$service"
    echo "Checking $name on port $port..."
    curl -f "http://localhost:$port/health" || curl -f "http://localhost:$port/" || {
        echo "ERROR: $name health check failed"
        exit 1
    }
done

# VCF Agent specific checks
echo "Testing VCF Agent functionality..."
curl -X POST http://localhost:8080/api/v1/validate \
    -H "Content-Type: application/json" \
    -d '{"vcf_content": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"}' || {
    echo "ERROR: VCF validation endpoint failed"
    exit 1
}

# Memory optimization check
memory_reduction=$(curl -s http://localhost:8080/metrics | grep vcf_agent_memory_optimization_reduction_percent | awk '{print $2}')
if (( $(echo "$memory_reduction < 40" | bc -l) )); then
    echo "WARNING: Memory optimization below 40% ($memory_reduction%)"
fi

echo "All health checks passed!"
```

### Performance Validation
```bash
#!/bin/bash
# scripts/performance-check.sh

echo "Running performance validation..."

# Check response times
response_time=$(curl -w "%{time_total}" -s -o /dev/null http://localhost:8080/health)
if (( $(echo "$response_time > 2.0" | bc -l) )); then
    echo "WARNING: Health endpoint response time high: ${response_time}s"
fi

# Check memory usage
memory_usage=$(docker stats --no-stream --format "table {{.Container}}\t{{.MemUsage}}" | grep vcf-agent)
echo "Memory usage: $memory_usage"

# Check error rate
error_rate=$(curl -s http://localhost:9090/api/v1/query?query=rate\(vcf_agent_http_requests_total\{status\=~\"5..\"\}\[5m\]\) | jq -r '.data.result[0].value[1] // "0"')
if (( $(echo "$error_rate > 0.05" | bc -l) )); then
    echo "WARNING: Error rate high: ${error_rate}"
fi

echo "Performance validation complete"
```

## Troubleshooting Guide

### Common Issues

#### 1. High Memory Usage
**Symptoms:** Memory usage >85%, OOM kills
**Diagnosis:**
```bash
# Check memory metrics
curl -s http://localhost:9090/api/v1/query?query=system_memory_utilization

# Check memory optimization effectiveness
curl -s http://localhost:8080/metrics | grep memory_optimization_reduction_percent
```
**Resolution:**
1. Verify memory optimization is enabled
2. Check embedding dimensions configuration
3. Restart service if memory leak suspected
4. Scale horizontally if needed

#### 2. Slow Embedding Generation
**Symptoms:** P95 latency >2 seconds
**Diagnosis:**
```bash
# Check embedding latency metrics
curl -s http://localhost:9090/api/v1/query?query=histogram_quantile\(0.95,rate\(vcf_agent_embedding_generation_duration_seconds_bucket\[5m\]\)\)

# Check AI provider status
curl -s http://localhost:8080/debug/ai-providers
```
**Resolution:**
1. Check AI provider API status
2. Verify API rate limits not exceeded
3. Consider switching to faster embedding model
4. Implement request batching

#### 3. Tracing Data Loss
**Symptoms:** Missing traces in Jaeger
**Diagnosis:**
```bash
# Check OpenTelemetry Collector status
docker-compose logs otel-collector

# Check export success rate
curl -s http://localhost:9090/api/v1/query?query=rate\(otelcol_exporter_sent_spans_total\[5m\]\)
```
**Resolution:**
1. Restart OpenTelemetry Collector
2. Check Jaeger connectivity
3. Verify sampling configuration
4. Check memory limits on collector

#### 4. Service Won't Start
**Symptoms:** Container exits immediately
**Diagnosis:**
```bash
# Check container logs
docker-compose logs vcf-agent

# Check configuration
docker-compose config

# Verify secrets
ls -la secrets/
```
**Resolution:**
1. Verify all required secrets exist
2. Check configuration file syntax
3. Ensure proper file permissions
4. Verify Docker image integrity

### Emergency Procedures

#### Immediate Service Recovery
```bash
# Quick restart
docker-compose -f docker-compose.production.yml restart vcf-agent

# Full stack restart
docker-compose -f docker-compose.production.yml down
docker-compose -f docker-compose.production.yml up -d

# Emergency rollback
git checkout v1.0.0  # Previous working version
docker-compose -f docker-compose.production.yml up -d --force-recreate
```

## Scaling Procedures

### Horizontal Scaling
```bash
# Scale VCF Agent instances
docker-compose -f docker-compose.production.yml up -d --scale vcf-agent=3

# Verify load distribution
curl -s http://localhost:9090/api/v1/query?query=up\{job\=\"vcf-agent\"\}
```

### Vertical Scaling
```yaml
# Update docker-compose.production.yml
services:
  vcf-agent:
    deploy:
      resources:
        limits:
          memory: 4G
          cpus: '2.0'
        reservations:
          memory: 2G
          cpus: '1.0'
```

### Auto-scaling Configuration
```yaml
# For Kubernetes deployment
apiVersion: autoscaling/v2
kind: HorizontalPodAutoscaler
metadata:
  name: vcf-agent-hpa
spec:
  scaleTargetRef:
    apiVersion: apps/v1
    kind: Deployment
    name: vcf-agent
  minReplicas: 2
  maxReplicas: 10
  metrics:
  - type: Resource
    resource:
      name: cpu
      target:
        type: Utilization
        averageUtilization: 70
  - type: Resource
    resource:
      name: memory
      target:
        type: Utilization
        averageUtilization: 80
```

## Rollback Procedures

### Automated Rollback
```bash
#!/bin/bash
# scripts/rollback.sh

PREVIOUS_VERSION=${1:-v1.0.0}

echo "Rolling back to version $PREVIOUS_VERSION..."

# Checkout previous version
git checkout $PREVIOUS_VERSION

# Rebuild and deploy
docker-compose -f docker-compose.production.yml build vcf-agent
docker-compose -f docker-compose.production.yml up -d --no-deps vcf-agent

# Validate rollback
sleep 30
./scripts/health-check.sh

echo "Rollback to $PREVIOUS_VERSION complete"
```

### Manual Rollback Steps
1. Identify last known good version
2. Stop current services
3. Checkout previous version
4. Rebuild containers
5. Start services
6. Validate functionality
7. Update monitoring dashboards

## Monitoring and Alerting

### Key Metrics to Monitor
- **Error Rate:** <5% (Critical: >10%)
- **Response Time:** P95 <2s (Warning: >2s)
- **Memory Usage:** <85% (Warning: >85%)
- **CPU Usage:** <80% (Warning: >80%)
- **Memory Optimization:** >40% (Critical: <40%)
- **Throughput:** >100 variants/sec (Warning: <100)

### Alert Escalation
1. **Info:** Log to monitoring system
2. **Warning:** Notify on-call engineer
3. **Critical:** Page on-call engineer + manager

### Dashboard URLs
- **Grafana:** http://localhost:3000/d/vcf-agent-overview
- **Prometheus:** http://localhost:9090
- **Jaeger:** http://localhost:16686

### Log Locations
- **Application Logs:** `docker-compose logs vcf-agent`
- **System Logs:** `/var/log/syslog`
- **Audit Logs:** `/var/log/audit/audit.log`

## Maintenance Procedures

### Daily Tasks
- [ ] Check service health status
- [ ] Review error logs
- [ ] Verify backup completion
- [ ] Monitor resource usage

### Weekly Tasks
- [ ] Review performance metrics
- [ ] Update security patches
- [ ] Rotate secrets if needed
- [ ] Clean up old logs

### Monthly Tasks
- [ ] Review and update runbooks
- [ ] Conduct disaster recovery test
- [ ] Update monitoring thresholds
- [ ] Security audit

## Contact Information

### On-Call Escalation
1. **Primary:** DevOps Engineer (Slack: @devops-oncall)
2. **Secondary:** Platform Team Lead (Phone: +1-555-0123)
3. **Escalation:** Engineering Manager (Email: manager@company.com)

### External Dependencies
- **OpenAI Support:** https://help.openai.com
- **Anthropic Support:** https://support.anthropic.com
- **Infrastructure Provider:** AWS Support (Case Priority: High)

---

**Document Version:** 1.0  
**Last Updated:** 2025-01-05  
**Next Review:** 2025-02-05 