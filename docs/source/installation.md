# Installation Guide

This guide covers the installation and setup of the VCF Analysis Agent for different environments and use cases.

## Quick Installation

### Prerequisites

- **Python**: 3.9 or higher
- **Docker**: For containerized deployment (recommended)
- **Git**: For source code access
- **System Memory**: Minimum 8GB RAM (16GB+ recommended for large VCF files)
- **Storage**: 10GB+ free space for databases and cache

### Option 1: Docker Installation (Recommended)

The fastest way to get started is using Docker:

```bash
# Clone the repository
git clone https://github.com/your-org/vcf-analysis-agent.git
cd vcf-analysis-agent

# Start with Docker Compose
docker-compose up -d

# Verify installation
docker-compose exec vcf-agent vcf-agent --version
```

### Option 2: Python Package Installation

```bash
# Install from PyPI (when available)
pip install vcf-analysis-agent

# Or install from source
git clone https://github.com/your-org/vcf-analysis-agent.git
cd vcf-analysis-agent
pip install -e .
```

### Option 3: Development Installation

For development and customization:

```bash
# Clone repository
git clone https://github.com/your-org/vcf-analysis-agent.git
cd vcf-analysis-agent

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
pip install -e .

# Install development dependencies
pip install -r requirements-dev.txt
```

## System Requirements

### Minimum Requirements

- **CPU**: 2 cores
- **RAM**: 8GB
- **Storage**: 10GB free space
- **OS**: Linux, macOS, or Windows with WSL2

### Recommended Requirements

- **CPU**: 4+ cores
- **RAM**: 16GB+
- **Storage**: 50GB+ SSD
- **OS**: Linux (Ubuntu 20.04+) or macOS

### For Large-Scale Processing

- **CPU**: 8+ cores
- **RAM**: 32GB+
- **Storage**: 100GB+ NVMe SSD
- **Network**: High-bandwidth for AI API calls

## Dependencies

### Core Dependencies

- **Python Libraries**: pandas, numpy, pydantic, typer
- **Database Systems**: LanceDB, Kuzu
- **AI Integration**: openai, litellm, ollama
- **Monitoring**: prometheus-client, opentelemetry

### External Tools

- **bcftools**: VCF file processing (auto-installed in Docker)
- **Ollama**: Local AI models (optional)

### AI Model Providers

Choose one or more:

1. **OpenAI**: Requires API key
2. **Cerebras**: Requires API key  
3. **Ollama**: Local models (no API key needed)

## Configuration

### Environment Variables

Create a `.env` file in the project root:

```bash
# AI Provider Configuration
OPENAI_API_KEY=your_openai_key_here
CEREBRAS_API_KEY=your_cerebras_key_here

# Database Configuration
LANCEDB_PATH=./lancedb
KUZU_DB_PATH=./kuzu_db

# Monitoring (Optional)
VCF_AGENT_METRICS_PORT=8000
VCF_AGENT_PUSHGATEWAY_URL=http://localhost:9091

# Logging
LOG_LEVEL=INFO
```

### Database Setup

The databases are automatically initialized on first run:

```bash
# Initialize databases
vcf-agent init

# Verify database connections
vcf-agent status
```

## Verification

### Test Installation

```bash
# Check version
vcf-agent --version

# Run basic validation
vcf-agent validate sample_data/example.vcf

# Test AI integration
vcf-agent analyze sample_data/example.vcf --provider ollama
```

### Performance Test

```bash
# Run performance validation
python scripts/performance_validation.py
```

Expected output:
```