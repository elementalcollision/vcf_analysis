# VCF Analysis Agent

A powerful and intelligent agent for analyzing, processing, and managing Variant Call Format (VCF) files in genomics research and clinical applications.

## Overview

The VCF Analysis Agent is an AI-powered tool that combines the robustness of bcftools with modern AI capabilities to provide intelligent analysis, validation, and processing of VCF files. It aims to streamline genomic variant analysis workflows while ensuring data quality and compliance with standards.

## Features

- Intelligent VCF file processing and validation
- AI-powered variant analysis and interpretation
- Integration with bcftools for core VCF operations
- Automated quality control and standardization
- Smart filtering and annotation capabilities
- Extensible plugin architecture for custom analyses

## LLM Provider Integration & Credential Management

The VCF Analysis Agent supports multiple Large Language Model (LLM) providers for flexible, secure, and cost-effective AI integration:

- **Ollama** (local, open-source; default)
- **OpenAI** (cloud, commercial)
- **Cerebras** (cloud, specialized; mock implementation for now)

### Credential Management

You can provide API credentials in two ways:

1. **Environment Variables (.env file)**
   - Place a `.env` file in your project root with:
     ```env
     OPENAI_API_KEY=sk-...
     CEREBRAS_API_KEY=csk-...
     ```
   - These will be loaded automatically at runtime.

2. **JSON Credentials File**
   - Create a file (e.g., `my_credentials.json`):
     ```json
     {
       "openai": { "api_key": "sk-..." },
       "cerebras": { "api_key": "csk-..." }
     }
     ```
   - Pass the path via CLI (`--credentials my_credentials.json`) or in your config/session.
   - The JSON file takes precedence over environment variables if both are present.

### Selecting LLM Providers

You can select the LLM provider in several ways:

- **CLI:**
  ```bash
  python -m vcf_agent.cli --model openai --credentials my_credentials.json "echo: Hello from OpenAI!"
  python -m vcf_agent.cli --model cerebras "echo: Hello from Cerebras!"
  python -m vcf_agent.cli --model ollama "echo: Hello from Ollama!"
  ```
- **SessionConfig (Python):**
  ```python
  from vcf_agent.config import SessionConfig
  from vcf_agent.agent import get_agent_with_session
  session_config = SessionConfig(model_provider="openai", credentials_file="my_credentials.json")
  agent = get_agent_with_session(session_config)
  ```
- **Default:** If not specified, the agent uses Ollama (local) by default.

### Architecture Notes
- The agent uses a unified API client interface for all providers (see `src/vcf_agent/api_clients.py`).
- Provider selection and credential loading are fully documented in the [API Integration Decision Record](.context/decisions/2025-05-21-api-integration-update.md).
- Cerebras integration is currently a mock for testing and future extension.

## Installation

```bash
# Clone the repository
git clone https://github.com/elementalcollision/vcf_analysis.git
cd vcf_analysis

# Set up Python environment (Python 3.11+ required)
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install uv (modern Python dependency manager)
pip install uv

# Install dependencies using uv and pyproject.toml
uv pip install -r requirements.txt
```

- To add or update dependencies, edit `pyproject.toml` and run:
  ```bash
  uv pip compile pyproject.toml --output-file=requirements.txt
  uv pip install -r requirements.txt
  ```

## Containerization (Docker & OrbStack)

The agent is fully containerized for production, local development, and CI/CD, following best practices for security, reproducibility, and multi-architecture support.

### Dockerfile Highlights
- **Multi-stage build:** Separates build and runtime for smaller, more secure images.
- **Pinned Python version:** Uses `python:3.11.8-slim-bookworm` for reproducibility and compatibility.
- **Non-root user:** Runs as `appuser` for security.
- **Layer caching:** Dependency install steps are optimized for fast rebuilds.
- **.dockerignore:** Ensures only necessary files are included in the image.

### Multi-Architecture Builds (Buildx)
- Supports both `linux/amd64` (x86) and `linux/arm64` (Apple Silicon, OrbStack, cloud ARM).
- Use Docker Buildx for cross-platform builds:
  ```bash
  docker buildx create --use
  docker buildx build --platform linux/amd64,linux/arm64 -t vcf-agent:latest .
  # For local dev (single arch):
  docker build -t vcf-agent:dev .
  ```
- OrbStack supports both ARM and x86 images; test with your local architecture as needed.

### Example: Build and Run Locally
```bash
# Build the image (single arch)
docker build -t vcf-agent:dev .
# Run the agent CLI
docker run --rm -it vcf-agent:dev "echo: Hello from Docker!"
```

### Example: Multi-Arch Build for CI/CD
```bash
docker buildx create --use
docker buildx build --platform linux/amd64,linux/arm64 -t vcf-agent:latest .
```

### OrbStack Compatibility
- OrbStack uses BuildKit and supports multi-arch images.
- Use volume mounts for fast local iteration.
- Images are accessible at `~/OrbStack/docker/images` for debugging.

### Security & Best Practices
- Runs as a non-root user.
- Only production dependencies and code are included in the runtime image.
- Use Trivy or Snyk in CI/CD to scan for vulnerabilities.
- Pin all dependencies in `requirements.txt` and `pyproject.toml`.

### Troubleshooting
- If you encounter architecture errors, ensure Buildx is enabled and the correct platform is specified.
- For missing dependencies, add them to the builder stage in the Dockerfile.
- Use `docker builder prune -a` to clear build cache if disk space is low.

### CI/CD Integration
- Use Buildx in your CI pipeline for reproducible, multi-arch builds.
- Add security scanning as a pipeline step.

For more details, see the Dockerfile and `.dockerignore` in the repo, and consult the [OrbStack Docker documentation](https://docs.orbstack.dev/docker/images).

## Development & Testing

The project maintains a comprehensive test suite to ensure code quality, stability, and correctness.

### Test Suite Structure

The test suite is organized into the following categories:

-   **Unit Tests (`tests/unit/`)**: Focus on individual modules and functions in isolation. Extensive mocking is used for external dependencies like KuzuDB, LanceDB, and OpenTelemetry.
-   **Integration Tests (`tests/integration/`)**: Test the interaction between different components of the agent, including CLI, API endpoints, and database integrations.
    -   **End-to-End (E2E) Tests (`tests/integration/e2e/`)**: Validate complete workflows from the user's perspective, such as CLI-based data ingestion and querying.
-   **Prompt Contract Tests (`tests/prompt_contracts/`)**: Verify the agent's interactions with LLMs for specific tasks, ensuring consistent output formats.
-   **Agent Tests (`tests/test_agent.py`, `tests/test_agent_llm_integration.py`)**: Test core agent logic, tool dispatch, and LLM integration.
-   **Tool & Utility Tests**: Dedicated tests for bcftools integration, VCF utilities, API clients, validation logic, etc.

### Running Tests

To run the full test suite:

```bash
pytest -v
```

To run tests with coverage for the `src/vcf_agent` directory and report missing lines:

```bash
pytest --cov=src/vcf_agent --cov-report=term-missing -v
```

HTML and XML coverage reports are also generated (`htmlcov/` and `coverage.xml`).

### Current Test Status (as of 2025-08-01)

-   **Total Tests**: ~258
-   **Passing**: ~231
-   **Skipped**: ~16
    -   ~5 due to unresolved OpenAI schema validation issues for array tool parameters (see "Known Issues").
    -   ~10 in `test_edgecase_compliance.py` due to missing VCF test files (see "Known Issues").
    -   ~1 performance test in `test_vcf_comparison_tool.py` due to a missing large VCF file (see "Known Issues").
-   **XFAIL (Expected Failures)**: 5 (related to `bcftools stats` permissiveness on certain VCFs).
-   **XPASS (Unexpected Passes)**: 5 (GATK validation tests for the same VCFs now pass, indicating stricter validation than `bcftools stats`).

### Coverage

-   **Overall Project Coverage**: ~73%
-   **`src/vcf_agent/tracing.py`**: ~95% (Remaining "missed" lines are typically comments or outdated references from the coverage tool for the current code state).
-   **`src/vcf_agent/metrics.py`**: ~76% (Remaining "missed" lines in original targets often correspond to helper functions that are not currently implemented. Existing functions are well-covered).

Significant effort (documented in `TASK-007.md`) was undertaken to resolve critical test failures and improve coverage for tracing and metrics.

### Known Issues & Future Testing Considerations

1.  **OpenAI Schema Validation for Array Types**: Several tests are skipped due to `litellm.exceptions.BadRequestError` related to array tool schemas for OpenAI. This requires further investigation, potentially with Strands/LiteLLM.
2.  **Missing Edge Case VCF Files**: Some tests in `test_edgecase_compliance.py` are skipped as they require specific VCF files not currently in `sample_test_data/`. These files need to be created or sourced.
3.  **Large VCF File for Performance Test**: The performance test `tests/test_vcf_comparison_tool.py::test_vcf_comparison_large_file_performance` requires a large VCF file (~129MB) and a corresponding reference FASTA. Due to its size, the VCF is not stored in the repository. Instructions for downloading the base 1000 Genomes chr22 VCF, and then preparing a `testready` version (e.g., `chr22.1kg.phase3.v5a.testready.vcf.gz`) with the correct contig header for use with a `22.fa` reference, are detailed in `sample_test_data/README.md`. This preparation is necessary because the `bcftools` version 1.21 on the primary test system (macOS via Homebrew) did not support the `--add-missing-contigs` option for `bcftools norm`.
4.  **Helper Function Coverage in `metrics.py`**: While core metrics definitions and existing helper functions are tested, some originally planned observer helper functions (e.g., for CLI commands, specific tool usage beyond AI) are not currently implemented in `metrics.py`. If these are added in the future, corresponding tests will be needed.

Consult `.context/tasks/completed/TASK-007.md` for a detailed history of recent test stabilization efforts.

## Observability Stack

The VCF Analysis Agent includes a comprehensive observability stack that provides monitoring, metrics collection, distributed tracing, and visualization capabilities. This stack consists of:

- **OpenTelemetry (OTel)**: Distributed tracing and metrics collection
- **Prometheus**: Metrics storage and querying  
- **Grafana**: Metrics visualization and dashboards
- **Jaeger**: Distributed trace visualization and analysis

### Quick Start

To run the complete observability stack locally:

```bash
# Start all services (Prometheus, Grafana, Jaeger)
docker-compose up -d

# Run the VCF agent (metrics will be automatically collected)
python -m vcf_agent.cli ask "What are the basic stats for sample_data/minimal.vcf.gz?"

# Access the UIs:
# - Grafana: http://localhost:3000 (admin/admin)
# - Prometheus: http://localhost:9090
# - Jaeger: http://localhost:16686
# - Agent Metrics: http://localhost:8000/metrics
```

### Key Features

- **Automatic Instrumentation**: OpenTelemetry automatically traces HTTP requests, logging, and asyncio operations
- **Custom Metrics**: Detailed metrics for AI interactions, tool usage, BCFtools operations, and CLI commands
- **Structured Logging**: JSON logs with trace context for correlation
- **Real-time Dashboards**: Pre-built Grafana dashboards for monitoring AI performance
- **Distributed Tracing**: End-to-end request tracing through Jaeger

### Metrics Collection

The agent exposes metrics in two ways:

1. **HTTP Endpoint**: Real-time metrics at `http://localhost:8000/metrics` for Prometheus scraping
2. **Pushgateway**: For short-lived CLI commands, metrics are pushed to Prometheus Pushgateway

**Key Metrics Categories:**
- **AI Interactions**: Request rates, response times, error rates, token usage
- **Tool Usage**: Performance and error tracking for all agent tools
- **BCFtools Operations**: Subprocess execution metrics
- **CLI Commands**: Command execution duration and success rates

**Environment Variables:**
- `VCF_AGENT_METRICS_PORT`: Metrics HTTP server port (default: 8000)
- `VCF_AGENT_PUSHGATEWAY_URL`: Pushgateway URL for CLI metrics

### Distributed Tracing

OpenTelemetry provides end-to-end tracing of requests through the system:

- **Automatic Spans**: HTTP requests, database operations, external API calls
- **Custom Spans**: Tool executions, AI model interactions, VCF processing
- **Trace Context**: Propagated through logs for correlation

**Environment Variables:**
- `OTEL_EXPORTER_OTLP_TRACES_ENDPOINT`: Jaeger collector endpoint (default: http://localhost:4317)
- `OTEL_EXPORTER_OTLP_TRACES_PROTOCOL`: Protocol (grpc or http/protobuf)
- `OTEL_PYTHON_LOG_LEVEL`: Set to "debug" for detailed OTel logging

### Developer Guide

#### Adding Custom Metrics

```python
from vcf_agent.metrics import http_registry
from prometheus_client import Counter

# Define a new metric
my_custom_counter = Counter(
    'my_custom_operations_total',
    'Description of my custom metric',
    ['operation_type', 'status'],
    registry=http_registry
)

# Use the metric
my_custom_counter.labels(operation_type='validation', status='success').inc()
```

#### Adding Custom Tracing

```python
from vcf_agent.tracing import init_tracer

# Get a tracer instance
tracer = init_tracer("my-service")

# Create custom spans
with tracer.start_as_current_span("my_operation") as span:
    span.set_attribute("operation.type", "validation")
    span.set_attribute("file.path", "/path/to/file.vcf")
    # Your operation code here
    span.set_attribute("operation.result", "success")
```

#### Observing AI Interactions

The agent automatically tracks AI interactions, but you can also manually record them:

```python
from vcf_agent.metrics import observe_ai_interaction

# Record an AI interaction
observe_ai_interaction(
    model_provider="ollama",
    endpoint_task="vcf_analysis", 
    duration_seconds=2.5,
    prompt_tokens=150,
    completion_tokens=75,
    success=True
)
```

### Configuration Files

The observability stack is configured through several files:

- **`docker-compose.yml`**: Defines all observability services
- **`prometheus.yml`**: Prometheus scraping configuration
- **`config/grafana/dashboards/`**: Pre-built Grafana dashboards

### Troubleshooting

**No metrics in Prometheus:**
- Check that the agent metrics server is running on port 8000
- Verify Prometheus can reach `host.docker.internal:8000/metrics`
- Check Docker networking configuration

**No traces in Jaeger:**
- Verify Jaeger collector is running on port 4317
- Check `OTEL_EXPORTER_OTLP_TRACES_ENDPOINT` environment variable
- Enable debug logging with `OTEL_PYTHON_LOG_LEVEL=debug`

**Grafana dashboard issues:**
- Ensure Prometheus data source is configured correctly
- Import dashboard JSON from `config/grafana/dashboards/`
- Check that metrics are being generated by running agent commands

For detailed information on configuring and using the observability stack, see the [complete observability documentation](docs/source/monitoring_with_prometheus.md).

## Usage

- The agent and its tools can be used via the CLI by providing a prompt string. For example:

### Validate a VCF file
```bash
python -m vcf_agent.cli "validate_vcf: sample_data/<your_sample.vcf.gz>"
```
This will run the validation tool on the specified file and print the result (valid/invalid and reason).

> **Note:** Replace `<your_sample.vcf.gz>` with the actual filename of your sample VCF file in the `sample_data/` directory.

### Echo text (example tool)
```bash
python -m vcf_agent.cli "echo: Hello, world!"
```
This will return the echoed text from the agent.

You can invoke any registered tool by using the tool name followed by a colon and its argument(s).

## Data Flow Overview

A typical workflow with the VCF Analysis Agent involves the following steps:
1. **Ingestion**: VCF data is processed. Variant information and AI-generated embeddings are stored in LanceDB. Genomic relationships (e.g., variant-to-sample observations) are stored in Kuzu.
2. **AI Enrichment (Optional)**: Variants can be further enriched using AI models for annotation, clinical significance prediction, etc. (This is an area of ongoing development).
3. **Querying**: Users can query LanceDB for similar variants using embeddings or filter by metadata. Kuzu can be queried for contextual information about variants and samples.
4. **Analysis**: The combined data from LanceDB, Kuzu, and AI enrichment steps enables comprehensive variant analysis.

## LanceDB Integration and Usage Examples

The agent supports LanceDB as a vector database for storing and searching variant embeddings. You can use the following CLI commands to interact with LanceDB:

### Initialize LanceDB Table
```bash
python -m vcf_agent.cli init-lancedb --db_path ./lancedb --table_name variants
```

### Add a Variant Record
```bash
python -m vcf_agent.cli add-variant --variant_id rs123 --chrom 1 --pos 12345 --ref A --alt G --embedding 0.1,0.2,0.3,... --clinical_significance Pathogenic
```

### Search by Embedding
```bash
python -m vcf_agent.cli search-embedding --embedding 0.1,0.2,0.3,... --limit 5 --filter_sql "clinical_significance = 'Pathogenic'"
```

### Update a Variant Record
```bash
python -m vcf_agent.cli update-variant --variant_id rs123 --updates '{"clinical_significance": "Benign", "pos": 12346}'
```
- The `--updates` argument must be a valid JSON string representing the fields to change and their new values.

### Delete Variant Records
```bash
python -m vcf_agent.cli delete-variants --filter_sql "chrom = '1' AND pos < 10000 AND clinical_significance = 'Uncertain'"
```
- The `--filter_sql` argument uses SQL-like syntax to specify which records to delete.

### Create Scalar Index
```bash
python -m vcf_agent.cli create-lancedb-index --db_path ./lancedb --table_name variants --column chrom --index_type BTREE --replace
```
Creates a scalar index on a column (e.g., `chrom`, `pos`, `clinical_significance`) to improve filter performance.

### Filter by Metadata
```bash
python -m vcf_agent.cli filter-lancedb --db_path ./lancedb --table_name variants --filter_sql "chrom = '1' AND pos > 1000" --select_columns variant_id,chrom,pos --limit 10
```
Returns variants matching the metadata filter, with options to select specific columns and limit results.

- The `--embedding` argument (for `add-variant` and `search-embedding`) should be a comma-separated list of floats matching your embedding dimension (e.g., 1024).

See the [LanceDB Decision Record](.context/decisions/2025-05-24-lancedb-integration.md) for more details on schema and integration best practices.

## Kuzu Graph Database Integration

The agent also integrates with Kuzu, a high-performance graph database, to store and query relationships between variants, samples, and other genomic entities. This allows for more complex contextual queries.

### Kuzu Schema
The Kuzu database uses a schema involving `Variant` nodes, `Sample` nodes, and `ObservedIn` relationships:
- **Variant Node**: `variant_id` (primary key), `chrom`, `pos`, `ref`, `alt`, `rs_id`.
- **Sample Node**: `sample_id` (primary key).
- **ObservedIn Relationship**: Connects `Variant` to `Sample`, with a `zygosity` property.

### Kuzu CLI Commands

#### Initialize Kuzu Database and Schema
```bash
python -m vcf_agent.cli init-kuzu --db_path ./kuzu_db
```
This command creates the necessary node and relationship tables in the Kuzu database.

#### Populate Kuzu from VCF
```bash
python -m vcf_agent.cli populate-kuzu-from-vcf --vcf_file sample_data/your_sample.vcf.gz --kuzu_db_path ./kuzu_db
```
This command parses a VCF file and populates the Kuzu graph with `Variant` nodes, `Sample` nodes (derived from the VCF header), and `ObservedIn` relationships.

#### Get Variant Context from Kuzu
```bash
python -m vcf_agent.cli get-variant-context --variant_ids "chr1-123-A-G,chrX-456-C-T" --kuzu_db_path ./kuzu_db
```
Retrieves sample and zygosity information for the specified variant IDs from the Kuzu graph.

## End-to-End (E2E) Testing Status

The VCF Analysis Agent now includes a comprehensive suite of End-to-End (E2E) tests located in `tests/integration/e2e/`. These tests validate:
- Core VCF ingestion workflows via both CLI and API.
- AI-driven enrichment placeholders (further implementation pending).
- Querying capabilities for both LanceDB (vector search, metadata filtering) and Kuzu (graph context retrieval).
- Robust error handling for scenarios such as database unavailability (LanceDB, Kuzu), read-only filesystems, unreadable input files, and malformed CLI arguments.
- Security checks, including basic SQL injection prevention for LanceDB queries.

The E2E tests ensure the reliability and resilience of the agent across its key functionalities and integrations.

## Developer Documentation

Comprehensive developer documentation is available in the Sphinx docs, including:

- API reference for all core modules (autodoc)
- CLI usage and examples for all LanceDB commands
- Concurrency model (with Mermaid diagrams)
- Security framework summary
- Variant schema and LanceDB table structure
- Complete observability stack documentation

**To build the docs:**
```bash
cd docs
make html
# Open docs/build/html/index.html in your browser
```
See `docs/source/lancedb_developer_guide.rst` for the LanceDB developer guide and `docs/source/monitoring_with_prometheus.md` for observability documentation. Mermaid diagrams are supported via `sphinxcontrib-mermaid`.

## CI/CD Workflows

All CI/CD workflows are located in `kestra/flows/`. Example: `python-pip-audit.yml` runs dependency vulnerability scanning and generates a CycloneDX SBOM.

To add new workflows, place your YAML files in `kestra/flows/`.

## Development

This project is under active development. See the [Project Requirements Document](PRD%20-%20VCF%20Analysis%20Agent.md) for detailed information about the project scope and roadmap.

### Development Setup

1. Install development dependencies:
   ```bash
   pip install -r requirements-dev.txt
   ```

2. Set up pre-commit hooks:
   ```bash
   pre-commit install
   ```

3. Run tests:
   ```bash
   pytest
   ```

### Setting up for Development
```bash
# (Ensure you have followed installation steps above)

# Activate virtual environment
source .venv/bin/activate

# Install development dependencies (if any, specified in pyproject.toml dev group)
# uv pip install -e .[dev] # If you have a [project.optional-dependencies] dev group
```

### Running Tests
```bash
pytest
```
This will run all unit and integration tests. For E2E tests specifically:
```bash
pytest tests/integration/e2e/
```

### Code Formatting and Linting
(Details about Black, Ruff, etc., if used, would go here. Assuming basic setup for now.)

## Known Issues and Workarounds

### Kuzu QueryResult Lifetime
A potential segmentation fault can occur with the Kuzu Python bindings if `QueryResult` objects are not explicitly deleted (`del query_result`) before their parent `Connection` or `Database` object is closed or goes out of scope. This is due to a suspected use-after-free issue in the Kuzu C++ layer.

**Workaround:** Always ensure `del` is called on `QueryResult` objects, followed by `gc.collect()`, as soon as they are no longer needed, especially in test fixtures and long-running processes.
**Details:** See the [Kuzu Bug Report](kuzu_bug_report.md) for a more detailed explanation and a reproducible example. The issue is tracked on GitHub at [kuzudb/kuzu#5457](https://github.com/kuzudb/kuzu/issues/5457).

## License

This project is licensed under the terms of the LICENSE file included in the repository.

## Contributing

Contributions are welcome! Please read our contributing guidelines (coming soon) before submitting pull requests.

## Project Status

This project is currently in active development. See the Projects tab for current progress and planned features.

- **Core Scaffolding:** Initial agent structure, CLI, and a basic "echo" tool are functional.
- **Local LLM:** The agent is configured to use `qwen3:latest` via Ollama (ensure your Ollama server at `<ORBSTACK_SERVER_IP_PLACEHOLDER>` or your configured IP is running with this model).
- **CI/CD (Kestra):** Basic Kestra workflow for build/lint/test is defined but on HOLD pending setup resolution.

## Project Structure

The codebase follows modern Python best practices for maintainability and extensibility.