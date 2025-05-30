# VCF Analysis Agent – Developer Guide

## 1. Project Overview
- **Purpose:** The VCF Analysis Agent automates VCF/BCF validation, annotation, and analysis, integrating bcftools, GATK, and AI.
- **Architecture:** Modular Python package with CLI, API, and extensible wrappers for external tools.

## 2. Installation & Configuration
- **Environment:** Python 3.11+, virtualenv, `uv` for dependencies.
- **Install:** Refer to the main `README.md` for detailed installation instructions using `uv`.
- **Config:** Core agent configuration is managed via `src/vcf_agent/config.py`. LLM provider selection and credentials can be set via environment variables, a JSON credentials file, or `SessionConfig` as described in the main `README.md`.

## 3. API Documentation (with Sphinx)
- **Docstring Style:** Use Google or NumPy style for all public functions/classes.
- **Example:**
  ```python
  def validate_vcf_file(filepath: str, tool: Optional[str] = None) -> Tuple[bool, Optional[str]]:
      """
      High-level validation: checks existence, format, index, and compliance (configurable backend).

      Args:
          filepath: Path to the VCF/BCF file.
          tool: Compliance tool to use ('bcftools', 'gatk', etc.).

      Returns:
          Tuple[bool, Optional[str]]: (is_valid, error_message_if_any)

      Example:
          >>> validate_vcf_file('sample_data/HG00098.vcf.gz')
          (True, None)
      """
  ```
- **Sphinx Setup & Generation:**
  - Ensure Sphinx and necessary extensions (`sphinx.ext.autodoc`, `sphinx.ext.napoleon`, `sphinx.ext.autosummary`, `sphinx_rtd_theme`, `sphinx_copybutton`, `sphinx_design`, `sphinx_external_toc`, `sphinx_tabs`, `sphinx_togglebutton`, `sphinxcontrib.mermaid`) are installed (see `pyproject.toml` or `requirements.txt`).
  - To build documentation:
    ```bash
    sphinx-build -M html docs/source docs/build/
    ```
  - For auto-regeneration of API docs from docstrings before building:
    ```bash
    sphinx-apidoc -f -o docs/source src/vcf_agent
    ```

## 3.5. Docker Development Environment

### Quick Start with Docker

The VCF Analysis Agent provides a complete Docker development environment with hot reloading, debugging tools, and full observability stack integration.

#### Development Environment Setup

```bash
# Start development environment
docker-compose --profile development up -d

# Access development container with mounted source code
docker-compose exec vcf-agent-dev bash

# Install additional development packages
docker-compose exec vcf-agent-dev pip install ipython jupyter

# Start Jupyter notebook for interactive development
docker-compose exec vcf-agent-dev jupyter notebook --ip=0.0.0.0 --port=8888
```

#### Development Workflow

```bash
# 1. Start development stack
docker-compose --profile development up -d

# 2. Make code changes (automatically reflected in container)
# Edit files in your local IDE - changes are mounted into container

# 3. Run tests in container
docker-compose exec vcf-agent-dev pytest tests/unit/ -v

# 4. Test CLI commands
docker-compose exec vcf-agent-dev python -m vcf_agent.cli --help

# 5. Debug with interactive Python
docker-compose exec vcf-agent-dev python -c "
import vcf_agent
from vcf_agent.cli import main
# Interactive debugging
"
```

### Building and Testing Images

#### Local Development Builds

```bash
# Build development image
./scripts/docker-build.sh --target development

# Build production image
./scripts/docker-build.sh --target runtime

# Build with custom version
./scripts/docker-build.sh --version dev-$(git rev-parse --short HEAD)
```

#### Multi-Architecture Development

```bash
# Build for multiple architectures
./scripts/docker-build.sh --platform linux/amd64,linux/arm64

# Test specific architecture
./scripts/docker-build.sh --platform linux/amd64 --target development
docker run --rm -it vcf-analysis-agent:dev-amd64 bash
```

#### Security and Quality Checks

```bash
# Build with security scanning
./scripts/docker-build.sh --scan

# Manual vulnerability scanning
trivy image vcf-analysis-agent:latest

# Check image layers and size
docker history vcf-analysis-agent:latest
docker images vcf-analysis-agent --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}"
```

### Development Best Practices

#### Container Development Guidelines

1. **Source Code Mounting**: Use volume mounts for hot reloading during development
2. **Database Persistence**: Use named volumes for database data to persist across container restarts
3. **Environment Variables**: Use `.env` files for local development configuration
4. **Resource Limits**: Set appropriate CPU and memory limits for development containers

#### Debugging in Containers

```bash
# Debug container startup issues
docker-compose logs vcf-agent-dev

# Check container health
docker-compose exec vcf-agent-dev python -c "import vcf_agent; print('OK')"

# Monitor resource usage
docker stats vcf_analysis_agent_dev

# Access container filesystem
docker-compose exec vcf-agent-dev ls -la /app/

# Check environment variables
docker-compose exec vcf-agent-dev env | grep VCF
```

#### Testing in Containers

```bash
# Run full test suite in container
docker-compose exec vcf-agent-dev pytest

# Run tests with coverage
docker-compose exec vcf-agent-dev pytest --cov=src/vcf_agent --cov-report=term-missing

# Run specific test categories
docker-compose exec vcf-agent-dev pytest tests/unit/ -v
docker-compose exec vcf-agent-dev pytest tests/integration/ -v

# Test with different Python versions (if multiple images available)
docker run --rm -v $(pwd):/app/src python:3.11-slim bash -c "cd /app/src && pip install -r requirements.txt && pytest"
```

### Production Deployment Preparation

#### Image Optimization

```bash
# Build optimized production image
./scripts/docker-build.sh --target runtime --no-cache

# Verify image size and layers
docker images vcf-analysis-agent:latest
docker history vcf-analysis-agent:latest --no-trunc

# Test production image functionality
docker run --rm vcf-analysis-agent:latest python -m vcf_agent.cli --help
```

#### Configuration Management

```bash
# Create production configuration
cp config/docker/production.env .env.production

# Test with production configuration
docker run --rm --env-file .env.production vcf-analysis-agent:latest

# Validate configuration
docker run --rm --env-file .env.production vcf-analysis-agent:latest \
  python -c "from vcf_agent.config import get_config; print(get_config())"
```

### Observability in Development

#### Monitoring Stack

```bash
# Start full observability stack
docker-compose up -d

# Access monitoring dashboards
open http://localhost:3000  # Grafana (admin/admin)
open http://localhost:9090  # Prometheus
open http://localhost:16686 # Jaeger

# Check service health
curl -f http://localhost:9090/-/healthy
curl -f http://localhost:3000/api/health
```

#### Development Metrics

```bash
# View application metrics
curl http://localhost:8000/metrics

# Test tracing integration
docker-compose exec vcf-agent python -m vcf_agent.cli ask "What are the basic stats for sample_data/minimal.vcf.gz?"

# View traces in Jaeger UI
open http://localhost:16686
```

### Troubleshooting

#### Common Development Issues

1. **Port Conflicts**: Ensure ports 3000, 8000, 9090, 16686 are available
2. **Volume Permissions**: Fix with `sudo chown -R $USER:$USER ./data ./lancedb ./kuzu_db`
3. **Memory Issues**: Increase Docker memory allocation in Docker Desktop settings
4. **Build Failures**: Clear Docker cache with `docker system prune -a`

#### Debug Commands

```bash
# Check Docker daemon
docker version
docker info

# Verify Docker Compose configuration
docker-compose config

# Check container logs
docker-compose logs --tail=50 vcf-agent-dev

# Inspect container configuration
docker inspect vcf_analysis_agent_dev

# Check network connectivity
docker-compose exec vcf-agent-dev ping prometheus
docker-compose exec vcf-agent-dev curl -f http://prometheus:9090/-/healthy
```

## 4. CLI Usage
- Refer to the main `README.md` for a comprehensive list of CLI commands, their arguments, and usage examples for both agent interaction and database management (LanceDB, Kuzu).
- **Entrypoint Example:**
  ```bash
  python -m vcf_agent.cli "echo: Hello, world!"
  ```
- **Error Handling:** CLI tools provide clear error messages and non-zero exit codes on failure.

## 5. Testing Strategy
- **Framework:** pytest
- **Coverage:** Aim for >=90% on core modules. View coverage reports in `htmlcov/` after running tests with `--cov` flags (configured in `pytest.ini`).
- **Test Types:**
  - Unit: Individual functions and classes (e.g., in `tests/unit/`).
  - Integration: Interactions between components (e.g., in `tests/integration/`).
  - End-to-End (E2E): Full workflows via CLI and API, including database interactions and error scenarios (in `tests/integration/e2e/`).
- **Running Tests:**
  ```bash
  pytest  # Runs all tests
  pytest tests/unit/  # Run only unit tests
  pytest tests/integration/e2e/ # Run only E2E tests
  ```
- **Edge Cases:** Use `.context/plan/EDGECASE_VCF_DESCRIPTIONS.md` as a source for defining challenging VCF scenarios to test.

## 6. Extension Guidelines
- **Adding a New CLI Command:**
  1. Define the command, its arguments, and help text in `src/vcf_agent/cli.py` using `argparse`.
  2. Implement the command's logic, typically by calling functions from other modules.
  3. Add corresponding E2E tests in `tests/integration/e2e/`.
  4. Document the new command in the main `README.md` and relevant Sphinx guides.
- **Adding a New Agent Tool:**
  1. Define the tool's functionality as a Python function.
  2. Decorate it with `@tool` and provide a Pydantic model for its arguments, as shown in `src/vcf_agent/agent.py`.
  3. Ensure the tool is registered with the agent.
  4. Add unit tests for the tool's logic.
  5. Optionally, add integration or E2E tests if it interacts with external systems.
- **Plugin Architecture:** Follow existing patterns for integrating new functionalities or external tools (e.g., how `bcftools_integration.py` or `lancedb_integration.py` are structured and used).

## 7. Onboarding Notes
- **Checklist:**
  - Set up Python environment and install dependencies as per `README.md`.
  - Run all tests: `pytest`
  - Build Sphinx docs: `sphinx-build -M html docs/source docs/build/` and review locally.
  - Familiarize yourself with the main `README.md` and the Sphinx documentation, particularly the developer guides for LanceDB and Kuzu if working with those components.
- **Debugging:**
  - If Sphinx import errors occur, ensure `pythonpath = src` is set in `pytest.ini` and that your `PYTHONPATH` environment variable (if set) doesn't conflict. The `conf.py` for Sphinx should correctly point to the `src` directory.
  - Use `sphinx-autobuild docs/source docs/build/html` for live doc preview during editing.

## 8. Maintenance & Contribution
- **Docs:** Keep all documentation (`README.md`, Sphinx docs, `DEVELOPER_GUIDE.md`) synchronized with code changes.
- **CI/CD:** Ensure CI pipeline (e.g., GitHub Actions) runs tests, linters, and builds documentation.
- **Changelog:** Update a changelog file for significant changes or releases.
- **Deprecation:** Clearly mark deprecated APIs in docstrings and documentation, providing alternatives and a removal timeline.
- **Branching Strategy (Example):** `main` for releases, `develop` for ongoing work, feature branches off `develop`.
- **Pull Requests:** Ensure PRs are reviewed, pass CI checks, and include necessary documentation updates before merging.

## 9. Appendices
- **Glossary:** VCF, BCF, bcftools, GATK, LanceDB, Kuzu, Embedding, Cypher, SQLi, PII, SBOM etc.
- **Release Checklist:** All tests pass, docs build and are up-to-date, changelog updated, version bumped.

## 10. Key Documentation Links

- Main Project Information: `README.md`
- Sphinx Documentation (Conceptual Guides, API Reference): `docs/build/html/index.html` (after building)
  - LanceDB Integration: :doc:`docs/source/lancedb_developer_guide`
  - Kuzu Integration: :doc:`docs/source/kuzu_developer_guide`
  - Security Overview: :doc:`docs/source/security`
  - Auditing Overview: :doc:`docs/source/audit`
- Known Issues:
  - Kuzu `QueryResult` segfault: `kuzu_bug_report.md` and `docs/source/kuzu_developer_guide.rst`

---

**Tips for Clarity & Maintainability:**
- Use consistent docstring style (Google/NumPy).
- Keep examples up to date and testable (doctest or pytest).
- Use Sphinx cross-references for classes/functions.
- Document extension points and plugin patterns.
- Onboard new contributors with a clear checklist and troubleshooting section. 