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

You can build and run the VCF Analysis Agent in a containerized environment using Docker or OrbStack.

### Build the Docker image
```bash
./scripts/build.sh
```

### Run the container (with OrbStack or Docker)
```bash
./scripts/run.sh
```
- This mounts your project directory to `/app` in the container.
- The container exposes port 8000 by default (edit if needed).
- The OrbStack server is set to `<ORBSTACK_SERVER_IP>` (edit `scripts/run.sh` if your server changes).

### OrbStack Configuration
- The `orbstack.yaml` file is provided for local development and multi-service orchestration.
- Ensure your Ollama service within OrbStack (or your local Ollama server) is accessible and running `qwen3:latest` if you intend to use the default agent configuration.
- You can use OrbStack's UI or CLI to manage and launch the stack.

## Usage

Documentation and usage examples coming soon.

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

## License

This project is licensed under the terms of the LICENSE file included in the repository.

## Contributing

Contributions are welcome! Please read our contributing guidelines (coming soon) before submitting pull requests.

## Project Status

This project is currently in active development. See the Projects tab for current progress and planned features.

- **Core Scaffolding:** Initial agent structure, CLI, and a basic "echo" tool are functional.
- **Local LLM:** The agent is configured to use `qwen3:latest` via Ollama (ensure your Ollama server at `<ORBSTACK_SERVER_IP_PLACEHOLDER>` or your configured IP is running with this model).
- **CI/CD (Kestra):** Basic Kestra workflow for build/lint/test is defined but on HOLD pending setup resolution.

## CI/CD with Kestra

This project includes a [Kestra](https://kestra.io/) workflow for CI/CD automation, located at `.context/ci/kestra_vcf_agent.yml`.

### Register the Workflow

**With Kestra CLI:**
```bash
kestra workflow create .context/ci/kestra_vcf_agent.yml
```

**With Kestra Web UI:**
1. Go to `http://<KESTRA_SERVER_IP>:8080/ui` (replace `<KESTRA_SERVER_IP>` with your actual server IP)
2. Log in and navigate to "Workflows"
3. Click "Create" or "Import" and upload `.context/ci/kestra_vcf_agent.yml`

### Run the Workflow

- From the UI: Select the workflow (`vcf-agent-ci` in namespace `vcf.agent`) and click "Run".
- From the CLI:
  ```bash
  kestra workflow trigger vcf.agent.vcf-agent-ci
  ```

### Monitor Results
- Use the Kestra UI to view workflow runs, logs, and results for build, lint, and test steps.

### Notes
- The workflow builds the Docker image, runs linting (pre-commit), and tests (pytest) in a reproducible environment.
- The agent currently defaults to using `qwen3:latest` via Ollama for local development. Update `src/vcf_agent/agent.py` if you wish to use a different model.
- You can extend the workflow as your CI/CD needs grow.

**Note:** Replace `<KESTRA_SERVER_IP>` and `<ORBSTACK_SERVER_IP>` with your actual server IP addresses in your local environment. Do not commit real IP addresses to public repositories.

## Project Structure

The codebase follows modern Python best practices for maintainability and extensibility:

```
vcf_analysis/
├── src/
│   └── vcf_agent/
│       ├── __init__.py                # Package marker
│       ├── agent.py                   # Main Strands agent scaffold
│       ├── cli.py                     # CLI entrypoint for agent interaction
│       ├── config.py                  # Configuration management (placeholder)
│       └── bcftools_integration.py    # Placeholder for bcftools integration
├── tests/
│   ├── __init__.py                    # Test package marker
│   └── test_agent.py                  # Basic agent tests (pytest)
```

- **src/vcf_agent/agent.py**: Defines the Strands agent, system prompt, and tool integration.
- **src/vcf_agent/cli.py**: Command-line interface for interacting with the agent.
- **src/vcf_agent/config.py**: Placeholder for future configuration logic.
- **src/vcf_agent/bcftools_integration.py**: Stub for integrating bcftools commands.
- **tests/**: Contains all tests, mirroring the main package structure.

This structure supports modular development, easy testing, and future extensibility. 