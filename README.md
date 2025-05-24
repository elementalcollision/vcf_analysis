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
python -m vcf_agent.cli delete-variants --filter_sql "variant_id = 'rs123'"
# or, for more complex deletions:
python -m vcf_agent.cli delete-variants --filter_sql "chrom = '1' AND pos < 10000 AND clinical_significance = 'Uncertain'"
```
- The `--filter_sql` argument uses SQL-like syntax to specify which records to delete.

- The `--embedding` argument should be a comma-separated list of floats matching your embedding dimension (e.g., 1024).

See the [Implementation Plan](.context/tasks/active/TASK-005-01.md) and [LanceDB Decision Record](.context/decisions/2025-05-24-lancedb-integration.md) for more details on schema and integration best practices.

### Deleting Variants

To delete variants matching a SQL filter:

```bash
python -m src.vcf_agent.cli delete-variants --db_path ./my_lancedb --table_name my_variants_table --filter_sql "clinical_significance = 'Uncertain significance'"
```

### Creating Scalar Indexes

To create a scalar index on a column (e.g., `clinical_significance`) to potentially speed up filtered queries:

```bash
python -m src.vcf_agent.cli create-lancedb-index --db_path ./my_lancedb --table_name my_variants_table --column clinical_significance
```

To create a specific type of scalar index (e.g., BITMAP) and replace it if it already exists:

```bash
python -m src.vcf_agent.cli create-lancedb-index --db_path ./my_lancedb --table_name my_variants_table --column chrom --index_type BITMAP --replace
```

### Advanced Filtering by Metadata

To filter variants based on complex SQL-like conditions on metadata fields:

```bash
python -m src.vcf_agent.cli filter-lancedb --db_path ./my_lancedb --table_name my_variants_table --filter_sql "chrom = '1' AND pos > 10000 AND clinical_significance = 'Pathogenic'"
```

To select specific columns and limit the number of results:

```bash
python -m src.vcf_agent.cli filter-lancedb --db_path ./my_lancedb --table_name my_variants_table --filter_sql "pos < 50000" --select_columns variant_id,chrom,pos --limit 5
```

### Important Security Considerations for LanceDB

While LanceDB provides powerful local vector storage, when dealing with sensitive VCF data and its embeddings, several security aspects require careful attention, especially since the Python client for local LanceDB relies heavily on the security of the underlying system and application-level practices. The VCF Agent has a defined [Auditable Security Framework for LanceDB Integration](.context/decisions/2025-05-24-lancedb-security-framework.md) which outlines comprehensive strategies. Key areas to be mindful of during deployment and operation include:

*   **Filesystem Auditing, ACLs, and RBAC:**
    *   **Access Control (ACLs):** LanceDB data is stored as files (e.g., in the `./lancedb` directory by default). It is critical to implement strict Access Control Lists (ACLs) on this directory and its contents. Only the VCF Agent service account/user should have read/write permissions. Unauthorized access to these files could lead to data breaches or tampering.
    *   **Filesystem Auditing:** Enable OS-level audit logging (e.g., `auditd` on Linux, or similar tools on macOS/Windows) for the LanceDB data directory. This helps track all access (reads, writes, permission changes) and aids in forensic analysis if an incident occurs.
    *   **Role-Based Access Control (RBAC):** For the VCF Agent application itself, ensure that any user-facing interfaces or APIs that trigger LanceDB operations enforce appropriate RBAC, so users can only perform actions and access data they are authorized for. Direct file-level access is not a substitute for application-level RBAC if multiple users or roles interact with the agent.

*   **Masking Sensitive Data in Logs and Queries:**
    *   **Logging:** The VCF Agent aims to log LanceDB operations. While parameters like `variant_id` or counts are logged, avoid logging full embedding vectors or potentially sensitive details from `filter_sql` clauses or `updates` dictionaries directly if they might contain PII or other sensitive identifiers. The application logging has been enhanced to summarize or log keys for such parameters. However, `filter_sql` parameters are logged and marked with a TODO for potential masking in production if they contain PII.
    *   **Query Construction:** Be cautious when constructing `filter_sql` strings, especially if they incorporate user-supplied input, to prevent injection-like vulnerabilities if the query language were more expressive (though LanceDB's SQL dialect is limited, caution is always good practice).

*   **Vulnerability Scanning and SBOM:**
    *   **Dependency Scanning:** Regularly scan all project dependencies, including `lancedb` itself and its transitive dependencies (like `pyarrow`), for known vulnerabilities using tools like `pip-audit`, Snyk, or Trivy. Keep dependencies updated.
    *   **Software Bill of Materials (SBOM):** Maintain an SBOM for the VCF Agent. This helps in tracking all components and quickly identifying if a newly discovered vulnerability in a dependency affects the agent.
    *   **Host OS:** Ensure the host operating system running the VCF Agent and LanceDB is also regularly scanned and patched.

*   **Encryption:**
    *   As outlined in the security framework, data at rest should be protected using filesystem-level encryption (e.g., LUKS, BitLocker, APFS encryption). This is crucial for local LanceDB deployments.

Implementing these measures, as detailed in the full security framework, is essential for maintaining the confidentiality, integrity, and availability of genomic data managed by the VCF Agent with LanceDB.

## Developer Documentation

Comprehensive developer documentation is available in the Sphinx docs, including:

- API reference for all core modules (autodoc)
- CLI usage and examples for all LanceDB commands
- Concurrency model (with Mermaid diagrams)
- Security framework summary
- Variant schema and LanceDB table structure

**To build the docs:**
```bash
cd docs
sphinx-build -b html source build/html
# Open docs/build/html/index.html in your browser
```
See `docs/source/lancedb_developer_guide.rst` for the LanceDB developer guide. Mermaid diagrams are supported via `sphinxcontrib-mermaid`.

## CI/CD Workflows

All CI/CD workflows are located in `kestra/flows/`. Example: `python-pip-audit.yml` runs dependency vulnerability scanning and generates a CycloneDX SBOM.

To add new workflows, place your YAML files in `kestra/flows/`.

## Development

This project is under active development. See the [Project Requirements Document](PRD%20-%20VCF%20Analysis%20Agent.md) for detailed information about the project scope and roadmap.

**Full developer documentation and API guide:** See [docs/DEVELOPER_GUIDE.md](docs/DEVELOPER_GUIDE.md)

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

## Project Structure

The codebase follows modern Python best practices for maintainability and extensibility:

## Agent Integration, Prompt Contracts, API Integration, and Testing

This project follows robust, open-source-aligned best practices for agent tool registration, prompt contract usage, API integration (with security), and onboarding/testing. Below is a summary of the framework and examples:

### 1. Agent Tool Registration
- Use the `@tool` decorator to register tools with clear docstrings and type hints.
- Document each tool's name, description, input/output schema, version, and dependencies.
- For complex tools, define a `TOOL_SPEC` dictionary with input schema and metadata.
- Register tools in a central list (see `src/vcf_agent/agent.py`).
- Reference: [Strands SDK Tool Registration](https://github.com/strands-agents/docs/blob/main/docs/user-guide/concepts/tools/tools_overview.md)

**Example:**
```python
from strands import Agent, tool

@tool
def echo(text: str) -> str:
    """Echoes the input text back to the user."""
    return f"Echo: {text}"

agent = Agent(tools=[echo])
```

### 2. Prompt Contract Usage
- Store prompt contracts as versioned YAML files in `prompts/`.
- Each contract includes: `id`, `version`, `compliance`, `changelog`, `role`, `context`, `instructions`, `constraints`, `evaluation`, `test_cases`, `json_schema`.
- Document contract structure and provide usage examples in `prompts/README.md`.
- Use the loader utility in `src/vcf_agent/prompt_templates.py` to render prompts.
- Require all AI outputs to match the contract's JSON schema for determinism and testability.
- Test contracts using pytest (see `tests/prompt_contracts/`).

**Example YAML:**
```yaml
id: vcf_analysis_summary_v1
version: 1.0.0
role: "Genomics Analyst"
context: "Summarize VCF for quality and compliance"
instructions: "Return only valid JSON matching the schema."
json_schema: {...}
test_cases:
  - input: sample_data/HG00098.vcf.gz
    expected_output: tests/golden/vcf_summary_HG00098.json
```

### 3. API Integration & Security
- Use secure credential management (env vars, JSON files; never hardcode keys).
- Document API endpoints, authentication methods, and required headers.
- Enforce TLS for all network communication.
- Rate-limit API calls and log all access for auditability.
- Provide code examples for API usage (see `src/vcf_agent/api_clients.py`).

**Example:**
```python
from vcf_agent.api_clients import OpenAIClient
client = OpenAIClient(api_key="sk-...")
response = client.chat_completion(messages=[...])
```

### 4. Onboarding & Testing Strategy
- Provide a clear onboarding checklist (see `docs/DEVELOPER_GUIDE.md`):
  - Environment setup, dependency installation, running tests, reviewing docstrings.
- Document test types: unit, integration, contract, edge case, and UAT.
- Use golden files and property-based tests for regression and compliance.
- Automate test runs in CI/CD (see Kestra workflow in `.context/ci/`).
- Document how to add new tools, prompt contracts, and tests.

**Example Checklist:**
```markdown
- [ ] Clone repo and set up Python environment
- [ ] Install dependencies (`uv pip install -r requirements.txt`)
- [ ] Run all tests (`pytest`)
- [ ] Review prompt contracts and schemas
- [ ] Add new tools using @tool decorator and update docs
```

### 5. Open Source & Community Standards
- Reference open-source agent frameworks (Strands, SuperAGI, Superagent, Relari Agent Contracts) for structure and compliance.
- Use semantic versioning and changelogs for all contracts and tools.
- Encourage contributions with clear guidelines and code examples.

For more details, see:
- `prompts/README.md` for prompt contract structure and examples
- `src/vcf_agent/agent.py` and `src/vcf_agent/prompt_templates.py` for tool and prompt integration
- `docs/DEVELOPER_GUIDE.md` for onboarding and extension guidelines
- [Strands Agents SDK Documentation](https://github.com/strands-agents/docs)
