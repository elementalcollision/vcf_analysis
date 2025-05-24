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
