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
├── sample_data/                       # Place sample VCFs and indexes here (not committed)
```

- **src/vcf_agent/agent.py**: Defines the Strands agent, system prompt, and tool integration.
- **src/vcf_agent/cli.py**: Command-line interface for interacting with the agent.
- **src/vcf_agent/config.py**: Placeholder for future configuration logic.
- **src/vcf_agent/bcftools_integration.py**: Stub for integrating bcftools commands.
- **tests/**: Contains all tests, mirroring the main package structure.

This structure supports modular development, easy testing, and future extensibility.

## bcftools Python Wrapper

- The agent includes a Python wrapper for `bcftools` commands (`view`, `query`, `filter`, `norm`, `stats`, `annotate`) in `src/vcf_agent/bcftools_integration.py`.
- Each command is exposed as a function: `bcftools_view`, `bcftools_query`, `bcftools_filter`, `bcftools_norm`, `bcftools_stats`, `bcftools_annotate`.
- Functions accept a list of arguments (as you would pass on the command line) and optional input data.
- Output and errors are captured and returned for error handling.
- The wrapper is easily extensible for additional `bcftools` commands.
- All wrapper functions are covered by automated tests using pytest.

**Example usage:**
```python
from vcf_agent.bcftools_integration import bcftools_stats

# Run 'bcftools stats --help'
returncode, stdout, stderr = bcftools_stats(["--help"])
if returncode == 0 or returncode == 1:
    print(stdout or stderr)
else:
    print(f"Error: {stderr}")
```

See the module docstrings for more details and usage patterns.

## Compliance Checking

The agent supports robust, configurable compliance checking for VCF/BCF files using multiple tools:

- **bcftools** (default): Fast, standard validation for VCF/BCF format and content.
- **GATK ValidateVariants**: Strict compliance checking, including reference validation (requires reference FASTA).
- **Manufacturer-specific tools**: Easily extensible for tools from Oxford Nanopore, Illumina, or others.

### Tool Selection and Configuration

- The compliance checker backend is selected via `src/vcf_agent/config.py`:
  ```python
  CONFIG = {
      "compliance_tool": "bcftools",  # Options: "bcftools", "gatk", or manufacturer/tool name
      "gatk_reference": "/path/to/ref.fasta",  # Required for GATK
      "manufacturer_tool": None,  # e.g., "ont_vcf_validator"
  }
  ```
- You can override the tool at runtime by passing a `tool` argument to validation functions.

### Example CLI Usage

Validate a VCF file using the default tool (set in config):
```bash
python -m vcf_agent.cli "validate_vcf: sample_data/<your_sample.vcf.gz>"
```

Validate using GATK (ensure reference is set in config):
```python
from vcf_agent.validation import validate_vcf_file
is_valid, error = validate_vcf_file("sample_data/<your_sample.vcf.gz>", tool="gatk")
```

### Adding Manufacturer-Specific Tools
- Implement a wrapper in `validation.py` or a new module.
- Set `compliance_tool` in config to the tool name.
- The compliance interface will route validation to your tool.

### Reference
- Official VCF/SAM/CRAM specifications: [hts-specs](https://github.com/samtools/hts-specs)

## Batch Compliance Checking & CI/CD Integration

The repository includes a batch compliance checker script for automated validation of VCF/BCF files:

### scripts/batch_compliance_check.py
- Scans a directory (default: `sample_data/`) for VCF/BCF files
- Optionally generates synthetic edge-case files with `--generate-edgecases`
- Runs compliance checks using the selected tool (`bcftools`, `gatk`, or as configured)
- Outputs a markdown summary report to stdout or a file
- Handles errors and missing tools gracefully

#### CLI Usage
```bash
# Generate edge-case files and run compliance check with default tool
python scripts/batch_compliance_check.py --generate-edgecases

# Run with GATK and output to a markdown file
python scripts/batch_compliance_check.py -t gatk -o compliance_report.md

# Run with bcftools on a custom directory
python scripts/batch_compliance_check.py -d my_vcfs/ -t bcftools
```

#### CLI Options
- `-d, --directory`: Directory to scan for VCF/BCF files (default: `sample_data`)
- `-t, --tool`: Compliance tool to use (`bcftools`, `gatk`, or as configured)
- `-o, --output`: Output markdown file (default: stdout)
- `--generate-edgecases`: Generate synthetic edge-case VCF files in the directory

#### CI/CD Integration (Kestra)
- The workflow `.context/ci/kestra_batch_compliance.yml` automates:
  - Edge-case file generation
  - Compliance checks with both `bcftools` and `gatk`
  - Archiving of markdown reports
  - Failing the workflow if any compliance check fails

## Testing and Configuration

- All wrapper functions are tested using pytest. To run the tests:
  ```bash
  pytest
  ```
- The project uses a `pytest.ini` file for configuration. This includes:
  ```ini
  [pytest]
  asyncio_default_fixture_loop_scope = function
  ```
  This setting is required for compatibility with pytest-asyncio and to avoid deprecation warnings in async test environments.

### CLI Subprocess Mocking for Tests

- CLI tests use a robust subprocess mocking strategy to ensure thorough coverage and reliability.
- When running CLI tests, the environment variable `VCF_AGENT_CLI_MOCK_RESPONSE` is set in the subprocess. If this variable is present, the CLI prints its value and exits, bypassing the real agent logic.
- This allows tests to assert CLI output deterministically, even in subprocesses, and ensures that CLI entrypoint code is fully covered.
- Example (from test):
  ```python
  import subprocess, sys, os
  result = subprocess.run(
      [sys.executable, '-m', 'vcf_agent.cli', 'echo: CLI test!'],
      capture_output=True, text=True,
      env={**os.environ, "VCF_AGENT_CLI_MOCK_RESPONSE": "Echo: CLI test!"}
  )
  assert result.returncode == 0
  assert 'Echo: CLI test!' in result.stdout
  ```
- This approach is used for all CLI tool tests, including validation and echo tools.

### Coverage

- The project maintains high test coverage, including for CLI and configuration modules.
- Coverage is measured using `pytest-cov`:
  ```bash
  PYTHONPATH=src pytest --cov=src/vcf_agent --cov-report=term-missing --cov-report=xml
  ```
- Coverage reports are available in the terminal and as `coverage.xml` for CI integration.

## Sample Data for Testing

- Sample VCF files for testing should be placed in the `sample_data/` directory (which is gitignored).
- A sample file is required for validation tests (see `tests/test_validation.py`).
- You can obtain example VCF files from: https://bioinformaticstools.mayo.edu/research/vcf-miner-sample-vcfs/
- Index files (`.csi`/`.tbi`) may be generated for these files using `bcftools index`.

## Golden-File Regression Tests

Golden-file tests ensure that the output of wrappers (e.g., bcftools_stats) on known input files remains consistent over time.

- **Inputs:** Place sample VCF/BCF files in `sample_data/` (e.g., `sample_data/<your_sample.vcf.gz>`).
- **Golden outputs:** Place expected outputs in `golden/` (e.g., `golden/<your_sample.stats.txt>`). These files are gitignored by default.
- **How to generate/update:**
  1. Run the wrapper on your sample file and save the output:
     ```bash
     python -m vcf_agent.cli "validate_vcf: sample_data/<your_sample.vcf.gz>" > golden/<your_sample.stats.txt>
     # Or use a Python script to call the wrapper and write output
     ```
  2. If the output changes intentionally, update the golden file and review the diff.
- **How to write a golden-file test:**
  ```python
  import os
  import pytest
  from vcf_agent.bcftools_integration import bcftools_stats

  def normalize(text):
      return '\n'.join(line.strip() for line in text.strip().splitlines() if line.strip())

  SAMPLE_VCF = "sample_data/<your_sample.vcf.gz>"
  GOLDEN_STATS = "golden/<your_sample.stats.txt>"

  @pytest.mark.skipif(not os.path.exists(SAMPLE_VCF), reason="Sample VCF file not found")
  @pytest.mark.skipif(not os.path.exists(GOLDEN_STATS), reason="Golden stats file not found")
  def test_bcftools_stats_golden():
      rc, out, err = bcftools_stats([SAMPLE_VCF])
      actual = out if out else err
      with open(GOLDEN_STATS) as f:
          golden = f.read()
      assert normalize(actual) == normalize(golden), "Output does not match golden file. Update the golden file if this is intentional."
  ```
- **Interpreting failures:** If a golden-file test fails, review the diff. If the change is expected, regenerate the golden file. If not, investigate for regressions. 

## References

- [SAM/BAM/CRAM/VCF Specifications (hts-specs)](https://github.com/samtools/hts-specs) 

## Edge-case documentation

Edge-case documentation is maintained in `.context/plan/EDGECASE_VCF_DESCRIPTIONS.md` (moved from `sample_data/` to ensure it is versioned and not masked by `.gitignore`). 