# VCF Analysis Agent – Developer Guide

## 1. Project Overview
- **Purpose:** The VCF Analysis Agent automates VCF/BCF validation, annotation, and analysis, integrating bcftools, GATK, and AI.
- **Architecture:** Modular Python package with CLI, API, and extensible wrappers for external tools.

## 2. Installation & Configuration
- **Environment:** Python 3.11+, virtualenv, `uv` for dependencies.
- **Install:**
  ```bash
  python -m venv .venv
  source .venv/bin/activate
  pip install uv
  uv pip install -r requirements.txt
  ```
- **Config:**
  Edit `src/vcf_agent/config.py` for compliance tool selection and references.

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
- **Sphinx Setup:**
  - Install: `pip install sphinx sphinx_rtd_theme`
  - Quickstart: `sphinx-quickstart docs`
  - Enable extensions in `docs/source/conf.py`:
    ```python
    extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.napoleon',
        'sphinx.ext.autosummary'
    ]
    ```
  - Generate API docs:
    ```bash
    sphinx-apidoc -f -o docs/source src/vcf_agent
    sphinx-build -M html docs/source docs/build/
    ```

## 4. CLI Usage
- **Entrypoint:**
  ```bash
  python -m vcf_agent.cli "validate_vcf: sample_data/<your_sample.vcf.gz>"
  ```
- **Example:**
  ```bash
  python -m vcf_agent.cli "echo: Hello, world!"
  ```
- **Error Handling:**
  All CLI tools return clear error messages and exit codes.

## 5. Testing Strategy
- **Framework:** pytest
- **Coverage:** Aim for ≥90% on core modules.
- **Test Types:**
  - Unit: Each wrapper and validation function.
  - Integration: End-to-end CLI and compliance checks.
  - Edge Cases: Use `.context/plan/EDGECASE_VCF_DESCRIPTIONS.md` as canonical source.
- **Example Test:**
  ```python
  def test_bcftools_stats_golden():
      rc, out, err = bcftools_stats(["sample_data/HG00098.vcf.gz"])
      assert rc == 0
      assert "lines" in out
  ```

## 6. Extension Guidelines
- **Adding a New bcftools Command:**
  1. Add a function in `bcftools_integration.py`:
     ```python
     def bcftools_newcmd(args: List[str], input_data: Optional[bytes] = None) -> Tuple[int, str, str]:
         return run_bcftools_command(["newcmd"] + args, input_data)
     ```
  2. Add tests in `tests/`.
  3. Document usage in the API docs and README.
- **Plugin Architecture:**
  Follow the pattern in `validation.py` for new tool integrations.

## 7. Onboarding Notes
- **Checklist:**
  - Set up Python environment and install dependencies.
  - Run tests: `pytest`
  - Review docstrings and Sphinx docs for API usage.
- **Debugging:**
  - If Sphinx import errors occur, add project root to `sys.path` in `conf.py`.
  - Use `sphinx-autobuild` for live doc preview.

## 8. Maintenance & Contribution
- **Docs:**
  - Keep `/docs` and code in sync.
  - Use Sphinx for HTML and API docs.
- **CI/CD:**
  - Add Sphinx build to CI pipeline.
- **Changelog:**
  - Update for every release.
- **Deprecation:**
  - Mark deprecated APIs in docstrings and docs.

## 9. Appendices
- **Glossary:**
  - VCF, BCF, bcftools, GATK, etc.
- **Release Checklist:**
  - All tests pass, docs build, changelog updated.

## 11. Agent Integration, Prompt Contracts, API Integration, and Testing

This project follows robust best practices for agent tool registration, prompt contract usage, API integration (with security), and onboarding/testing:

- **Tool Registration:** Use decorators and schemas to register tools with clear input/output definitions. See `src/vcf_agent/agent.py`.
- **Prompt Contracts:** Store versioned YAML contracts in `prompts/`, each with required fields, schemas, and test cases. See `prompts/README.md`.
- **API Integration & Security:** Use secure credential management, document endpoints and authentication, and enforce TLS. See `src/vcf_agent/api_clients.py`.
- **Onboarding & Testing:** Provide a clear checklist, document test types, use golden files, and automate CI/CD. See earlier sections and `docs/DEVELOPER_GUIDE.md`.
- **Open Source Standards:** Reference frameworks like Strands, SuperAGI, and Relari Agent Contracts for structure and compliance.

For full details and examples, see the 'Agent Integration, Prompt Contracts, API Integration, and Testing' section in the project README.

---

**Tips for Clarity & Maintainability:**
- Use consistent docstring style (Google/NumPy).
- Keep examples up to date and testable (doctest or pytest).
- Use Sphinx cross-references for classes/functions.
- Document extension points and plugin patterns.
- Onboard new contributors with a clear checklist and troubleshooting section. 