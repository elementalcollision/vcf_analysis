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