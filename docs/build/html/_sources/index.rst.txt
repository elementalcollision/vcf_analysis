.. VCF Analysis Agent documentation master file, created by
   sphinx-quickstart on Wed May 21 23:24:28 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

VCF Analysis Agent documentation
================================

Welcome to the VCF Analysis Agent documentation. This project provides an AI-powered tool for intelligent analysis, validation, and processing of Variant Call Format (VCF) files.

Core functionalities include:
- VCF file processing and validation.
- AI-driven variant analysis and interpretation (ongoing development).
- Integration with LanceDB for vector-based similarity searches of variants.
- Integration with Kuzu graph database for contextual variant and sample relationship queries.
- Robust error handling and a comprehensive End-to-End testing suite.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self # Includes content below
   vcf_agent
   lancedb_developer_guide
   kuzu_developer_guide
   monitoring_with_prometheus
   security
   audit
   modules

Agent Integration, Prompt Contracts, and API Usage
===================================================

The VCF Analysis Agent follows robust best practices for agent tool registration, prompt contract usage, and API integration (with security):

- **Tool Registration:** Uses decorators and schemas to register tools with clear input/output definitions. See `src/vcf_agent/agent.py`.
- **Prompt Contracts:** Stores versioned YAML contracts in `prompts/`, each with required fields, schemas, and test cases. See `prompts/README.md`.
- **API Integration & Security:** Employs secure credential management for LLM providers. See the main project `README.md` and `docs/source/security.rst` for details.
- **Onboarding & Testing:** Provides a clear checklist, documents test types (unit, integration, E2E), uses golden files, and automates CI/CD. See `docs/DEVELOPER_GUIDE.md`.

Containerization: Docker, Multi-Arch, and OrbStack
===================================================

The VCF Analysis Agent is containerized for production, local development, and CI/CD. Key features include multi-stage builds, non-root user execution, and multi-architecture support (`linux/amd64`, `linux/arm64`). For full details on building, running, and best practices, please refer to the 'Containerization (Docker & OrbStack)' section in the main project `README.md`.

Security and Auditing
=====================

Detailed information on security best practices, credential management, data handling, and auditing can be found in:

- :doc:`security`
- :doc:`audit`

Further specific guidance related to database security is available in the respective developer guides for LanceDB and Kuzu.

LanceDB Developer Guide
=======================

.. toctree::
   :maxdepth: 2

   lancedb_developer_guide.rst

