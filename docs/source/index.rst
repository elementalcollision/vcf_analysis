.. VCF Analysis Agent documentation master file, created by
   sphinx-quickstart on Wed May 21 23:24:28 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

VCF Analysis Agent documentation
================================

Add your content using ``reStructuredText`` syntax. See the
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
documentation for details.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

Agent Integration, Prompt Contracts, API Integration, and Testing
===============================================================

The VCF Analysis Agent follows robust best practices for agent tool registration, prompt contract usage, API integration (with security), and onboarding/testing:

- **Tool Registration:** Use decorators and schemas to register tools with clear input/output definitions. See `src/vcf_agent/agent.py`.
- **Prompt Contracts:** Store versioned YAML contracts in `prompts/`, each with required fields, schemas, and test cases. See `prompts/README.md`.
- **API Integration & Security:** Use secure credential management, document endpoints and authentication, and enforce TLS. See `src/vcf_agent/api_clients.py`.
- **Onboarding & Testing:** Provide a clear checklist, document test types, use golden files, and automate CI/CD. See `docs/DEVELOPER_GUIDE.md`.
- **Open Source Standards:** Reference frameworks like Strands, SuperAGI, and Relari Agent Contracts for structure and compliance.

For full details and examples, see the 'Agent Integration, Prompt Contracts, API Integration, and Testing' section in the project README.

