vcf\_agent package
==================

Submodules
----------

vcf\_agent.agent module
-----------------------

.. automodule:: vcf_agent.agent
   :members:
   :undoc-members:
   :show-inheritance:

vcf\_agent.bcftools\_integration module
---------------------------------------

.. automodule:: vcf_agent.bcftools_integration
   :members:
   :undoc-members:
   :show-inheritance:

vcf\_agent.cli module
---------------------

.. automodule:: vcf_agent.cli
   :members:
   :undoc-members:
   :show-inheritance:

vcf\_agent.config module
------------------------

.. automodule:: vcf_agent.config
   :members:
   :undoc-members:
   :show-inheritance:

vcf\_agent.gatk\_integration module
-----------------------------------

.. automodule:: vcf_agent.gatk_integration
   :members:
   :undoc-members:
   :show-inheritance:

vcf\_agent.validation module
----------------------------

.. automodule:: vcf_agent.validation
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: vcf_agent
   :members:
   :undoc-members:
   :show-inheritance:

Agent Integration, Prompt Contracts, API Integration, and Testing
---------------------------------------------------------------

The VCF Analysis Agent implements:

- **Tool Registration:** Decorators and schemas for clear tool definitions (see `src/vcf_agent/agent.py`).
- **Prompt Contracts:** Versioned YAML contracts with schemas and test cases (`prompts/`).
- **API Integration & Security:** Secure credential management, endpoint documentation, and TLS (`src/vcf_agent/api_clients.py`).
- **Onboarding & Testing:** Checklist, test types, golden files, and CI/CD (`docs/DEVELOPER_GUIDE.md`).
- **Open Source Standards:** Structure and compliance inspired by Strands, SuperAGI, and Relari Agent Contracts.

See the project README for full details and examples.
