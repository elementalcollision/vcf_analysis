cli_documentation_style_guide
=============================

.. py:module:: cli_documentation_style_guide

.. autoapi-nested-parse::

   VCF Analysis Agent CLI Documentation Style Guide

   This module provides comprehensive documentation standards, templates, and examples
   for documenting CLI command handlers in the VCF Analysis Agent project.

   Based on:
   - PEP 257 (Docstring Conventions)
   - Google Python Style Guide
   - Click CLI Documentation Best Practices
   - argparse Integration Requirements

   Version: 1.0
   Created: 2025-01-05
   Author: VCF Analysis Agent Team



Classes
-------

.. autoapisummary::

   cli_documentation_style_guide.CLIDocstringValidator
   cli_documentation_style_guide.CLIDocumentationStandards
   cli_documentation_style_guide.DocstringValidationResult


Functions
---------

.. autoapisummary::

   cli_documentation_style_guide.cli_command_template
   cli_documentation_style_guide.compliance_command_template
   cli_documentation_style_guide.example_ingest_vcf_docstring
   cli_documentation_style_guide.example_init_lancedb_docstring
   cli_documentation_style_guide.extract_cli_commands
   cli_documentation_style_guide.generate_command_docstring
   cli_documentation_style_guide.get_style_guide_summary
   cli_documentation_style_guide.lancedb_command_template
   cli_documentation_style_guide.validate_module_documentation
   cli_documentation_style_guide.vcf_processing_template


Module Contents
---------------

.. py:class:: CLIDocstringValidator(standards = None)

   Validates CLI command handler docstrings against style guide standards.

   Provides comprehensive validation including section presence, format
   compliance, content quality, and integration requirements.

   Initialize validator with documentation standards.

   :param standards: Documentation standards to validate against.
                     Default: CLIDocumentationStandards()


   .. py:method:: validate_docstring(func)

      Validate a CLI command handler's docstring.

      :param func: Function to validate docstring for.

      :returns: Detailed validation results.
      :rtype: DocstringValidationResult

      .. rubric:: Examples

      >>> validator = CLIDocstringValidator()
      >>> result = validator.validate_docstring(some_cli_command)
      >>> if not result.is_valid:
      ...     for issue in result.issues:
      ...         print(f"Issue: {issue}")



   .. py:attribute:: standards


.. py:class:: CLIDocumentationStandards

   Core standards for CLI command handler documentation.

   These standards ensure consistency, completeness, and integration
   with both developer tools (help(), IDEs) and end-user tools (--help).


   .. py:attribute:: MAX_LINE_LENGTH
      :value: 72



   .. py:attribute:: MINIMUM_COVERAGE
      :value: 95.0



   .. py:attribute:: OPTIONAL_SECTIONS
      :value: ['Notes', 'See Also', 'References', 'Environment Variables', 'Files']



   .. py:attribute:: REQUIRED_SECTIONS
      :value: ['Summary', 'Description', 'Arguments', 'Options', 'Returns', 'Raises', 'Examples', 'Exit Codes']



.. py:class:: DocstringValidationResult

   Result of docstring validation check.


   .. py:attribute:: coverage_score
      :type:  float


   .. py:attribute:: is_valid
      :type:  bool


   .. py:attribute:: issues
      :type:  List[str]


   .. py:attribute:: missing_sections
      :type:  List[str]


   .. py:attribute:: suggestions
      :type:  List[str]


.. py:function:: cli_command_template()

   Standard template for CLI command handler docstrings.


.. py:function:: compliance_command_template()

   Template for compliance validation command handlers.


.. py:function:: example_ingest_vcf_docstring()

   Example of properly documented ingest-vcf command handler.


.. py:function:: example_init_lancedb_docstring()

   Example of properly documented init-lancedb command handler.


.. py:function:: extract_cli_commands(module_path)

   Extract CLI command handler function names from a module.

   :param module_path: Path to Python module to analyze.

   :returns: Function names that appear to be CLI command handlers.
   :rtype: List[str]

   .. rubric:: Examples

   Extract commands from the CLI module:

   >>> commands = extract_cli_commands('src/vcf_agent/cli.py')
   >>> isinstance(commands, list)
   True


.. py:function:: generate_command_docstring(command_type, **kwargs)

   Generate a docstring template for a specific command type.

   :param command_type: Type of command ('lancedb', 'vcf_processing', 'compliance', 'general').
   :param \*\*kwargs: Template substitution values.

   :returns: Generated docstring template.
   :rtype: str

   .. rubric:: Examples

   Generate a LanceDB command docstring:

   >>> docstring = generate_command_docstring('lancedb')
   >>> 'lancedb' in docstring.lower()
   True


.. py:function:: get_style_guide_summary()

   Get a summary of the CLI documentation style guide.


.. py:function:: lancedb_command_template()

   Template for LanceDB operation command handlers.


.. py:function:: validate_module_documentation(module_path)

   Validate documentation for all CLI commands in a module.

   :param module_path: Path to Python module to validate.

   :returns: Validation results by function name.
   :rtype: Dict[str, DocstringValidationResult]

   .. rubric:: Examples

   Validate all CLI commands in a module:

   >>> results = validate_module_documentation('src/vcf_agent/cli.py')
   >>> isinstance(results, dict)
   True


.. py:function:: vcf_processing_template()

   Template for VCF file processing command handlers.


