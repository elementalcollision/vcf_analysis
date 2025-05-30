validate_cli_documentation
==========================

.. py:module:: validate_cli_documentation

.. autoapi-nested-parse::

   CLI Documentation Validation Script

   Validates that the CLI module docstring accurately documents all implemented commands.
   This prevents documentation drift and ensures developers can trust the module documentation.

   Usage:
       python scripts/validate_cli_documentation.py [--check-completeness] [--verbose]

   Exit codes:
       0: All documentation is complete and accurate
       1: Documentation discrepancies found
       2: Script execution error



Classes
-------

.. autoapisummary::

   validate_cli_documentation.CLIDocumentationValidator
   validate_cli_documentation.Command
   validate_cli_documentation.ValidationResult


Functions
---------

.. autoapisummary::

   validate_cli_documentation.main
   validate_cli_documentation.print_validation_report


Module Contents
---------------

.. py:class:: CLIDocumentationValidator(cli_module_path = 'src/vcf_agent/cli.py')

   Validates CLI module documentation against implementation.


   .. py:method:: extract_commands_from_docstring()

      Extract command names from the module docstring.



   .. py:method:: extract_commands_from_implementation()

      Extract top-level command names from argparse implementation.



   .. py:method:: extract_samspec_subcommands()

      Extract samspec subcommands for special handling.



   .. py:method:: validate()

      Perform complete validation of CLI documentation.



   .. py:method:: validate_command_examples()

      Validate that command examples in docstring use actual commands.



   .. py:method:: validate_samspec_subcommands_documented()

      Validate that samspec subcommands are mentioned in docstring.



   .. py:attribute:: cli_module_path


.. py:class:: Command

   Represents a CLI command with its metadata.


   .. py:attribute:: category
      :type:  Optional[str]
      :value: None



   .. py:attribute:: description
      :type:  str


   .. py:attribute:: is_subcommand
      :type:  bool
      :value: False



   .. py:attribute:: name
      :type:  str


   .. py:attribute:: parent_command
      :type:  Optional[str]
      :value: None



.. py:class:: ValidationResult

   Results of CLI documentation validation.


   .. py:attribute:: documented_but_not_implemented
      :type:  Set[str]


   .. py:attribute:: documented_commands
      :type:  Set[str]


   .. py:attribute:: errors
      :type:  List[str]


   .. py:attribute:: implemented_commands
      :type:  Set[str]


   .. py:attribute:: is_valid
      :type:  bool


   .. py:attribute:: missing_from_docs
      :type:  Set[str]


   .. py:attribute:: warnings
      :type:  List[str]


.. py:function:: main()

   Main function for CLI documentation validation.


.. py:function:: print_validation_report(result, verbose = False)

   Print a formatted validation report.


