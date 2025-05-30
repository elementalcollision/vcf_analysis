cli_enhanced_validation
=======================

.. py:module:: cli_enhanced_validation

.. autoapi-nested-parse::

   Enhanced CLI Documentation Validation Engine

   Comprehensive validation engine that integrates existing CLI documentation tools
   with advanced AST parsing, multi-format docstring support, and CI/CD integration.

   This engine serves as the unified entry point for all CLI documentation validation,
   providing:
   - Integration with existing validation tools
   - Advanced AST-based analysis
   - Multi-format docstring parsing
   - Performance optimization with caching
   - CI/CD integration with proper exit codes and reporting

   Usage:
       python scripts/cli_enhanced_validation.py [--mode=comprehensive] [--format=console] [--verbose]

   Exit codes:
       0: All validation passed
       1: Validation failures found
       2: Script execution error



Attributes
----------

.. autoapisummary::

   cli_enhanced_validation.DOCSTRING_PARSER_AVAILABLE
   cli_enhanced_validation.PYDANTIC_SETTINGS_AVAILABLE
   cli_enhanced_validation.PYDANTIC_SETTINGS_AVAILABLE
   cli_enhanced_validation.TOML_SUPPORT
   cli_enhanced_validation.TOML_SUPPORT


Classes
-------

.. autoapisummary::

   cli_enhanced_validation.ASTAnalyzer
   cli_enhanced_validation.CacheManager
   cli_enhanced_validation.CodeDefinition
   cli_enhanced_validation.CoverageReport
   cli_enhanced_validation.EnhancedCLIValidator
   cli_enhanced_validation.ExceptionInfo
   cli_enhanced_validation.FileValidationResult
   cli_enhanced_validation.FunctionSignature
   cli_enhanced_validation.GitIntegration
   cli_enhanced_validation.MultiFormatDocstringParser
   cli_enhanced_validation.OutputFormat
   cli_enhanced_validation.Parameter
   cli_enhanced_validation.ParsedDocstring
   cli_enhanced_validation.PerformanceOptimizer
   cli_enhanced_validation.ReturnInfo
   cli_enhanced_validation.StructureValidation
   cli_enhanced_validation.ValidationConfig
   cli_enhanced_validation.ValidationIssue
   cli_enhanced_validation.ValidationLocation
   cli_enhanced_validation.ValidationMode
   cli_enhanced_validation.ValidationReport


Functions
---------

.. autoapisummary::

   cli_enhanced_validation.main
   cli_enhanced_validation.setup_logging


Module Contents
---------------

.. py:class:: ASTAnalyzer

   Advanced AST analyzer for comprehensive code analysis.


   .. py:method:: analyze_docstring_coverage(definitions)

      Calculate comprehensive docstring coverage metrics.



   .. py:method:: extract_all_definitions(file_path)

      Extract all functions, methods, classes from a Python file.



   .. py:method:: extract_function_signature(node)

      Extract complete function signature including parameter and return types.

      :param node: AST FunctionDef node

      :returns: Complete signature information
      :rtype: FunctionSignature



   .. py:attribute:: logger


.. py:class:: CacheManager(cache_dir = '.validation-cache', ttl_hours = 24, max_size_mb = 100)

   File-based caching system for validation results with content-hash invalidation.

   Provides efficient caching of AST parsing results and validation outcomes
   with automatic cache invalidation when file content changes.


   .. py:method:: get(file_path, validation_type = 'ast')

      Retrieve cached validation result for a file.

      :param file_path: Path to the file being validated
      :param validation_type: Type of validation (ast, module, style)

      :returns: Cached result if valid, None if cache miss



   .. py:method:: get_stats()

      Get cache performance statistics.



   .. py:method:: invalidate(file_path = None)

      Invalidate cache entries.

      :param file_path: Specific file to invalidate, or None to clear all

      :returns: Number of cache entries removed



   .. py:method:: set(file_path, result, validation_type = 'ast')

      Store validation result in cache.

      :param file_path: Path to the file being validated
      :param result: Validation result to cache
      :param validation_type: Type of validation (ast, module, style)

      :returns: True if successfully cached, False otherwise



   .. py:attribute:: cache_dir


   .. py:attribute:: hits
      :value: 0



   .. py:attribute:: lock


   .. py:attribute:: logger


   .. py:attribute:: max_size_bytes
      :value: 104857600



   .. py:attribute:: misses
      :value: 0



   .. py:attribute:: ttl_seconds
      :value: 86400



.. py:class:: CodeDefinition

   Represents a code definition (function, method, class) found in source code.


   .. py:attribute:: args
      :type:  List[str]
      :value: []



   .. py:attribute:: docstring
      :type:  Optional[str]


   .. py:attribute:: file_path
      :type:  str


   .. py:attribute:: is_cli_handler
      :type:  bool
      :value: False



   .. py:attribute:: line_number
      :type:  int


   .. py:attribute:: name
      :type:  str


   .. py:attribute:: parent_class
      :type:  Optional[str]
      :value: None



   .. py:attribute:: type
      :type:  str


.. py:class:: CoverageReport

   Docstring coverage analysis report.


   .. py:attribute:: coverage_percentage
      :type:  float


   .. py:attribute:: documented_definitions
      :type:  int


   .. py:attribute:: total_definitions
      :type:  int


   .. py:attribute:: undocumented_definitions
      :type:  int


   .. py:attribute:: undocumented_items
      :type:  List[CodeDefinition]
      :value: []



.. py:class:: EnhancedCLIValidator(config_file = None)

   Enhanced CLI documentation validation engine with AST analysis, multi-format parsing, and CI/CD integration.

   Initialize the enhanced validator with all components.


   .. py:method:: clear_cache()

      Clear validation cache and return statistics.

      :returns: Dictionary with cache clearing statistics



   .. py:method:: generate_ci_report(results)

      Generate JSON output for CI/CD integration with GitHub Actions annotations.

      :param results: ValidationReport to convert

      :returns: Dictionary suitable for JSON serialization



   .. py:method:: get_performance_stats()

      Get comprehensive performance statistics.

      :returns: Dictionary with performance metrics



   .. py:method:: validate_comprehensive(target_paths, mode = 'comprehensive')

      Complete validation using all available tools with caching and performance optimization.

      :param target_paths: List of paths to validate (files or directories)
      :param mode: Validation mode (quick, comprehensive, strict)

      :returns: ValidationReport with comprehensive results



   .. py:method:: validate_incremental(base_branch = 'main')

      Performance-optimized validation for CI/CD workflows using Git integration.

      :param base_branch: Base branch to compare against for changed files

      :returns: ValidationReport with results for changed files only



   .. py:method:: validate_staged_files()

      Validate only staged files for pre-commit hook integration.

      :returns: ValidationReport with results for staged files only



   .. py:attribute:: ast_analyzer


   .. py:attribute:: cache_manager


   .. py:attribute:: config


   .. py:attribute:: docstring_parser


   .. py:attribute:: git_integration


   .. py:attribute:: logger


   .. py:attribute:: module_validator


   .. py:attribute:: performance_optimizer


   .. py:attribute:: style_validator


.. py:class:: ExceptionInfo

   Represents exception information from docstring.


   .. py:attribute:: description
      :type:  Optional[str]
      :value: None



   .. py:attribute:: type_name
      :type:  str


.. py:class:: FileValidationResult

   Results of validating a single file.


   .. py:attribute:: coverage_percentage
      :type:  float


   .. py:attribute:: definitions
      :type:  List[CodeDefinition]


   .. py:attribute:: documented_count
      :type:  int


   .. py:attribute:: execution_time
      :type:  float


   .. py:attribute:: file_path
      :type:  str


   .. py:attribute:: issues
      :type:  List[ValidationIssue]


   .. py:attribute:: metadata
      :type:  Dict[str, Any]


   .. py:attribute:: total_definitions
      :type:  int


.. py:class:: FunctionSignature

   Represents a function signature extracted from AST.


   .. py:attribute:: parameters
      :type:  List[Parameter]


   .. py:attribute:: return_type
      :type:  Optional[str]
      :value: None



.. py:class:: GitIntegration

   Git integration for incremental validation and change detection.

   Provides functionality to detect changed files, staged files, and compare
   against different branches for efficient incremental validation.


   .. py:method:: get_changed_files(base_branch = 'main', include_untracked = False)

      Get list of Python files changed since base branch.

      :param base_branch: Base branch to compare against
      :param include_untracked: Whether to include untracked files

      :returns: List of changed Python file paths



   .. py:method:: get_staged_files()

      Get list of staged Python files for pre-commit validation.

      :returns: List of staged Python file paths



   .. py:method:: is_file_ignored(file_path)

      Check if a file is ignored by git.

      :param file_path: Path to check

      :returns: True if file is ignored by git



   .. py:attribute:: logger


.. py:class:: MultiFormatDocstringParser

   Multi-format docstring parser with structured validation capabilities.

   Supports Google, NumPy, ReST, and Epydoc docstring formats using the
   docstring_parser library with graceful fallback for unsupported formats.


   .. py:method:: parse_docstring(docstring_text, format_hint = 'auto')

      Parse a docstring into structured components.

      :param docstring_text: Raw docstring text to parse
      :param format_hint: Format hint ('google', 'numpy', 'sphinx', 'epydoc', 'auto')

      :returns: Structured representation of the docstring
      :rtype: ParsedDocstring



   .. py:method:: validate_structure(parsed_docstring, function_signature)

      Validate parsed docstring structure against function signature.

      :param parsed_docstring: Parsed docstring components
      :param function_signature: Function signature from AST analysis

      :returns: Validation results with suggestions
      :rtype: StructureValidation



   .. py:attribute:: logger


   .. py:attribute:: supported_formats
      :value: ['google', 'numpy', 'sphinx', 'epydoc', 'auto']



.. py:class:: OutputFormat

   Bases: :py:obj:`enum.Enum`


   Generic enumeration.

   Derive from this class to define new enumerations.


   .. py:attribute:: CONSOLE
      :value: 'console'



   .. py:attribute:: GITHUB
      :value: 'github'



   .. py:attribute:: JSON
      :value: 'json'



.. py:class:: Parameter

   Represents a function parameter from docstring or signature.


   .. py:attribute:: default_value
      :type:  Optional[str]
      :value: None



   .. py:attribute:: description
      :type:  Optional[str]
      :value: None



   .. py:attribute:: is_optional
      :type:  bool
      :value: False



   .. py:attribute:: name
      :type:  str


   .. py:attribute:: type_hint
      :type:  Optional[str]
      :value: None



.. py:class:: ParsedDocstring

   Represents a parsed docstring with structured components.


   .. py:attribute:: description
      :type:  str


   .. py:attribute:: detected_format
      :type:  str


   .. py:attribute:: examples
      :type:  List[str]
      :value: []



   .. py:attribute:: parameters
      :type:  List[Parameter]
      :value: []



   .. py:attribute:: parsing_errors
      :type:  List[str]
      :value: []



   .. py:attribute:: parsing_successful
      :type:  bool
      :value: True



   .. py:attribute:: raises
      :type:  List[ExceptionInfo]
      :value: []



   .. py:attribute:: raw_text
      :type:  str


   .. py:attribute:: returns
      :type:  Optional[ReturnInfo]
      :value: None



   .. py:attribute:: summary
      :type:  str


.. py:class:: PerformanceOptimizer(config)

   Performance optimization features including parallel processing and progress reporting.

   Provides multi-threaded validation, memory optimization, and progress tracking
   for efficient processing of large codebases.


   .. py:method:: optimize_memory_usage()

      Optimize memory usage by cleaning up caches and forcing garbage collection.



   .. py:method:: validate_files_parallel(file_paths, validation_func)

      Validate multiple files in parallel using thread pool.

      :param file_paths: List of file paths to validate
      :param validation_func: Function to call for each file

      :returns: List of validation results



   .. py:attribute:: config


   .. py:attribute:: logger


.. py:class:: ReturnInfo

   Represents return value information from docstring.


   .. py:attribute:: description
      :type:  Optional[str]
      :value: None



   .. py:attribute:: type_hint
      :type:  Optional[str]
      :value: None



.. py:class:: StructureValidation

   Results of structured docstring validation.


   .. py:attribute:: extra_parameters
      :type:  List[str]
      :value: []



   .. py:attribute:: format_issues
      :type:  List[str]
      :value: []



   .. py:attribute:: is_valid
      :type:  bool


   .. py:attribute:: missing_parameters
      :type:  List[str]
      :value: []



   .. py:attribute:: missing_return_doc
      :type:  bool
      :value: False



   .. py:attribute:: suggestions
      :type:  List[str]
      :value: []



   .. py:attribute:: type_mismatches
      :type:  List[str]
      :value: []



   .. py:attribute:: undocumented_exceptions
      :type:  List[str]
      :value: []



.. py:class:: ValidationConfig(config_file = None)

   Configuration management with multiple source support.

   Supports loading configuration from:
   - .cli-validation.yml
   - pyproject.toml
   - Environment variables
   - Direct initialization


   .. py:method:: get_config_summary()

      Get a summary of current configuration.



   .. py:attribute:: cache_dir
      :value: '.validation-cache'



   .. py:attribute:: cache_enabled
      :value: True



   .. py:attribute:: cache_max_size_mb
      :value: 100



   .. py:attribute:: cache_ttl_hours
      :value: 24



   .. py:attribute:: coverage_threshold
      :value: 95.0



   .. py:attribute:: docstring_formats
      :value: ['google', 'numpy', 'sphinx']



   .. py:attribute:: fail_on_warnings
      :value: False



   .. py:attribute:: github_annotations
      :value: True



   .. py:attribute:: incremental_mode
      :value: True



   .. py:attribute:: json_output
      :value: False



   .. py:attribute:: logger


   .. py:attribute:: parallel_workers
      :value: 4



   .. py:attribute:: pre_commit_quick_mode
      :value: True



   .. py:attribute:: require_examples
      :value: True



   .. py:attribute:: staged_files_only
      :value: True



   .. py:attribute:: timeout_seconds
      :value: 300



   .. py:attribute:: validate_parameter_types
      :value: True



.. py:class:: ValidationIssue

   Represents a validation issue found during analysis (updated structure).


   .. py:attribute:: location
      :type:  ValidationLocation


   .. py:attribute:: message
      :type:  str


   .. py:attribute:: severity
      :type:  str


   .. py:attribute:: suggestions
      :type:  List[str]
      :value: []



   .. py:attribute:: type
      :type:  str


.. py:class:: ValidationLocation

   Represents the location of a validation issue.


   .. py:attribute:: column_number
      :type:  Optional[int]
      :value: None



   .. py:attribute:: file_path
      :type:  str


   .. py:attribute:: line_number
      :type:  int


.. py:class:: ValidationMode

   Bases: :py:obj:`enum.Enum`


   Generic enumeration.

   Derive from this class to define new enumerations.


   .. py:attribute:: COMPREHENSIVE
      :value: 'comprehensive'



   .. py:attribute:: QUICK
      :value: 'quick'



   .. py:attribute:: STRICT
      :value: 'strict'



.. py:class:: ValidationReport

   Comprehensive validation report.


   .. py:attribute:: coverage_percentage
      :type:  float


   .. py:attribute:: execution_time
      :type:  float


   .. py:attribute:: issues
      :type:  List[ValidationIssue]


   .. py:attribute:: metadata
      :type:  Dict[str, Any]


   .. py:attribute:: total_definitions
      :type:  int


   .. py:attribute:: total_files
      :type:  int


.. py:function:: main()

   Main CLI interface for enhanced CLI documentation validation.


.. py:function:: setup_logging(verbose = False)

   Setup logging configuration.


.. py:data:: DOCSTRING_PARSER_AVAILABLE
   :value: True


.. py:data:: PYDANTIC_SETTINGS_AVAILABLE
   :value: False


.. py:data:: PYDANTIC_SETTINGS_AVAILABLE
   :value: True


.. py:data:: TOML_SUPPORT
   :value: False


.. py:data:: TOML_SUPPORT
   :value: True


