generate_golden_files
=====================

.. py:module:: generate_golden_files

.. autoapi-nested-parse::

   Golden File Generation Script

   Generates golden files by running validated VCF operations on known-good test data.
   This script should be run manually after verifying that the operations produce
   correct outputs.



Attributes
----------

.. autoapisummary::

   generate_golden_files.project_root


Functions
---------

.. autoapisummary::

   generate_golden_files.generate_analysis_golden_files
   generate_golden_files.generate_bcftools_golden_files
   generate_golden_files.generate_comparison_golden_files
   generate_golden_files.generate_validation_golden_files
   generate_golden_files.main


Module Contents
---------------

.. py:function:: generate_analysis_golden_files()

   Generate golden files for agent analysis operations.


.. py:function:: generate_bcftools_golden_files()

   Generate golden files for bcftools operations.


.. py:function:: generate_comparison_golden_files()

   Generate golden files for VCF comparison operations.


.. py:function:: generate_validation_golden_files()

   Generate golden files for VCF validation operations.


.. py:function:: main()

   Main function to generate all golden files.


.. py:data:: project_root

