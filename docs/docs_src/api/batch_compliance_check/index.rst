batch_compliance_check
======================

.. py:module:: batch_compliance_check

.. autoapi-nested-parse::

   Batch Compliance Checker for VCF/BCF Files
   - Scans a directory for VCF/BCF files
   - Optionally generates synthetic edge-case files
   - Runs compliance checks (bcftools, GATK, or as configured)
   - Outputs a markdown summary report



Attributes
----------

.. autoapisummary::

   batch_compliance_check.EDGE_CASES


Functions
---------

.. autoapisummary::

   batch_compliance_check.generate_color_markdown_report
   batch_compliance_check.generate_edgecases
   batch_compliance_check.generate_html_report
   batch_compliance_check.generate_markdown_report
   batch_compliance_check.is_vcf_or_bcf
   batch_compliance_check.main


Module Contents
---------------

.. py:function:: generate_color_markdown_report(results, directory)

.. py:function:: generate_edgecases(directory)

.. py:function:: generate_html_report(results, directory)

.. py:function:: generate_markdown_report(results, directory)

.. py:function:: is_vcf_or_bcf(filename)

.. py:function:: main()

.. py:data:: EDGE_CASES

