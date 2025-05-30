demo_phase5_iggy
================

.. py:module:: demo_phase5_iggy

.. autoapi-nested-parse::

   Phase 5.1 Apache Iggy Integration Demo
   =====================================

   Demonstrates the Apache Iggy integration for ultra-high-performance
   genomic variant streaming with the hybrid architecture.

   Prerequisites:
   1. Apache Iggy server running: docker run --rm -p 8080:8080 -p 3000:3000 -p 8090:8090 iggyrs/iggy:0.4.21
   2. Dependencies installed: pip install iggy-py>=0.4.0

   Usage:
       python scripts/demo_phase5_iggy.py --variants 1000 --show-performance



Attributes
----------

.. autoapisummary::

   demo_phase5_iggy.exit_code
   demo_phase5_iggy.logger


Functions
---------

.. autoapisummary::

   demo_phase5_iggy.check_iggy_availability
   demo_phase5_iggy.create_test_variants
   demo_phase5_iggy.demo_iggy_processor
   demo_phase5_iggy.demo_performance_comparison
   demo_phase5_iggy.demo_streaming_coordinator
   demo_phase5_iggy.main
   demo_phase5_iggy.print_performance_summary


Module Contents
---------------

.. py:function:: check_iggy_availability()
   :async:


   Check if Iggy server is available.


.. py:function:: create_test_variants(count)

   Create test VCF variants for demonstration.


.. py:function:: demo_iggy_processor(variants)
   :async:


   Demonstrate direct Iggy processor usage.


.. py:function:: demo_performance_comparison(variants)
   :async:


   Compare performance between different processing methods.


.. py:function:: demo_streaming_coordinator(variants)
   :async:


   Demonstrate the hybrid streaming coordinator.


.. py:function:: main()
   :async:


   Main demonstration function.


.. py:function:: print_performance_summary(results, variant_count)

   Print performance summary.


.. py:data:: exit_code

.. py:data:: logger

