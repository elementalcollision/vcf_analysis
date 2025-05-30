validate_phase5_2
=================

.. py:module:: validate_phase5_2

.. autoapi-nested-parse::

   Phase 5.2 Production Validation Script
   =====================================

   Comprehensive validation of the hybrid Apache Iggy + Kafka streaming
   architecture with intelligent routing, circuit breaker patterns, and
   exactly-once delivery semantics.

   Features Validated:
   - Dual-platform coordination and intelligent routing
   - Circuit breaker patterns and automatic failover
   - Message deduplication and exactly-once semantics
   - Performance monitoring and health-based decisions
   - End-to-end processing with synthetic VCF data

   Usage:
       python scripts/validate_phase5_2.py



Attributes
----------

.. autoapisummary::

   validate_phase5_2.logger


Classes
-------

.. autoapisummary::

   validate_phase5_2.Phase5_2Validator


Functions
---------

.. autoapisummary::

   validate_phase5_2.main


Module Contents
---------------

.. py:class:: Phase5_2Validator

   Comprehensive validator for Phase 5.2 dual-platform architecture.

   Demonstrates production-ready patterns including:
   - Intelligent routing based on platform health
   - Circuit breaker patterns for automatic failover
   - Message deduplication for exactly-once semantics
   - Performance monitoring and alerting

   Initialize validator with test configuration.


   .. py:method:: create_synthetic_variants(count = 1000)

      Create synthetic VCF variants for testing.



   .. py:method:: run_validation_suite()
      :async:


      Run complete Phase 5.2 validation suite.



   .. py:method:: validate_circuit_breaker_patterns()
      :async:


      Validate circuit breaker implementation for platform health management.



   .. py:method:: validate_end_to_end_processing()
      :async:


      Validate end-to-end dual-platform processing.



   .. py:method:: validate_intelligent_routing()
      :async:


      Validate intelligent platform routing based on health metrics.



   .. py:method:: validate_message_deduplication()
      :async:


      Validate exactly-once semantics with message deduplication.



   .. py:method:: validate_performance_monitoring()
      :async:


      Validate performance monitoring and platform health tracking.



   .. py:attribute:: config


   .. py:attribute:: results


.. py:function:: main()
   :async:


   Main validation entry point.


.. py:data:: logger

