comprehensive_testing
=====================

.. py:module:: comprehensive_testing

.. autoapi-nested-parse::

   Comprehensive Testing Suite for VCF Analysis Agent

   This script performs comprehensive testing to validate optimizations and ensure production readiness:
   - Performance regression testing
   - Functionality validation
   - Load testing
   - Memory leak detection
   - Error handling validation



Attributes
----------

.. autoapisummary::

   comprehensive_testing.logger


Classes
-------

.. autoapisummary::

   comprehensive_testing.ComprehensiveTestSuite


Functions
---------

.. autoapisummary::

   comprehensive_testing.main
   comprehensive_testing.save_test_results


Module Contents
---------------

.. py:class:: ComprehensiveTestSuite

   Comprehensive testing suite for VCF Analysis Agent.


   .. py:method:: cleanup_test_environment()

      Clean up test environment.



   .. py:method:: generate_test_summary()

      Generate comprehensive test summary.



   .. py:method:: print_test_results()

      Print formatted test results.



   .. py:method:: run_all_tests()

      Run all comprehensive tests.



   .. py:method:: setup_test_environment()

      Set up isolated test environment.



   .. py:method:: test_error_handling()

      Test error handling and recovery.



   .. py:method:: test_functionality_validation()

      Validate core functionality still works after optimizations.



   .. py:method:: test_load_performance()

      Test performance under load.



   .. py:method:: test_memory_leaks()

      Test for memory leaks.



   .. py:method:: test_performance_regression()

      Test for performance regressions.



   .. py:attribute:: results


   .. py:attribute:: temp_dir
      :value: None



.. py:function:: main()

   Main function to run comprehensive tests.


.. py:function:: save_test_results(results)

   Save test results to file.


.. py:data:: logger

