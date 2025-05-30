performance_analysis
====================

.. py:module:: performance_analysis

.. autoapi-nested-parse::

   Performance Analysis Script for VCF Analysis Agent

   This script performs comprehensive performance analysis including:
   - Memory usage profiling
   - Database query performance
   - AI model response times
   - Code hotspot identification
   - Resource utilization analysis



Classes
-------

.. autoapisummary::

   performance_analysis.PerformanceAnalyzer


Functions
---------

.. autoapisummary::

   performance_analysis.main


Module Contents
---------------

.. py:class:: PerformanceAnalyzer(output_dir = 'performance_reports')

   Comprehensive performance analyzer for VCF Agent.


   .. py:method:: analyze_ai_performance()

      Analyze AI model performance.



   .. py:method:: analyze_code_hotspots()

      Analyze code performance hotspots using cProfile.



   .. py:method:: analyze_database_performance()

      Analyze database operation performance.



   .. py:method:: analyze_memory_usage()

      Analyze memory usage patterns.



   .. py:method:: analyze_resource_utilization()

      Analyze system resource utilization.



   .. py:method:: generate_optimization_recommendations()

      Generate optimization recommendations based on analysis.



   .. py:method:: print_summary()

      Print analysis summary.



   .. py:method:: run_full_analysis()

      Run complete performance analysis.



   .. py:method:: save_results(filename = None)

      Save analysis results to file.



   .. py:attribute:: output_dir


   .. py:attribute:: results


   .. py:attribute:: start_time


.. py:function:: main()

   Main function to run performance analysis.


