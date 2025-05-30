test-documentation
==================

.. py:module:: test-documentation

.. autoapi-nested-parse::

   Documentation Testing Framework for CLI Enhanced Validation Engine
   Comprehensive testing suite for documentation quality, links, performance, and content accuracy.



Attributes
----------

.. autoapisummary::

   test-documentation.logger


Classes
-------

.. autoapisummary::

   test-documentation.DocumentationTester


Functions
---------

.. autoapisummary::

   test-documentation.main


Module Contents
---------------

.. py:class:: DocumentationTester(base_dir = None)

   Main documentation testing framework.


   .. py:method:: run_all_tests()

      Run all documentation tests and return comprehensive results.



   .. py:method:: test_accessibility()

      Test accessibility compliance.



   .. py:method:: test_content_validation()

      Validate content structure and completeness.



   .. py:method:: test_link_validation()

      Validate all internal and external links.



   .. py:method:: test_mkdocs_build()

      Test MkDocs builds successfully without errors.



   .. py:method:: test_performance()

      Test page load performance and build times.



   .. py:method:: test_seo()

      Test SEO optimization.



   .. py:method:: test_sphinx_build()

      Test Sphinx API documentation generates correctly.



   .. py:attribute:: base_dir


   .. py:attribute:: base_url
      :value: 'http://localhost:8000'



   .. py:attribute:: docs_dir


   .. py:attribute:: external_link_timeout
      :value: 10



   .. py:attribute:: max_concurrent_links
      :value: 10



   .. py:attribute:: performance_thresholds


   .. py:attribute:: site_dir


   .. py:attribute:: sphinx_build_dir


.. py:function:: main()

   Main test runner.


.. py:data:: logger

