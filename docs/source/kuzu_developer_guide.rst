\
Kuzu Developer Guide
====================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Overview
--------

Kuzu is an embedded, high-performance property graph database management system (GDBMS).
In the VCF Agent project, Kuzu serves as a critical component for storing and querying
relationships between genetic variants and samples. Its primary role is to provide
rich graph context for variant IDs that are initially retrieved through vector
searches in LanceDB. This allows the agent to understand how variants are observed
across different samples, their zygosity, and other relational information.

Kuzu's embedded nature means it runs within the agent's process, simplifying
deployment and reducing latency for graph queries. It supports the Cypher query
language, making it relatively easy to define schemas and query graph data.

Core Module: ``src.vcf_agent.graph_integration``
-------------------------------------------------

The primary interaction with Kuzu is managed through the
``src.vcf_agent.graph_integration`` module. This module provides a suite of
functions for database connection, schema management, data ingestion, and querying.

.. automodule:: src.vcf_agent.graph_integration
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

While ``automodule`` above lists all members, some key functions include:

*   ``get_kuzu_db_connection(db_path: str, read_only: bool)``: Establishes a connection to the Kuzu database (on-disk or in-memory).
*   ``create_schema(conn: kuzu.Connection)``: Defines and creates the necessary graph schema.
*   ``add_variant(conn: kuzu.Connection, variant_data: Dict[str, Any])``: Adds a variant node.
*   ``add_sample(conn: kuzu.Connection, sample_data: Dict[str, Any])``: Adds a sample node.
*   ``link_variant_to_sample(conn: kuzu.Connection, variant_id: str, sample_id: str, properties: Dict[str, Any])``: Links a variant to a sample.
*   ``execute_query(conn: kuzu.Connection, cypher_query: str, params: Optional[Dict[str, Any]])``: Executes an arbitrary Cypher query.
*   ``get_variant_context(conn: kuzu.Connection, variant_ids: List[str])``: Retrieves graph context for a list of variant IDs (crucial for post-LanceDB enrichment).

Schema Details
--------------

The Kuzu graph database within the VCF Agent uses the following schema:

**Nodes:**

1.  **Variant**
    *   Properties:
        *   ``variant_id`` (STRING, PRIMARY KEY): A unique identifier for the variant (e.g., "chr1-12345-A-T").
        *   ``chrom`` (STRING): Chromosome (e.g., "1", "X").
        *   ``pos`` (INT64): Position on the chromosome.
        *   ``ref`` (STRING): Reference allele.
        *   ``alt`` (STRING): Alternative allele.
        *   ``rs_id`` (STRING, optional): dbSNP Reference SNP identifier (e.g., "rs12345").

2.  **Sample**
    *   Properties:
        *   ``sample_id`` (STRING, PRIMARY KEY): A unique identifier for the sample.

**Relationships:**

1.  **ObservedIn**
    *   Direction: (Variant) -[:ObservedIn]-> (Sample)
    *   Properties:
        *   ``zygosity`` (STRING): The zygosity of the variant in the sample (e.g., "HET" for heterozygous, "HOM" for homozygous).

Usage Pattern with LanceDB
--------------------------

The typical workflow involving Kuzu is as follows:

1.  **VCF Parsing**: Variants and sample information are parsed from a VCF file.
2.  **LanceDB Ingestion**: Variant information (potentially embeddings or key features) is ingested into LanceDB for fast vector similarity searches.
3.  **Kuzu Ingestion**: Variant nodes, sample nodes, and their ObservedIn relationships (including zygosity) are ingested into Kuzu.
4.  **Vector Search (LanceDB)**: The user queries the agent, which performs a semantic or similarity search in LanceDB, retrieving a list of relevant ``variant_id``s.
5.  **Graph Context Enrichment (Kuzu)**: These ``variant_id``s are passed to the ``get_variant_context`` function in the ``graph_integration`` module. Kuzu then efficiently queries the graph to find all samples linked to these variants and the properties of those links (e.g., zygosity).
6.  **Result Aggregation**: The agent combines the information from LanceDB and Kuzu to provide a comprehensive answer to the user.

This two-database approach leverages the strengths of both: LanceDB for fast vector search and Kuzu for complex relational queries and graph traversal.

Example (Conceptual):
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Assume 'lancedb_conn' and 'kuzu_conn' are established connections
   # And 'user_query_embedding' is the embedding of the user's query

   # 1. Search in LanceDB
   # relevant_variant_ids = lancedb_module.search_variants(lancedb_conn, user_query_embedding)
   # Example: relevant_variant_ids = ["chr1-100-A-T", "chrX-5000-C-G"]

   # 2. Get graph context from Kuzu
   # variant_contexts = graph_integration.get_variant_context(kuzu_conn, relevant_variant_ids)

   # 3. Process and present results
   # for v_id, context in variant_contexts.items():
   #     print(f"Variant: {v_id}")
   #     for sample_info in context.get('samples', []):
   #         print(f"  Observed in: {sample_info['sample_id']}, Zygosity: {sample_info['zygosity']}")


This guide provides a developer-focused overview of Kuzu integration. For
more detailed API documentation, refer to the auto-generated module documentation
and the docstrings within ``src/vcf_agent/graph_integration.py``. 