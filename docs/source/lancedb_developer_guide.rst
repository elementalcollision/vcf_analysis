LanceDB Developer Guide
=======================

.. contents:: Table of Contents
   :depth: 2
   :local:

Overview
--------
This guide provides comprehensive documentation for the LanceDB integration in the VCF Analysis Agent, including API usage, CLI commands, concurrency model, security considerations, and schema details.

API Reference
-------------
.. automodule:: vcf_agent.lancedb_integration
   :members:
   :undoc-members:
   :show-inheritance:

CLI Usage Guide
---------------

The VCF Analysis Agent CLI provides several commands for interacting with LanceDB. Below are the main commands, their arguments, and usage examples:

**Initialize LanceDB Table**

.. code-block:: bash

    python -m vcf_agent.cli init-lancedb --db_path ./lancedb --table_name variants

Initializes a LanceDB table for storing variants. If the table already exists, it will be opened.

**Add a Variant**

.. code-block:: bash

    python -m vcf_agent.cli add-variant --db_path ./lancedb --table_name variants \
        --variant_id chr1-123-A-G --chrom 1 --pos 123 --ref A --alt G \
        --embedding 0.1,0.2,... --clinical_significance Pathogenic

Adds a new variant record to the LanceDB table. The embedding should be a comma-separated list of floats.

**Search by Embedding**

.. code-block:: bash

    python -m vcf_agent.cli search-embedding --db_path ./lancedb --table_name variants \
        --embedding 0.1,0.2,... --limit 5 --filter_sql "chrom = '1'"

Searches for variants similar to the provided embedding. Optionally filter results using SQL-like syntax.

**Update a Variant**

.. code-block:: bash

    python -m vcf_agent.cli update-variant --db_path ./lancedb --table_name variants \
        --variant_id chr1-123-A-G --updates '{"clinical_significance": "Benign"}'

Updates fields for a specific variant. The updates argument must be a valid JSON string.

**Delete Variants**

.. code-block:: bash

    python -m vcf_agent.cli delete-variants --db_path ./lancedb --table_name variants \
        --filter_sql "chrom = '1' AND pos < 1000"

Deletes variants matching the filter.

**Create Scalar Index**

.. code-block:: bash

    python -m vcf_agent.cli create-lancedb-index --db_path ./lancedb --table_name variants \
        --column chrom --index_type BTREE --replace

Creates a scalar index on a column to improve filter performance.

**Filter by Metadata**

.. code-block:: bash

    python -m vcf_agent.cli filter-lancedb --db_path ./lancedb --table_name variants \
        --filter_sql "chrom = '1' AND pos > 1000" --select_columns variant_id,chrom,pos --limit 10

Returns variants matching the metadata filter.

Data Flow Overview
------------------

.. mermaid::

   graph LR
     CLI[CLI/API Call] --> LDB[LanceDB Table]
     LDB -->|Embedding Search| Kuzu[Kuzu Graph Enrichment]
     Kuzu -->|Context| LDB
     LDB --> Result[Result to User]

Schema Definition
-----------------

The LanceDB table for variants uses the following schema (see also the `Variant` class in the API Reference):

+------------------------+---------------------+-----------------------------------------------+
| Field                  | Type                | Description                                   |
+========================+=====================+===============================================+
| variant_id             | str                 | Unique identifier (e.g., 'chr1-123-A-G')      |
+------------------------+---------------------+-----------------------------------------------+
| embedding              | Vector(1024)        | 1024-dimensional embedding vector             |
+------------------------+---------------------+-----------------------------------------------+
| chrom                  | str                 | Chromosome (e.g., '1', 'X')                   |
+------------------------+---------------------+-----------------------------------------------+
| pos                    | int                 | Genomic position                              |
+------------------------+---------------------+-----------------------------------------------+
| ref                    | str                 | Reference allele                              |
+------------------------+---------------------+-----------------------------------------------+
| alt                    | str                 | Alternate allele                              |
+------------------------+---------------------+-----------------------------------------------+
| clinical_significance  | Optional[str]       | Clinical significance annotation              |
+------------------------+---------------------+-----------------------------------------------+

**Example Record:**

.. code-block:: json

    {
        "variant_id": "chr1-123-A-G",
        "embedding": [0.1, 0.2, ...],
        "chrom": "1",
        "pos": 123,
        "ref": "A",
        "alt": "G",
        "clinical_significance": "Pathogenic"
    }

Variant Table Schema (ER Diagram)
---------------------------------

.. mermaid::

   erDiagram
     VARIANT {
       string variant_id PK
       vector embedding
       string chrom
       int pos
       string ref
       string alt
       string clinical_significance
     }

Concurrency Model
-----------------

.. mermaid::

   graph TD
     A[Thread 1: Add Variant] -- waits for lock --> L[lancedb_write_lock]
     B[Thread 2: Update Variant] -- waits for lock --> L
     C[Thread 3: Delete Variant] -- waits for lock --> L
     L -- grants lock --> A
     L -- grants lock --> B
     L -- grants lock --> C

The LanceDB integration in the VCF Analysis Agent is designed to be thread-safe for local file-based databases. This is achieved using a global `threading.Lock` called `lancedb_write_lock`.

**Why is a lock needed?**
- LanceDB uses local file storage by default. Concurrent write operations (add, update, delete, index creation) from multiple threads can cause race conditions or data corruption.
- The global lock serializes all write-related operations, ensuring only one thread can modify the database at a time within a single process.

**How is the lock used?**
- All write operations (`add_variants`, `update_variant`, `delete_variants`, `create_scalar_index`) are wrapped in a `with lancedb_write_lock:` block.
- This ensures that only one thread can perform a write at any given time.

**Example:**
.. code-block:: python

    lancedb_write_lock = threading.Lock()

    def add_variants(table, variants):
        with lancedb_write_lock:
            table.add(variants)

**Limitations:**
- This locking strategy only protects against concurrent writes within a single Python process.
- For multi-process or distributed deployments, a more robust distributed locking mechanism or a database designed for high concurrency (such as LanceDB Cloud) is recommended.

**Read Operations:**
- Read operations (e.g., search, filter) are not locked, as they do not modify the database and are safe to run concurrently.

**Summary:**
- The use of `threading.Lock` ensures safe, reliable operation for local development and most single-node deployments.
- Developers extending the system for distributed or high-concurrency use cases should consider additional synchronization mechanisms. 

Security Notes
--------------

The LanceDB integration follows the Auditable Security Framework to ensure data confidentiality, integrity, and auditability. Key practices include:

- **Data Encryption:**
  - Use filesystem-level encryption for LanceDB storage (e.g., LUKS, BitLocker, APFS encryption).
  - If using LanceDB Enterprise, enable built-in AES-256 encryption and integrate with a KMS if available.
- **Access Control:**
  - Restrict file and directory permissions so only the VCF Agent service account can access LanceDB data.
  - Run the agent under a dedicated, unprivileged service account.
- **Audit Logging:**
  - All critical operations (connect, add, search, update, delete, index creation, filter) are logged with timestamps and status.
  - Sensitive data in logs is masked (e.g., SQL literals, variant IDs if needed).
  - Enable OS-level audit logging for file access in the LanceDB directory.
- **Data Masking & Anonymization:**
  - Mask or anonymize PII/PHI in VCF data and embeddings where feasible.
  - Mask sensitive fields in logs and outputs unless required for authorized use.
- **Secure Configuration:**
  - Harden the OS and restrict network access to LanceDB files.
  - Use secure defaults for all configuration parameters.
- **Vulnerability Management:**
  - Regularly scan dependencies and the OS for vulnerabilities.
  - Maintain an SBOM and promptly patch critical issues.
- **Incident Response:**
  - Maintain and test an incident response plan for LanceDB and VCF data.
- **Compliance:**
  - Align practices with relevant regulations (e.g., HIPAA, GDPR) and internal governance.

For full details, see the "LanceDB Auditable Security Framework" decision document. 