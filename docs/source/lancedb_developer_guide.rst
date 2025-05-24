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
.. note::
   This section will provide usage instructions and examples for all LanceDB-related CLI commands.

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
.. note::
   This section will summarize key points from the Auditable Security Framework relevant to developers, including logging and data handling best practices.

Schema Definition
-----------------
.. note::
   This section will detail the Variant schema and LanceDB table structure. 