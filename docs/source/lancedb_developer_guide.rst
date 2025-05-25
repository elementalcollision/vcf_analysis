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

The LanceDB integration in the VCF Analysis Agent is designed with security in mind. For comprehensive details on security best practices, credential management, data handling, and auditing, please refer to the main project documentation:

- :doc:`security`
- :doc:`audit`

Key security considerations specific to LanceDB local file-based usage include:

- **Filesystem Permissions & ACLs**: Strict access controls must be applied to the LanceDB data directory (e.g., `./lancedb`). Only the agent process should have read/write access.
- **Input Validation**: The agent includes basic checks against SQL injection in `filter_sql` clauses. Always sanitize or validate user-provided inputs if they are to be part of database queries.
- **Data Encryption**: Utilize filesystem-level encryption (e.g., LUKS, BitLocker, APFS encryption) for the directory storing LanceDB data to protect data at rest.
- **Logging & Auditing**: Sensitive query parameters are masked in logs. OS-level filesystem auditing should be enabled for the LanceDB data directory.

Refer to the main :doc:`security` and :doc:`audit` documents for the complete framework.

Performance Considerations
------------------------

The LanceDB integration in the VCF Analysis Agent is designed to be efficient and scalable. Key considerations include:

- **Filter Performance**: The agent uses scalar indexes and efficient query processing to quickly retrieve relevant variants.
- **Memory Usage**: The agent manages memory efficiently by only loading necessary data into memory.
- **Concurrency**: The agent is designed to handle multiple concurrent operations, ensuring smooth performance even under high load.

Vulnerability Scanning with pip-audit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To ensure that your Python dependencies are free from known vulnerabilities, the VCF Agent project uses `pip-audit` for automated vulnerability scanning.

Purpose
^^^^^^^
- Identify and remediate known vulnerabilities in Python dependencies before deployment.
- Support compliance and secure development lifecycle requirements.

CI/CD Integration
^^^^^^^^^^^^^^^^^
- The CI/CD pipeline (see `kestra/flows/python-pip-audit.yml`) runs `pip-audit` on every build.
- If vulnerabilities are found, the build fails and a report is generated.
- A CycloneDX SBOM is also generated and uploaded as a build artifact for compliance and supply chain security.

Local Usage
^^^^^^^^^^^
- Developers can run `pip-audit` locally to check for vulnerabilities before committing or pushing code.
- Install `pip-audit` (if not already available):

  .. code-block:: bash

     pip install pip-audit

- Run a vulnerability scan on your current environment or requirements file:

  .. code-block:: bash

     pip-audit -r requirements.txt

- To generate a CycloneDX SBOM as well:

  .. code-block:: bash

     pip-audit -r requirements.txt -f cyclonedx-json -o sbom.json

- Review the output for any vulnerabilities and address them before merging or deploying.

Interpreting Results
^^^^^^^^^^^^^^^^^^^^
- If vulnerabilities are found, `pip-audit` will list the affected packages, versions, and CVE identifiers.
- The CI/CD build will fail if any vulnerabilities are detected, enforcing a secure baseline.
- For more details, see the [pip-audit documentation](https://pypi.org/project/pip-audit/).

Related: See the SBOM generation section for more on supply chain security.

Software Bill of Materials (SBOM) Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Software Bill of Materials (SBOM) is a comprehensive inventory of all packages and dependencies used in your project. SBOMs are critical for supply chain security, vulnerability management, and regulatory compliance.

Purpose
^^^^^^^
- Provide transparency into all software components and dependencies.
- Support vulnerability management, license compliance, and incident response.
- Meet requirements for many security frameworks and regulations.

CI/CD Integration
^^^^^^^^^^^^^^^^^
- The CI/CD pipeline (see `kestra/flows/python-pip-audit.yml`) generates a CycloneDX SBOM using both `pip-audit` and `cyclonedx-bom`.
- The SBOMs are uploaded as build artifacts (e.g., `sbom.json`, `sbom-cyclonedx.json`) for later review, audit, or compliance reporting.

Local Usage
^^^^^^^^^^^
- Developers can generate an SBOM locally for their environment or requirements file.
- Using pip-audit:

  .. code-block:: bash

     pip-audit -r requirements.txt -f cyclonedx-json -o sbom.json

- Using cyclonedx-bom:

  .. code-block:: bash

     pip install cyclonedx-bom
     cyclonedx-py requirements -i requirements.txt -o sbom-cyclonedx.json

- The resulting SBOM files can be shared with auditors, security teams, or customers as needed.

SBOM Formats
^^^^^^^^^^^^
- CycloneDX is a widely adopted, machine-readable SBOM format supported by many security and compliance tools.
- For more details, see the [CycloneDX documentation](https://cyclonedx.org/docs/).

Best Practices
^^^^^^^^^^^^^^
- Always generate and store SBOMs for each release or deployment.
- Review SBOMs for unexpected or outdated dependencies.
- Use SBOMs to support vulnerability management and compliance audits.

Filesystem ACLs and Secure Permissions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To minimize the risk of unauthorized access or modification, always apply the principle of least privilege to the LanceDB data directory and files.

Principle of Least Privilege
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Grant only the minimum permissions needed for the VCF Agent to function.
- Use groups/roles for access control, not individual users.
- Store sensitive files in dedicated, locked-down directories.
- Avoid world-readable/writable permissions (never use 777 or 'Everyone:Full Control').

Linux/macOS
^^^^^^^^^^^
- **Recommended Permissions:**

  .. code-block:: bash

     # Set owner to 'lancedbuser' and restrict access
     groupadd lancedb
     useradd -g lancedb lancedbuser
     chown -R lancedbuser:lancedb /var/lib/lancedb
     chmod 700 /var/lib/lancedb

  # For files:
     chmod 600 /var/lib/lancedb/*

- **Advanced (ACLs):**

  .. code-block:: bash

     setfacl -m u:username:rwx /var/lib/lancedb

- **Check Permissions:**

  .. code-block:: bash

     ls -ld /var/lib/lancedb
     ls -l /var/lib/lancedb

Windows
^^^^^^^
- **Recommended Permissions:**
  - Remove 'Everyone' and 'Users' groups from sensitive directories.
  - Grant 'Full Control' only to the service account or Administrators.
- **Set via PowerShell:**

  .. code-block:: powershell

     $folder = "C:\\path\\to\\lancedb"
     $acl = Get-Acl $folder
     $rule = New-Object System.Security.AccessControl.FileSystemAccessRule("LanceDBUser","FullControl","Allow")
     $acl.SetAccessRule($rule)
     Set-Acl $folder $acl

- **Check Permissions:**

  .. code-block:: bash

     icacls C:\path\to\lancedb

Review and Monitoring
^^^^^^^^^^^^^^^^^^^^^
- Schedule permission reviews (quarterly or after role changes).
- Enable auditing and monitor for suspicious or failed access attempts.
- See :doc:`audit/secure_permissions` for a detailed, step-by-step guide.

References
^^^^^^^^^^
- [OWASP File Permission Guide](https://owasp.org/www-project-web-security-testing-guide/latest/4-Web_Application_Security_Testing/02-Configuration_and_Deployment_Management_Testing/09-Test_File_Permission)
- [Microsoft File and Folder Permissions](https://docs.microsoft.com/en-us/windows/security/threat-protection/security-policy-settings/access-control-permissions)
- [Linux chmod, chown, setfacl](https://man7.org/linux/man-pages/man1/chmod.1.html) 