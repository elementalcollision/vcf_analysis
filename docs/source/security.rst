Security
========

This document outlines the security considerations and best practices for the VCF Analysis Agent.

General Principles
------------------
- **Least Privilege**: Processes and components should run with the minimum necessary permissions.
- **Defense in Depth**: Employ multiple layers of security controls.
- **Secure Defaults**: Configure components with security in mind from the outset.
- **Regular Updates**: Keep all dependencies, OS, and tools patched and up-to-date.

Credential Management
---------------------

Secure management of API keys and other credentials is critical, especially for accessing Large Language Models (LLMs) and other cloud services.

- **LLM API Keys**: 
    - Refer to the main project `README.md` for detailed instructions on providing credentials via `.env` files or JSON configuration files.
    - Avoid hardcoding credentials in source code.
    - Use environment variables or secure configuration files for production deployments.
    - Restrict access to credential files using appropriate filesystem permissions.

- **Database Credentials (if applicable)**:
    - If Kuzu or LanceDB were to be deployed as remote services requiring authentication (currently they are used locally), connection strings and credentials must be managed securely, similar to LLM API keys.

Data Handling for VCF Files
---------------------------
Variant Call Format (VCF) files can contain sensitive genomic information and Personally Identifiable Information (PII).

- **Access Control**: Implement strict access controls on VCF files and any derived data (e.g., LanceDB/Kuzu databases). Only authorized personnel and processes should have access.
- **Encryption**: 
    - **At Rest**: Genomic data and database files stored locally should be protected by full-disk encryption (e.g., LUKS, BitLocker, APFS encryption) on the host system.
    - **In Transit**: If the agent communicates with remote services or databases (not current for LanceDB/Kuzu core usage but applicable for LLM APIs), ensure TLS/SSL is used for all communications.
- **Anonymization/Pseudonymization**: If possible and appropriate for the use case, consider anonymizing or pseudonymizing VCF data before processing, especially in research settings.
- **Data Minimization**: Only collect and retain the data necessary for the intended analysis.

Database Security
-----------------

Specific security measures for the integrated LanceDB and Kuzu databases are crucial.

### LanceDB Security
LanceDB is used locally, so its security heavily relies on host and application-level measures.

- **Filesystem Permissions**: As detailed in the :doc:`lancedb_developer_guide`, the directory storing LanceDB data (e.g., `./lancedb`) must have strict ACLs. Only the VCF Agent service account/user should have read/write permissions.
- **Input Validation**: The agent implements checks to prevent basic SQL injection attempts in `filter_sql` queries for LanceDB. See :doc:`lancedb_developer_guide` for more details.
- **Auditing**: Filesystem-level auditing for the LanceDB data directory is recommended. See :doc:`audit`.

### Kuzu Security
Kuzu is also used locally, and similar host/application-level security measures apply.

- **Filesystem Permissions**: The Kuzu database directory (e.g., `./kuzu_db`) requires strict ACLs, granting access only to the VCF Agent process.
- **Input Validation**: While Cypher queries for Kuzu are typically constructed internally, any user-influenced parameters should be handled carefully to prevent injection if applicable (though less common than SQLi).
- **Auditing**: Filesystem-level auditing for the Kuzu data directory is recommended. See :doc:`audit`.

Application Security
--------------------

- **Dependency Management**: 
    - Use `uv` or `pip` to manage dependencies as defined in `pyproject.toml` and `requirements.txt`.
    - Regularly update dependencies to patch known vulnerabilities.
- **Vulnerability Scanning**: 
    - Employ tools like `pip-audit` (as integrated via `make audit-dependencies`) or commercial scanners (Snyk, Trivy) to check for vulnerabilities in Python packages.
    - Integrate these scans into the CI/CD pipeline.
- **Software Bill of Materials (SBOM)**: 
    - Generate and maintain an SBOM to track all software components and their versions. This aids in quickly identifying impacts of newly discovered vulnerabilities.

Secure Development Practices
----------------------------
- Follow secure coding guidelines.
- Perform code reviews with a focus on security aspects.
- Implement comprehensive testing, including E2E tests for security-related error conditions.

Incident Response
-----------------
- Have a plan for responding to potential security incidents, including identification, containment, eradication, recovery, and lessons learned. 