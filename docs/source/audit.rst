Auditing
========

Effective auditing is essential for monitoring the VCF Analysis Agent, ensuring security, and aiding in troubleshooting and compliance.

Logging Strategy
----------------

The VCF Analysis Agent employs a structured logging strategy using Python's built-in `logging` module. Logs are designed to capture key operations, parameters, and outcomes while being mindful of sensitive data.

- **What is Logged**:
    - Agent startup and configuration.
    - Invocation of CLI commands and their primary arguments.
    - Calls to core functions within modules (e.g., `lancedb_integration`, `graph_integration`, `vcf_utils`).
    - Database operations (connection attempts, table creation, additions, searches, updates, deletions) with key parameters.
    - Success or failure status of major operations.
    - Errors and exceptions, including stack traces where appropriate.
    - Key steps in AI-driven enrichment processes (if applicable).

- **Log Levels**: Standard log levels (DEBUG, INFO, WARNING, ERROR, CRITICAL) are used to control verbosity.
    - `INFO`: General operational messages.
    - `DEBUG`: Detailed information for troubleshooting.
    - `WARNING`: Potential issues or unexpected situations that don't halt operation.
    - `ERROR`/`CRITICAL`: Failures and critical problems.

- **PII and Sensitive Data Masking**:
    - **SQL Filters**: For LanceDB `filter_sql` parameters, a `mask_sensitive_sql` function (in `lancedb_integration.py`) is used to mask literal values in logs, especially for columns known to be sensitive (e.g., `patient_id`, `email`) or matching common PII patterns. Raw, unmasked filter strings are used for actual query execution but not logged directly if they are likely to contain sensitive data.
    - **Embedding Vectors**: Full embedding vectors are generally not logged; summaries (e.g., shape, length) may be logged instead.
    - **Update Dictionaries**: For `update_variant` operations, only the keys of the update dictionary are logged, not the full values, to avoid exposing sensitive data that might be part of an update.
    - **VCF Data**: Raw VCF record content is generally not logged unless in DEBUG mode for specific, limited troubleshooting.

- **Log Format**: Logs typically include a timestamp, log level, module name, and the message.

- **Log Output**: By default, logs are output to standard error/output. For production deployments, configure a logging handler to write to persistent files or a centralized logging system (e.g., ELK stack, Splunk).

Filesystem Auditing
-------------------

For locally stored data, particularly LanceDB and Kuzu database files, OS-level filesystem auditing is highly recommended.

- **Purpose**: To track all access (reads, writes, permission changes) to sensitive data directories. This is crucial for detecting unauthorized access or modifications and for forensic analysis.

- **Recommendations**:
    - **LanceDB Data Directory** (e.g., `./lancedb`):
        - Implement strict filesystem Access Control Lists (ACLs) allowing access only to the VCF Agent service user/process.
        - Enable filesystem audit logging for this directory. Tools include:
            - **Linux**: `auditd`. Configure rules to watch the LanceDB directory for read, write, attribute change events.
            - **macOS**: The built-in security auditing framework (`auditd` like, can be configured via `audit_control` and `audit_user` files).
            - **Windows**: Group Policy settings for Object Access Auditing.
    - **Kuzu Data Directory** (e.g., `./kuzu_db`):
        - Apply the same ACL and filesystem auditing principles as for the LanceDB directory.

- **Review Audit Logs**: Regularly review filesystem audit logs for suspicious activity. Integrate these logs into a Security Information and Event Management (SIEM) system if available.

Application-Level Auditing
--------------------------

Beyond basic logging:
- **Decision Records**: Key architectural and operational decisions are documented in `.context/decisions/`. These serve as an audit trail for design choices.
- **Task Management**: The Aegis framework's task files in `.context/tasks/` track the planning, execution, and completion of development and operational tasks, providing a form of process audit.
- **Version Control**: Git commit history serves as an audit trail for all code and configuration changes.

Review and Compliance
---------------------
- Periodically review logging and auditing configurations to ensure they meet current security and operational requirements.
- Align auditing practices with any relevant compliance standards (e.g., HIPAA, GDPR if applicable to the data being processed). 