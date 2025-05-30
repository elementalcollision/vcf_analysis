# Filesystem Auditing Guide for LanceDB Data Directory

## Purpose
Filesystem auditing helps detect unauthorized or suspicious access to the LanceDB data directory, supporting compliance, incident response, and forensic analysis. This guide provides platform-specific instructions and best practices.

---

## macOS
### Tools
- **fs_usage**: Real-time file access monitoring.
- **audit**: Built-in auditing framework (OpenBSM).

### Example: Monitor LanceDB Directory with `fs_usage`
```sh
sudo fs_usage -w -f filesys | grep /path/to/lancedb
```

### Enabling File Auditing with OpenBSM
1. Edit `/etc/security/audit_control` to include the LanceDB directory:
   ```
   dir:/path/to/lancedb
   ``
2. Restart auditd:
   ```sh
   sudo launchctl stop com.apple.auditd
   sudo launchctl start com.apple.auditd
   ```
3. Review logs in `/var/audit/`.

---

## Linux
### Tools
- **auditd**: Linux Audit Daemon

### Example: Audit LanceDB Directory
1. Add a rule (replace `/path/to/lancedb`):
   ```sh
   sudo auditctl -w /path/to/lancedb -p rwxa -k lancedb_audit
   ```
2. Review logs:
   ```sh
   sudo ausearch -k lancedb_audit
   sudo aureport -f
   ```
3. Make rules persistent by adding to `/etc/audit/rules.d/lancedb.rules`:
   ```
   -w /path/to/lancedb -p rwxa -k lancedb_audit
   ```

---

## Windows
### Tools
- **Group Policy**: Enable Object Access auditing
- **PowerShell**: Configure auditing on folders

### Example: Enable Auditing via PowerShell
```powershell
# Enable auditing for a folder
$folder = "C:\\path\\to\\lancedb"
$acl = Get-Acl $folder
$audit = New-Object System.Security.AccessControl.FileSystemAuditRule("Everyone","FullControl","Success,Failure")
$acl.AddAuditRule($audit)
Set-Acl $folder $acl
```

### Review Events
- Use Event Viewer: Windows Logs > Security (look for Event ID 4663)

---

## Best Practices
- **Restrict permissions**: Only allow necessary users access.
- **Monitor logs**: Set up alerts for suspicious activity.
- **Log retention**: Store logs securely and retain per compliance policy.
- **Review regularly**: Periodically review audit logs for anomalies.

## Example: Identifying Suspicious Access
- Multiple failed access attempts
- Access from unexpected users or times
- Unusual file modifications or deletions

## References
- [Apple OpenBSM Audit](https://www.apple.com/business/docs/IT_Security_Guide.pdf)
- [Linux Audit Documentation](https://linux.die.net/man/8/auditd)
- [Microsoft File and Folder Auditing](https://docs.microsoft.com/en-us/windows/security/threat-protection/auditing/basic-audit-object-access)

---

> **Note:** This guide provides general recommendations. Tailor auditing to your organization's compliance, privacy, and operational requirements. 