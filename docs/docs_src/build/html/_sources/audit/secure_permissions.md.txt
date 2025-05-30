# Secure File and Directory Permissions Guide

## Purpose
This guide provides actionable steps to secure application data directories and files, minimizing the risk of unauthorized access or modification. It is intended for developers, sysadmins, and security teams.

---

## Summary Checklist
- **Principle of Least Privilege:** Grant only the minimum permissions needed.
- **Use Groups/Roles:** Assign permissions via groups, not individuals.
- **Restrict Sensitive Data:** Store sensitive files in dedicated, locked-down directories.
- **Regularly Review Permissions:** Audit and adjust as roles or needs change.
- **Monitor Access:** Enable logging/auditing for sensitive directories.
- **Avoid World-Readable/Writable:** Never use 777 (Linux/macOS) or 'Everyone:Full Control' (Windows).

---

## Linux/macOS
### Recommended Permissions
- **Directories:** `chmod 700 /path/to/dir` (owner only)
- **Files:** `chmod 600 /path/to/file` (owner only)
- **Set Owner:** `chown user:group /path/to/dir_or_file`
- **Advanced (ACLs):** `setfacl -m u:username:rwx /path/to/dir`

### Example: Secure a LanceDB Data Directory
```sh
# Set owner to 'lancedbuser' and restrict access
groupadd lancedb
useradd -g lancedb lancedbuser
chown -R lancedbuser:lancedb /var/lib/lancedb
chmod 700 /var/lib/lancedb
```

### Check Permissions
```sh
ls -ld /var/lib/lancedb
ls -l /var/lib/lancedb
```

---

## Windows
### Recommended Permissions
- **NTFS Permissions:**
  - Remove 'Everyone' and 'Users' groups from sensitive directories.
  - Grant 'Full Control' only to the service account or Administrators.
- **Set via GUI:**
  - Right-click folder > Properties > Security > Edit
- **Set via PowerShell:**
```powershell
$folder = "C:\\path\\to\\lancedb"
$acl = Get-Acl $folder
$rule = New-Object System.Security.AccessControl.FileSystemAccessRule("LanceDBUser","FullControl","Allow")
$acl.SetAccessRule($rule)
Set-Acl $folder $acl
```

### Check Permissions
- Use `icacls C:\path\to\lancedb`

---

## Review and Monitoring
- **Audit Regularly:** Schedule permission reviews (quarterly or after role changes).
- **Enable Auditing:** See [Filesystem Auditing Guide](filesystem_auditing.md).
- **Log Access:** Monitor for suspicious or failed access attempts.

---

## References
- [OWASP File Permission Guide](https://owasp.org/www-project-web-security-testing-guide/latest/4-Web_Application_Security_Testing/02-Configuration_and_Deployment_Management_Testing/09-Test_File_Permission)
- [Microsoft File and Folder Permissions](https://docs.microsoft.com/en-us/windows/security/threat-protection/security-policy-settings/access-control-permissions)
- [Linux chmod, chown, setfacl](https://man7.org/linux/man-pages/man1/chmod.1.html)

> **Note:** Always tailor permissions to your organization's policy and application requirements. 