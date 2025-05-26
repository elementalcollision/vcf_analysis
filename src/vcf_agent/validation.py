"""
VCF/BCF File I/O and Validation Logic
- File existence and format checks
- Index detection (bgzip, CSI/TBI)
- Integration with bcftools stats/validate
- Error reporting
"""

import os
from typing import Tuple, Optional
from .bcftools_integration import bcftools_stats
from .gatk_integration import gatk_validatevariants
from .config import CONFIG


def file_exists(filepath: str) -> bool:
    """
    Check if the file exists and is readable.

    Args:
        filepath (str): Path to the file.

    Returns:
        bool: True if the file exists and is readable, False otherwise.

    Example:
        >>> file_exists('sample_data/HG00098.vcf.gz')
        True
    """
    return os.path.isfile(filepath) and os.access(filepath, os.R_OK)


def is_vcf_or_bcf(filepath: str) -> bool:
    """
    Check if the file has a VCF, VCF.GZ, or BCF extension.

    Args:
        filepath (str): Path to the file.

    Returns:
        bool: True if the file has a recognized VCF/BCF extension, False otherwise.

    Example:
        >>> is_vcf_or_bcf('sample_data/HG00098.vcf.gz')
        True
    """
    ext = filepath.lower()
    return ext.endswith('.vcf') or ext.endswith('.vcf.gz') or ext.endswith('.bcf')


def has_index(filepath: str) -> bool:
    """
    Check for the presence of a CSI or TBI index for bgzipped VCF/BCF files.

    Args:
        filepath (str): Path to the file.

    Returns:
        bool: True if index is present or not required, False otherwise.

    Example:
        >>> has_index('sample_data/HG00098.vcf.gz')
        True
    """
    if filepath.endswith('.vcf.gz') or filepath.endswith('.bcf'):
        csi = filepath + '.csi'
        tbi = filepath + '.tbi'
        return os.path.exists(csi) or os.path.exists(tbi)
    return True  # Uncompressed VCF does not require index


def validate_with_bcftools_stats(filepath: str) -> Tuple[int, str, str]:
    """
    Run bcftools stats on the file.

    Args:
        filepath (str): Path to the VCF/BCF file.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Example:
        >>> validate_with_bcftools_stats('sample_data/HG00098.vcf.gz')
        (0, '...', '')
    """
    return bcftools_stats([filepath])


def validate_with_gatk(filepath: str, reference: Optional[str] = None) -> Tuple[int, str, str]:
    """
    Run GATK ValidateVariants on the file.

    Args:
        filepath (str): Path to the VCF/BCF file.
        reference (Optional[str]): Path to reference FASTA. If None, uses config.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Raises:
        ValueError: If reference FASTA is not provided.

    Example:
        >>> validate_with_gatk('sample_data/HG00098.vcf.gz', reference='ref.fasta')
        (0, '...', '')
    """
    if reference is None:
        reference = CONFIG.get("gatk_reference")
    if not reference:
        return 1, "", "Reference FASTA required for GATK ValidateVariants"
    return gatk_validatevariants(["-V", filepath, "-R", reference])


def validate_with_manufacturer_tool(filepath: str, tool_name: str) -> Tuple[int, str, str]:
    """
    Placeholder for manufacturer-specific compliance tool integration.

    Args:
        filepath (str): Path to the VCF/BCF file.
        tool_name (str): Name of the manufacturer tool.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Example:
        >>> validate_with_manufacturer_tool('sample_data/HG00098.vcf.gz', 'ont_vcf_validator')
        (1, '', "Manufacturer tool 'ont_vcf_validator' integration not implemented.")
    """
    return 1, "", f"Manufacturer tool '{tool_name}' integration not implemented."


def compliance_check(filepath: str, tool: Optional[str] = None) -> Tuple[int, str, str]:
    """
    Unified compliance checker for VCF/BCF files.
    Selects backend based on config or argument.

    Args:
        filepath (str): Path to the VCF/BCF file.
        tool (Optional[str]): Compliance tool to use ('bcftools', 'gatk', etc.).

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Example:
        >>> compliance_check('sample_data/HG00098.vcf.gz', tool='bcftools')
        (0, '...', '')
    """
    if tool is None:
        tool = CONFIG.get("compliance_tool", "bcftools")
    if tool == "bcftools":
        return validate_with_bcftools_stats(filepath)
    elif tool == "gatk":
        return validate_with_gatk(filepath)
    else:
        # Ensure tool_name is a string
        tool_name = tool if tool is not None else "unknown_tool"
        return validate_with_manufacturer_tool(filepath, tool_name)


def validate_vcf_file(filepath: str, tool: Optional[str] = None) -> Tuple[bool, Optional[str]]:
    """
    High-level validation: checks existence, format, index, and compliance (configurable backend).

    Args:
        filepath (str): Path to the VCF/BCF file.
        tool (Optional[str]): Compliance tool to use ('bcftools', 'gatk', etc.).

    Returns:
        Tuple[bool, Optional[str]]: (is_valid, error_message_if_any)

    Example:
        >>> validate_vcf_file('sample_data/HG00098.vcf.gz')
        (True, None)
    """
    if not file_exists(filepath):
        return False, f"File not found or not readable: {filepath}"
    if not is_vcf_or_bcf(filepath):
        return False, f"File does not have a recognized VCF/BCF extension: {filepath}"
    if not has_index(filepath):
        return False, f"Missing CSI or TBI index for bgzipped VCF/BCF file: {filepath}"
    
    rc, out, err = compliance_check(filepath, tool)
    
    # Determine the tool name used for logging/messaging
    used_tool_name = tool or CONFIG.get('compliance_tool', 'unknown_tool')

    if rc != 0:
        # If the command itself failed, report error from stderr
        error_message = err.strip() if err else "Unknown error (non-zero exit code with no stderr)"
        return False, f"Compliance check ({used_tool_name}) failed: {error_message}"

    # For bcftools, even with rc=0, check stdout for warnings
    # GATK and other tools might indicate errors differently (e.g., non-zero rc is primary)
    if used_tool_name == "bcftools":
        # Combine stdout and stderr for bcftools warnings, as some might go to stdout
        # and actual errors (if rc was non-zero) would have been in err.
        # Here, rc is 0, so we are looking for warnings in 'out' or 'err'.
        bcftools_output_warnings = []
        if out:
            for line in out.splitlines():
                if line.startswith("[W::") or line.startswith("[E::"):
                    bcftools_output_warnings.append(line)
        if err: # Also check stderr in case warnings ended up there with rc=0
            for line in err.splitlines():
                if line.startswith("[W::") or line.startswith("[E::"):
                    bcftools_output_warnings.append(line)
        
        if bcftools_output_warnings:
            # If warnings found in bcftools stats output, consider it a validation failure for strict checks
            return False, f"Compliance check (bcftools) passed with rc=0 but contained warnings: {'; '.join(bcftools_output_warnings)}"

    # If no non-zero rc, and for bcftools no specific warnings in stdout/stderr, then it's valid
    return True, None 