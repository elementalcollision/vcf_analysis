"""
Python wrapper for bcftools view, query, and filter commands.
Uses subprocess for robust execution and output/error capture.
Extensible for additional bcftools commands.
"""

import subprocess
from typing import List, Optional, Tuple

def run_bcftools_command(command: List[str], input_data: Optional[bytes] = None) -> Tuple[int, str, str]:
    """
    Run a bcftools command and capture output and errors.

    Args:
        command (List[str]): List of command arguments (e.g., ["view", "-h", "file.vcf"])
        input_data (Optional[bytes]): Optional bytes to pass to stdin.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Raises:
        FileNotFoundError: If bcftools is not installed.

    Example:
        >>> run_bcftools_command(["view", "-h", "sample_data/HG00098.vcf.gz"])
        (0, '...', '')
    """
    full_cmd = ["bcftools"] + command
    try:
        result = subprocess.run(
            full_cmd,
            input=input_data,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )
        return result.returncode, result.stdout.decode(), result.stderr.decode()
    except FileNotFoundError as e:
        raise FileNotFoundError("bcftools is not installed or not in PATH") from e


def bcftools_view(args: List[str], input_data: Optional[bytes] = None) -> Tuple[int, str, str]:
    """
    Run bcftools view with the given arguments.

    Args:
        args (List[str]): Arguments for bcftools view.
        input_data (Optional[bytes]): Optional bytes to pass to stdin.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Example:
        >>> bcftools_view(["-h", "sample_data/HG00098.vcf.gz"])
        (0, '...', '')
    """
    return run_bcftools_command(["view"] + args, input_data)


def bcftools_query(args: List[str], input_data: Optional[bytes] = None) -> Tuple[int, str, str]:
    """
    Run bcftools query with the given arguments.

    Args:
        args (List[str]): Arguments for bcftools query.
        input_data (Optional[bytes]): Optional bytes to pass to stdin.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Example:
        >>> bcftools_query(["-f", "%CHROM\t%POS\n", "sample_data/HG00098.vcf.gz"])
        (0, '...', '')
    """
    return run_bcftools_command(["query"] + args, input_data)


def bcftools_filter(args: List[str], input_data: Optional[bytes] = None) -> Tuple[int, str, str]:
    """
    Run bcftools filter with the given arguments.

    Args:
        args (List[str]): Arguments for bcftools filter.
        input_data (Optional[bytes]): Optional bytes to pass to stdin.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Example:
        >>> bcftools_filter(["-e", "QUAL<20", "sample_data/HG00098.vcf.gz"])
        (0, '...', '')
    """
    return run_bcftools_command(["filter"] + args, input_data)


def bcftools_norm(args: List[str], input_data: Optional[bytes] = None) -> Tuple[int, str, str]:
    """
    Run bcftools norm with the given arguments.

    Args:
        args (List[str]): Arguments for bcftools norm.
        input_data (Optional[bytes]): Optional bytes to pass to stdin.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Example:
        >>> bcftools_norm(["-m", "-any", "sample_data/HG00098.vcf.gz"])
        (0, '...', '')
    """
    return run_bcftools_command(["norm"] + args, input_data)


def bcftools_stats(args: List[str], input_data: Optional[bytes] = None) -> Tuple[int, str, str]:
    """
    Run bcftools stats with the given arguments.

    Args:
        args (List[str]): Arguments for bcftools stats.
        input_data (Optional[bytes]): Optional bytes to pass to stdin.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Example:
        >>> bcftools_stats(["sample_data/HG00098.vcf.gz"])
        (0, '...', '')
    """
    return run_bcftools_command(["stats"] + args, input_data)


def bcftools_annotate(args: List[str], input_data: Optional[bytes] = None) -> Tuple[int, str, str]:
    """
    Run bcftools annotate with the given arguments.

    Args:
        args (List[str]): Arguments for bcftools annotate.
        input_data (Optional[bytes]): Optional bytes to pass to stdin.

    Returns:
        Tuple[int, str, str]: (returncode, stdout, stderr)

    Example:
        >>> bcftools_annotate(["-a", "annots.vcf", "sample_data/HG00098.vcf.gz"])
        (0, '...', '')
    """
    return run_bcftools_command(["annotate"] + args, input_data)

# Example extensibility: add more bcftools commands as needed 