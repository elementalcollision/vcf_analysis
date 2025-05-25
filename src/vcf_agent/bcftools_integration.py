"""
Python wrapper for bcftools view, query, and filter commands.
Uses subprocess for robust execution and output/error capture.
Extensible for additional bcftools commands.
"""

import subprocess
from typing import List, Optional, Tuple
import tempfile
import shutil
import time # For timing
from . import metrics # Import metrics module

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
    bcftools_subcommand = command[0] if command else "unknown_subcommand"
    start_time = time.time()
    status = "success" # Assume success initially
    return_code = -1 # Default to -1 if subprocess.run fails before setting result.returncode

    try:
        result = subprocess.run(
            full_cmd,
            input=input_data,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )
        return_code = result.returncode
        if return_code != 0:
            status = "error"
        return result.returncode, result.stdout.decode(), result.stderr.decode()
    except FileNotFoundError as e:
        status = "error"
        # Log specific error for bcftools not found, but metrics will capture generic command error
        metrics.log.error("bcftools_not_found", command_used=full_cmd, error=str(e))
        raise FileNotFoundError("bcftools is not installed or not in PATH") from e
    except Exception as e: # Catch any other unexpected errors during subprocess.run
        status = "error"
        metrics.log.error("bcftools_run_unexpected_error", command_used=full_cmd, error=str(e))
        # Re-raise to allow higher-level handlers to deal with it, 
        # but metrics will be recorded in finally.
        raise
    finally:
        duration = time.time() - start_time
        metrics.VCF_AGENT_BCFTOOLS_DURATION_SECONDS.labels(
            bcftools_subcommand=bcftools_subcommand, status=status
        ).observe(duration)
        metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels(
            bcftools_subcommand=bcftools_subcommand, status=status
        ).inc()
        
        if status == "error":
            # We use return_code !=0 as the primary indicator for bcftools errors for the metric
            # FileNotFoundError and other exceptions during subprocess.run also result in status="error"
            metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL.labels(
                bcftools_subcommand=bcftools_subcommand
            ).inc()


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

    Example::

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

def bcftools_isec(file1: str, file2: str) -> dict:
    """
    Run bcftools isec to compare two VCF files.
    Returns paths to unique and common VCFs in a temp directory.
    """
    temp_dir = tempfile.mkdtemp(prefix="bcftools_isec_")
    try:
        rc, out, err = run_bcftools_command(["isec", file1, file2, "-p", temp_dir])
        if rc != 0:
            raise RuntimeError(f"bcftools isec failed: {err}")
        return {
            "temp_dir": temp_dir,
            "unique_to_file1": f"{temp_dir}/0000.vcf",
            "unique_to_file2": f"{temp_dir}/0001.vcf",
            "concordant": f"{temp_dir}/0002.vcf"
        }
    except Exception as e:
        shutil.rmtree(temp_dir)
        raise e

def count_variants_in_vcf(vcf_path: str) -> int:
    """
    Count the number of variant records in a VCF file (ignores header lines).
    """
    count = 0
    with open(vcf_path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                count += 1
    return count

def parse_vcf_variants(vcf_path):
    """
    Parse a VCF file and return a list of variant dicts (CHROM, POS, REF, ALT).
    Skips header lines.
    """
    variants = []
    try:
        with open(vcf_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue  # skip malformed lines
                chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
                variants.append({
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt
                })
    except Exception as e:
        raise RuntimeError(f"Failed to parse VCF {vcf_path}: {e}")
    return variants

def vcf_compare(file1: str, file2: str) -> dict:
    """
    Compare two VCF files using bcftools isec and return concordant/discordant/unique counts and lists.
    Returns a JSON-ready dict. Always returns all required fields, even on error.
    """
    try:
        isec = bcftools_isec(file1, file2)
        unique1 = parse_vcf_variants(isec["unique_to_file1"])
        unique2 = parse_vcf_variants(isec["unique_to_file2"])
        concordant = parse_vcf_variants(isec["concordant"])
        result = {
            "concordant_variant_count": len(concordant),
            "discordant_variant_count": len(unique1) + len(unique2),
            "unique_to_file_1": unique1,
            "unique_to_file_2": unique2,
            "quality_metrics": {}
        }
        return result
    except Exception as e:
        return {
            "concordant_variant_count": 0,
            "discordant_variant_count": 0,
            "unique_to_file_1": [],
            "unique_to_file_2": [],
            "quality_metrics": {},
            "error": str(e)
        }
    finally:
        # Clean up temp directory if it exists
        try:
            if 'isec' in locals() and "temp_dir" in isec:
                shutil.rmtree(isec["temp_dir"], ignore_errors=True)
        except Exception:
            pass
    # Fallback return to satisfy linter
    return {
        "concordant_variant_count": 0,
        "discordant_variant_count": 0,
        "unique_to_file_1": [],
        "unique_to_file_2": [],
        "quality_metrics": {},
        "error": "Unknown error in vcf_compare"
    }

# Example extensibility: add more bcftools commands as needed 