"""
vcf_utils.py
Module for extracting and formatting VCF data for LLM-powered analysis.
Uses cyvcf2 for efficient VCF parsing.
"""

from cyvcf2 import VCF
from typing import Dict, Any, List
import os


def extract_variant_summary(vcf_path: str) -> Dict[str, Any]:
    """
    Extracts summary statistics from a VCF file.
    Returns a dictionary with variant counts, types, and sample stats.
    """
    if not os.path.isfile(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    vcf = VCF(vcf_path)
    variant_count = 0
    variant_types = {}
    samples = vcf.samples
    sample_counts = {sample: 0 for sample in samples}

    for record in vcf:
        variant_count += 1
        vtype = record.var_type  # e.g., 'snp', 'indel', etc.
        variant_types[vtype] = variant_types.get(vtype, 0) + 1
        # Count variants per sample (GT field)
        for i, gt in enumerate(record.genotypes):
            # gt is a list like [0, 1, True] (last is phased)
            if any(allele > 0 for allele in gt[:2]):
                sample_counts[samples[i]] += 1

    return {
        "variant_count": variant_count,
        "variant_types": variant_types,
        "samples": samples,
        "sample_variant_counts": sample_counts,
    }


def extract_comparable_features(vcf_path: str) -> Dict[str, Any]:
    """
    Extracts features from a VCF file for comparison tasks.
    Returns a dictionary with key features for side-by-side comparison.
    """
    summary = extract_variant_summary(vcf_path)
    # Add more features as needed for comparison
    return summary


def get_sample_names(vcf_path: str) -> List[str]:
    """
    Returns the list of sample names in the VCF file.
    """
    vcf = VCF(vcf_path)
    return vcf.samples 