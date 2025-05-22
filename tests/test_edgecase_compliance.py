import os
import pytest
from vcf_agent.validation import validate_vcf_file

EDGE_CASES = [
    ("edgecase_missing_header.vcf", "missing header"),
    ("edgecase_bad_info.vcf", "bad info field"),
    ("edgecase_nonstandard_chrom.vcf", "nonstandard chrom"),
    ("edgecase_symbolic_allele.vcf", "symbolic allele"),
    ("edgecase_multiallelic.vcf", "multiallelic site"),
    ("edgecase_missing_format.vcf", "missing format"),
    ("edgecase_invalid_qual.vcf", "invalid QUAL value"),
    ("edgecase_missing_filter.vcf", "missing FILTER field"),
    ("edgecase_info_trailing_semicolon.vcf", "INFO trailing semicolon"),
    ("edgecase_duplicate_info.vcf", "duplicate INFO key"),
    ("edgecase_inconsistent_format.vcf", "inconsistent FORMAT/sample"),
    ("edgecase_empty_alt.vcf", "empty ALT field"),
]

@pytest.mark.parametrize("filename,desc", EDGE_CASES)
def test_edgecase_bcftools(filename, desc):
    path = os.path.join("sample_data", filename)
    if not os.path.exists(path):
        pytest.skip(f"{filename} not found")
    is_valid, error = validate_vcf_file(path, tool="bcftools")
    assert not is_valid, f"{filename} should be invalid ({desc})"
    assert error, f"{filename} should return an error message"

@pytest.mark.parametrize("filename,desc", EDGE_CASES)
def test_edgecase_gatk(filename, desc):
    path = os.path.join("sample_data", filename)
    if not os.path.exists(path):
        pytest.skip(f"{filename} not found")
    is_valid, error = validate_vcf_file(path, tool="gatk")
    assert not is_valid, f"{filename} should be invalid ({desc})"
    assert error, f"{filename} should return an error message" 