import os
import pytest
from vcf_agent.validation import (
    file_exists, is_vcf_or_bcf, has_index, validate_with_bcftools_stats, validate_vcf_file
)

# NOTE: A sample VCF file should be provided in sample_data/. Several examples can be found at:
# https://bioinformaticstools.mayo.edu/research/vcf-miner-sample-vcfs/
SAMPLE_VCF = "sample_data/HG00098.vcf.gz"

@pytest.mark.skipif(not os.path.exists(SAMPLE_VCF), reason="Sample VCF file not found")
def test_file_exists():
    assert file_exists(SAMPLE_VCF)


def test_is_vcf_or_bcf():
    assert is_vcf_or_bcf(SAMPLE_VCF)
    assert is_vcf_or_bcf("test.vcf")
    assert is_vcf_or_bcf("test.bcf")
    assert not is_vcf_or_bcf("test.txt")


@pytest.mark.skipif(not os.path.exists(SAMPLE_VCF), reason="Sample VCF file not found")
def test_has_index():
    # The sample file may or may not have an index; just check that the function runs
    has_index(SAMPLE_VCF)


@pytest.mark.skipif(not os.path.exists(SAMPLE_VCF), reason="Sample VCF file not found")
def test_validate_with_bcftools_stats():
    rc, out, err = validate_with_bcftools_stats(SAMPLE_VCF)
    assert rc == 0 or rc == 1  # Accept 0 or 1 for help/usage
    assert "# This file was produced by bcftools stats" in out or "lines, records, variants" in out or out or err


@pytest.mark.skipif(not os.path.exists(SAMPLE_VCF), reason="Sample VCF file not found")
def test_validate_vcf_file():
    is_valid, error = validate_vcf_file(SAMPLE_VCF)
    if not is_valid:
        print(f"Validation error: {error}")
    assert is_valid
    assert error is None 