import os
import tempfile
import pytest
from vcf_agent.validation import (
    file_exists, is_vcf_or_bcf, has_index, validate_vcf_file
)
from vcf_agent.bcftools_integration import (
    bcftools_view, bcftools_query, bcftools_filter, bcftools_stats
)

# Nonexistent file
NONEXISTENT = "sample_data/does_not_exist.vcf.gz"

# Unsupported extension
UNSUPPORTED = "sample_data/unsupported.txt"

def test_file_exists_nonexistent():
    assert not file_exists(NONEXISTENT)

def test_is_vcf_or_bcf_unsupported():
    assert not is_vcf_or_bcf(UNSUPPORTED)
    assert not is_vcf_or_bcf("foo.bam")

def test_validate_vcf_file_nonexistent():
    is_valid, error = validate_vcf_file(NONEXISTENT)
    assert not is_valid
    assert error is not None
    assert ("not found" in error or "not readable" in error)

def test_validate_vcf_file_unsupported():
    # Create a temp .txt file
    with tempfile.NamedTemporaryFile(suffix=".txt") as tmp:
        is_valid, error = validate_vcf_file(tmp.name)
        assert not is_valid
        assert error is not None
        assert "recognized VCF/BCF extension" in error

def test_validate_vcf_file_missing_index():
    # Create a temp .vcf.gz file with no index
    with tempfile.NamedTemporaryFile(suffix=".vcf.gz") as tmp:
        tmp.write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        tmp.flush()
        is_valid, error = validate_vcf_file(tmp.name)
        assert not is_valid
        assert error is not None
        assert "Missing CSI or TBI index" in error

def test_bcftools_view_invalid_args():
    # Should fail with invalid args
    rc, out, err = bcftools_view(["--notarealoption"]) 
    assert rc != 0
    assert err

def test_bcftools_stats_nonexistent_file():
    rc, out, err = bcftools_stats([NONEXISTENT])
    assert rc != 0
    assert ("No such file" in err or "cannot open" in err or err) 