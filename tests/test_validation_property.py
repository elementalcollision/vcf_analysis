import os
import tempfile
import pytest
from hypothesis import given, strategies as st, assume
from vcf_agent.validation import is_vcf_or_bcf, file_exists, validate_vcf_file
from vcf_agent.bcftools_integration import bcftools_view

# Safe ASCII for file suffixes and args
safe_ascii = st.text(alphabet=list("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.-_"), min_size=1, max_size=5)
safe_arg = st.text(alphabet=list("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_=."), min_size=1, max_size=10)

# Property-based test for is_vcf_or_bcf with random extensions
@given(st.text(min_size=1) | st.just(""))
def test_is_vcf_or_bcf_random_ext(ext):
    filename = f"foo.{ext}"
    result = is_vcf_or_bcf(filename)
    # Should only be True for .vcf, .vcf.gz, .bcf (case-insensitive)
    valid = filename.lower().endswith('.vcf') or filename.lower().endswith('.vcf.gz') or filename.lower().endswith('.bcf')
    assert result == valid

# Property-based test for file_exists with random file names
@given(safe_ascii)
def test_file_exists_random_name(name):
    # Should be False for random names that don't exist
    assert not file_exists(f"sample_data/{name}")

# Property-based test for file_exists with temp files
@given(safe_ascii)
def test_file_exists_tempfile(suffix):
    with tempfile.NamedTemporaryFile(suffix=suffix) as tmp:
        assert file_exists(tmp.name)

# Property-based test for bcftools_view with random arguments
@given(st.lists(safe_arg, min_size=1, max_size=3))
def test_bcftools_view_random_args(args):
    # Avoid null bytes and empty args
    assume(all(arg and "\x00" not in arg for arg in args))
    rc, out, err = bcftools_view(args)
    # Should not crash; rc may be nonzero for invalid args
    assert isinstance(rc, int)
    assert isinstance(out, str)
    assert isinstance(err, str)

# Property-based test for validate_vcf_file with random file names/extensions
@given(safe_ascii)
def test_validate_vcf_file_random_name(name):
    is_valid, error = validate_vcf_file(f"sample_data/{name}")
    assert not is_valid or error is None  # Should be invalid unless a real file 