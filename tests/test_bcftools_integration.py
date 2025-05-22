import pytest
from vcf_agent.bcftools_integration import (
    bcftools_view, bcftools_query, bcftools_filter,
    bcftools_norm, bcftools_stats, bcftools_annotate
)
import shutil

bcftools_exists = shutil.which("bcftools") is not None

pytestmark = pytest.mark.skipif(not bcftools_exists, reason="bcftools is not installed")

def test_bcftools_view_help():
    returncode, stdout, stderr = bcftools_view(["--help"])
    assert (returncode == 0 or returncode == 1)
    assert "Usage:   bcftools view" in stdout or "Usage:   bcftools view" in stderr

def test_bcftools_query_help():
    returncode, stdout, stderr = bcftools_query(["-h"])
    assert (returncode == 0 or returncode == 1)
    assert "Usage:   bcftools query" in stdout or "Usage:   bcftools query" in stderr

def test_bcftools_filter_help():
    returncode, stdout, stderr = bcftools_filter(["-h"])
    assert (returncode == 0 or returncode == 1)
    assert "Usage:   bcftools filter" in stdout or "Usage:   bcftools filter" in stderr

def test_bcftools_norm_help():
    returncode, stdout, stderr = bcftools_norm(["--help"])
    assert (returncode == 0 or returncode == 1)
    assert "Usage:   bcftools norm" in stdout or "Usage:   bcftools norm" in stderr

def test_bcftools_stats_help():
    returncode, stdout, stderr = bcftools_stats(["--help"])
    assert (returncode == 0 or returncode == 1)
    assert "Usage:   bcftools stats" in stdout or "Usage:   bcftools stats" in stderr

def test_bcftools_annotate_help():
    returncode, stdout, stderr = bcftools_annotate(["--help"])
    assert (returncode == 0 or returncode == 1)
    assert "Usage:   bcftools annotate" in stdout or "Usage:   bcftools annotate" in stderr 