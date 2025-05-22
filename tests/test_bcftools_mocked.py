import pytest
from unittest.mock import patch, MagicMock
from vcf_agent.bcftools_integration import (
    bcftools_view, bcftools_query, bcftools_filter,
    bcftools_norm, bcftools_stats, bcftools_annotate
)

SUCCESS = MagicMock(return_value=MagicMock(returncode=0, stdout=b"success output", stderr=b""))
FAILURE = MagicMock(return_value=MagicMock(returncode=1, stdout=b"", stderr=b"error output"))

@patch("subprocess.run", SUCCESS)
def test_bcftools_view_success():
    rc, out, err = bcftools_view(["--help"])
    assert rc == 0
    assert out == "success output"
    assert err == ""

@patch("subprocess.run", FAILURE)
def test_bcftools_view_failure():
    rc, out, err = bcftools_view(["--badarg"])
    assert rc == 1
    assert out == ""
    assert err == "error output"

@patch("subprocess.run", SUCCESS)
def test_bcftools_query_success():
    rc, out, err = bcftools_query(["--help"])
    assert rc == 0
    assert out == "success output"
    assert err == ""

@patch("subprocess.run", FAILURE)
def test_bcftools_query_failure():
    rc, out, err = bcftools_query(["--badarg"])
    assert rc == 1
    assert out == ""
    assert err == "error output"

@patch("subprocess.run", SUCCESS)
def test_bcftools_filter_success():
    rc, out, err = bcftools_filter(["--help"])
    assert rc == 0
    assert out == "success output"
    assert err == ""

@patch("subprocess.run", FAILURE)
def test_bcftools_filter_failure():
    rc, out, err = bcftools_filter(["--badarg"])
    assert rc == 1
    assert out == ""
    assert err == "error output"

@patch("subprocess.run", SUCCESS)
def test_bcftools_norm_success():
    rc, out, err = bcftools_norm(["--help"])
    assert rc == 0
    assert out == "success output"
    assert err == ""

@patch("subprocess.run", FAILURE)
def test_bcftools_norm_failure():
    rc, out, err = bcftools_norm(["--badarg"])
    assert rc == 1
    assert out == ""
    assert err == "error output"

@patch("subprocess.run", SUCCESS)
def test_bcftools_stats_success():
    rc, out, err = bcftools_stats(["--help"])
    assert rc == 0
    assert out == "success output"
    assert err == ""

@patch("subprocess.run", FAILURE)
def test_bcftools_stats_failure():
    rc, out, err = bcftools_stats(["--badarg"])
    assert rc == 1
    assert out == ""
    assert err == "error output"

@patch("subprocess.run", SUCCESS)
def test_bcftools_annotate_success():
    rc, out, err = bcftools_annotate(["--help"])
    assert rc == 0
    assert out == "success output"
    assert err == ""

@patch("subprocess.run", FAILURE)
def test_bcftools_annotate_failure():
    rc, out, err = bcftools_annotate(["--badarg"])
    assert rc == 1
    assert out == ""
    assert err == "error output" 