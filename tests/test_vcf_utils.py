import os
import sys
import pytest

# Ensure src/ is in the path for import
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from vcf_agent import vcf_utils

SAMPLE_DIR = os.path.join(os.path.dirname(__file__), '..', 'sample_data')

MINIMAL_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'edgecase_minimal.vcf'))
MULTIALLELIC_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'edgecase_multiallelic.vcf'))


def test_extract_variant_summary_minimal():
    summary = vcf_utils.extract_variant_summary(MINIMAL_VCF)
    assert isinstance(summary, dict)
    assert 'variant_count' in summary
    assert 'variant_types' in summary
    assert 'samples' in summary
    assert 'sample_variant_counts' in summary
    assert summary['variant_count'] > 0
    assert isinstance(summary['samples'], list)
    assert isinstance(summary['sample_variant_counts'], dict)


def test_extract_variant_summary_multiallelic():
    summary = vcf_utils.extract_variant_summary(MULTIALLELIC_VCF)
    assert isinstance(summary, dict)
    assert summary['variant_count'] > 0
    assert isinstance(summary['variant_types'], dict)


def test_get_sample_names():
    samples = vcf_utils.get_sample_names(MINIMAL_VCF)
    assert isinstance(samples, list)
    assert len(samples) > 0


def test_extract_comparable_features():
    features = vcf_utils.extract_comparable_features(MINIMAL_VCF)
    assert isinstance(features, dict)
    assert 'variant_count' in features
    assert 'variant_types' in features


def test_extract_comparable_features_multiallelic():
    features = vcf_utils.extract_comparable_features(MULTIALLELIC_VCF)
    assert isinstance(features, dict)
    assert 'variant_count' in features
    assert 'variant_types' in features
    assert features['variant_count'] > 0
    assert isinstance(features['variant_types'], dict) 