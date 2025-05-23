import os
import sys
import json
import pytest
import time

# Ensure src/ is in the path for import
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from vcf_agent.agent import vcf_comparison_tool

SAMPLE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../sample_data'))

# Example VCFs for different test cases
MINIMAL_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'edgecase_minimal.vcf.gz'))
MULTIALLELIC_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'edgecase_multiallelic.vcf.gz'))
EMPTY_ALT_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'edgecase_empty_alt.vcf.gz'))
INCONSISTENT_FORMAT_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'edgecase_inconsistent_format.vcf.gz'))
BAD_INFO_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'edgecase_bad_info.vcf.gz'))

# Placeholder for gold-standard benchmarking (replace with real files if available)
GOLD_STANDARD_VCF = MINIMAL_VCF
TEST_VCF = MULTIALLELIC_VCF

LARGE_VCF_1 = os.path.join(SAMPLE_DIR, '1KG.chr22.anno.vcf.gz')
LARGE_VCF_2 = os.path.join(SAMPLE_DIR, '1KG.chr22.anno.vcf.gz')  # Self-comparison for performance

REFERENCE_FASTA = os.path.join(SAMPLE_DIR, 'chr22.fa')

# Utility: skip if reference FASTA not available
pytestmark = pytest.mark.skipif(
    not os.path.isfile(REFERENCE_FASTA),
    reason='Reference FASTA required for normalization not found.'
)

def test_vcf_comparison_normalization():
    result = vcf_comparison_tool(MINIMAL_VCF, MINIMAL_VCF, REFERENCE_FASTA)
    data = json.loads(result)
    assert data['concordant_variant_count'] >= 0
    assert data['discordant_variant_count'] == 0
    assert isinstance(data['unique_to_file_1'], list)
    assert isinstance(data['unique_to_file_2'], list)
    assert isinstance(data['quality_metrics'], dict)
    assert 'per_sample_concordance' in data
    assert isinstance(data['per_sample_concordance'], dict)


def test_vcf_comparison_complex_variants():
    result = vcf_comparison_tool(MINIMAL_VCF, MULTIALLELIC_VCF, REFERENCE_FASTA)
    data = json.loads(result)
    assert 'concordant_variant_count' in data
    assert 'discordant_variant_count' in data
    assert isinstance(data['unique_to_file_1'], list)
    assert isinstance(data['unique_to_file_2'], list)
    assert isinstance(data['quality_metrics'], dict)
    assert 'per_sample_concordance' in data
    assert isinstance(data['per_sample_concordance'], dict)


def test_vcf_comparison_edge_cases():
    result = vcf_comparison_tool(MINIMAL_VCF, EMPTY_ALT_VCF, REFERENCE_FASTA)
    data = json.loads(result)
    assert 'concordant_variant_count' in data
    assert 'discordant_variant_count' in data
    assert isinstance(data['unique_to_file_1'], list)
    assert isinstance(data['unique_to_file_2'], list)
    assert isinstance(data['quality_metrics'], dict)
    assert 'per_sample_concordance' in data
    assert isinstance(data['per_sample_concordance'], dict)
    # Compare inconsistent format VCFs
    result2 = vcf_comparison_tool(MINIMAL_VCF, INCONSISTENT_FORMAT_VCF, REFERENCE_FASTA)
    data2 = json.loads(result2)
    assert 'concordant_variant_count' in data2
    assert 'discordant_variant_count' in data2
    assert 'per_sample_concordance' in data2
    assert isinstance(data2['per_sample_concordance'], dict)


def test_vcf_comparison_quality_metrics():
    result = vcf_comparison_tool(MINIMAL_VCF, BAD_INFO_VCF, REFERENCE_FASTA)
    data = json.loads(result)
    assert 'quality_metrics' in data
    assert isinstance(data['quality_metrics'], dict)
    assert 'per_sample_concordance' in data
    assert isinstance(data['per_sample_concordance'], dict)


def test_vcf_comparison_gold_standard():
    # Placeholder: compare gold-standard to test VCF (replace with real files for true benchmarking)
    result = vcf_comparison_tool(GOLD_STANDARD_VCF, TEST_VCF, REFERENCE_FASTA)
    data = json.loads(result)
    assert 'concordant_variant_count' in data
    assert 'discordant_variant_count' in data
    assert isinstance(data['unique_to_file_1'], list)
    assert isinstance(data['unique_to_file_2'], list)
    assert isinstance(data['quality_metrics'], dict)
    assert 'per_sample_concordance' in data
    assert isinstance(data['per_sample_concordance'], dict)


@pytest.mark.slow
def test_vcf_comparison_large_file_performance():
    if not os.path.isfile(LARGE_VCF_1):
        pytest.skip('Large VCF file not available')
    start = time.time()
    result = vcf_comparison_tool(LARGE_VCF_1, LARGE_VCF_2, REFERENCE_FASTA)
    elapsed = time.time() - start
    data = json.loads(result)
    assert 'concordant_variant_count' in data
    assert 'discordant_variant_count' in data
    assert isinstance(data['unique_to_file_1'], list)
    assert isinstance(data['unique_to_file_2'], list)
    assert isinstance(data['quality_metrics'], dict)
    assert 'per_sample_concordance' in data
    assert isinstance(data['per_sample_concordance'], dict)
    assert elapsed < 30  # seconds, adjust as needed for your environment 