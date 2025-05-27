import os
import sys
import json
import pytest
import time

# Ensure src/ is in the path for import
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from vcf_agent.agent import vcf_comparison_tool

# Define the sample directory relative to this test file
# Assumes sample VCFs are in a 'sample_data' directory parallel to 'tests'
SAMPLE_DIR = os.path.join(os.path.dirname(__file__), '..', 'sample_test_data')

# Example VCFs for different test cases
MINIMAL_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'minimal.vcf.gz'))
MULTIALLELIC_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'multiallelic.vcf.gz'))
EMPTY_ALT_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'empty_alt.vcf.gz'))
INCONSISTENT_FORMAT_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'inconsistent_format.vcf.gz'))
BAD_INFO_VCF = os.path.abspath(os.path.join(SAMPLE_DIR, 'bad_info.vcf.gz'))

# Placeholder for gold-standard benchmarking (replace with real files if available)
GOLD_STANDARD_VCF = MINIMAL_VCF
TEST_VCF = MULTIALLELIC_VCF

LARGE_VCF_FILENAME = 'chr22.1kg.phase3.v5a.testready.vcf.gz'
LARGE_VCF_PATH = os.path.join(SAMPLE_DIR, LARGE_VCF_FILENAME)

REFERENCE_FASTA = os.path.join(SAMPLE_DIR, 'chr22.fa')
LARGE_VCF_REFERENCE_FASTA = os.path.join(SAMPLE_DIR, '22.fa')

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
    if not os.path.isfile(LARGE_VCF_PATH):
        pytest.skip(f'Large VCF file ({LARGE_VCF_FILENAME}) not found in {SAMPLE_DIR}. See sample_test_data/README.md for download and preparation instructions.')
    if not os.path.isfile(LARGE_VCF_REFERENCE_FASTA):
        pytest.skip(f'Reference FASTA for large VCF ({os.path.basename(LARGE_VCF_REFERENCE_FASTA)}) not found in {SAMPLE_DIR}. This should be created by test setup or manually as per extended instructions.')
    start = time.time()
    result = vcf_comparison_tool(LARGE_VCF_PATH, LARGE_VCF_PATH, LARGE_VCF_REFERENCE_FASTA)
    elapsed = time.time() - start
    data = json.loads(result)
    assert 'concordant_variant_count' in data
    assert 'discordant_variant_count' in data
    assert isinstance(data['unique_to_file_1_count'], int)
    assert isinstance(data['unique_to_file_2_count'], int)
    assert 'per_sample_concordance' in data
    assert isinstance(data['quality_metrics'], dict)
    assert elapsed < 300  # seconds, initial baseline for a ~130MB VCF file 