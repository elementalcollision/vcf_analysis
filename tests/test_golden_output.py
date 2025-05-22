import os
import pytest
from vcf_agent.bcftools_integration import bcftools_stats

def normalize(text):
    return '\n'.join(line.strip() for line in text.strip().splitlines() if line.strip())

SAMPLE_VCF = "sample_data/HG00098.vcf.gz"
GOLDEN_STATS = "golden/HG00098.stats.txt"

@pytest.mark.skipif(not os.path.exists(SAMPLE_VCF), reason="Sample VCF file not found")
@pytest.mark.skipif(not os.path.exists(GOLDEN_STATS), reason="Golden stats file not found")
def test_bcftools_stats_golden():
    rc, out, err = bcftools_stats([SAMPLE_VCF])
    # Use output or error if output is empty
    actual = out if out else err
    with open(GOLDEN_STATS) as f:
        golden = f.read()
    assert normalize(actual) == normalize(golden), "bcftools_stats output does not match golden file. If this is intentional, update the golden file." 