"""
Golden File Tests

Tests to validate that current VCF processing outputs match expected golden files.
These tests ensure that changes to the codebase don't introduce regressions.
"""

import pytest
import json
from pathlib import Path

from vcf_agent.bcftools_integration import (
    bcftools_view,
    bcftools_stats,
    bcftools_query
)
from vcf_agent.validation import validate_vcf_file
from tests.utils.golden_file_utils import (
    assert_golden_match,
    GoldenFileManager,
    GoldenFileError
)


class TestBcftoolsGoldenFiles:
    """Test bcftools operations against golden files."""
    
    def test_bcftools_view_header_small_valid(self):
        """Test bcftools view header output against golden file."""
        rc, stdout, stderr = bcftools_view(["-h", "sample_test_data/small_valid.vcf"])
        assert rc == 0, f"bcftools view failed: {stderr}"
        
        assert_golden_match(
            actual=stdout,
            golden_path="bcftools/view/view_header_small_valid.txt",
            content_type="text"
        )
    
    def test_bcftools_view_header_sample1(self):
        """Test bcftools view header output for sample1."""
        rc, stdout, stderr = bcftools_view(["-h", "sample_test_data/sample1.vcf"])
        assert rc == 0, f"bcftools view failed: {stderr}"
        
        assert_golden_match(
            actual=stdout,
            golden_path="bcftools/view/view_header_sample1.txt",
            content_type="text"
        )
    
    def test_bcftools_stats_small_valid(self):
        """Test bcftools stats output against golden file."""
        rc, stdout, stderr = bcftools_stats(["sample_test_data/small_valid.vcf"])
        assert rc == 0, f"bcftools stats failed: {stderr}"
        
        assert_golden_match(
            actual=stdout,
            golden_path="bcftools/stats/stats_small_valid.txt",
            content_type="text"
        )
    
    def test_bcftools_stats_sample1(self):
        """Test bcftools stats output for sample1."""
        rc, stdout, stderr = bcftools_stats(["sample_test_data/sample1.vcf"])
        assert rc == 0, f"bcftools stats failed: {stderr}"
        
        assert_golden_match(
            actual=stdout,
            golden_path="bcftools/stats/stats_sample1.txt",
            content_type="text"
        )
    
    def test_bcftools_query_basic_small_valid(self):
        """Test bcftools query basic format against golden file."""
        rc, stdout, stderr = bcftools_query(["-f", "%CHROM\t%POS\t%REF\t%ALT\n", "sample_test_data/small_valid.vcf"])
        assert rc == 0, f"bcftools query failed: {stderr}"
        
        assert_golden_match(
            actual=stdout,
            golden_path="bcftools/query/query_basic_small_valid.txt",
            content_type="text"
        )
    
    def test_bcftools_query_basic_sample1(self):
        """Test bcftools query basic format for sample1."""
        rc, stdout, stderr = bcftools_query(["-f", "%CHROM\t%POS\t%REF\t%ALT\n", "sample_test_data/sample1.vcf"])
        assert rc == 0, f"bcftools query failed: {stderr}"
        
        assert_golden_match(
            actual=stdout,
            golden_path="bcftools/query/query_basic_sample1.txt",
            content_type="text"
        )


class TestValidationGoldenFiles:
    """Test VCF validation operations against golden files."""
    
    def test_validation_small_valid(self):
        """Test VCF validation output against golden file."""
        validation_result = validate_vcf_file("sample_test_data/small_valid.vcf")
        
        # Convert to JSON for comparison
        result_json = json.dumps(validation_result, indent=2, sort_keys=True)
        
        assert_golden_match(
            actual=result_json,
            golden_path="validation/valid_files/validation_small_valid.json",
            content_type="json"
        )
    
    def test_validation_sample1(self):
        """Test VCF validation output for sample1."""
        validation_result = validate_vcf_file("sample_test_data/sample1.vcf")
        
        # Convert to JSON for comparison
        result_json = json.dumps(validation_result, indent=2, sort_keys=True)
        
        assert_golden_match(
            actual=result_json,
            golden_path="validation/valid_files/validation_sample1.json",
            content_type="json"
        )


class TestGoldenFileFramework:
    """Test the golden file framework itself."""
    
    def test_golden_file_manager_initialization(self):
        """Test that GoldenFileManager initializes correctly."""
        manager = GoldenFileManager()
        assert manager.golden_dir.exists()
        assert manager.golden_dir.name == "golden"
    
    def test_golden_file_exists(self):
        """Test golden file existence checking."""
        manager = GoldenFileManager()
        
        # Test existing file
        assert manager.golden_file_exists("bcftools/view/view_header_small_valid.txt")
        
        # Test non-existing file
        assert not manager.golden_file_exists("nonexistent/file.txt")
    
    def test_load_golden_file(self):
        """Test loading golden files."""
        manager = GoldenFileManager()
        
        content = manager.load_golden_file("bcftools/view/view_header_small_valid.txt")
        assert content is not None
        assert len(content) > 0
        assert "##fileformat=VCFv4.2" in content
    
    def test_load_nonexistent_golden_file(self):
        """Test loading non-existent golden file raises error."""
        manager = GoldenFileManager()
        
        with pytest.raises(GoldenFileError, match="Golden file not found"):
            manager.load_golden_file("nonexistent/file.txt")
    
    def test_save_and_load_golden_file(self, tmp_path):
        """Test saving and loading golden files."""
        # Use temporary directory for this test
        manager = GoldenFileManager(golden_dir=tmp_path)
        
        test_content = "Test content\nLine 2\n"
        test_path = "test/test_file.txt"
        
        # Save file
        manager.save_golden_file(test_path, test_content)
        
        # Verify file exists
        assert manager.golden_file_exists(test_path)
        
        # Load and verify content
        loaded_content = manager.load_golden_file(test_path)
        assert loaded_content == test_content


class TestGoldenFileComparison:
    """Test golden file comparison functionality."""
    
    def test_exact_text_match(self):
        """Test exact text matching."""
        from tests.utils.golden_file_utils import compare_text_content
        
        text1 = "Line 1\nLine 2\nLine 3"
        text2 = "Line 1\nLine 2\nLine 3"
        
        is_match, diff = compare_text_content(text1, text2)
        assert is_match
        assert diff == ""
    
    def test_text_mismatch(self):
        """Test text mismatch detection."""
        from tests.utils.golden_file_utils import compare_text_content
        
        text1 = "Line 1\nLine 2\nLine 3"
        text2 = "Line 1\nLine 2 CHANGED\nLine 3"
        
        is_match, diff = compare_text_content(text1, text2)
        assert not is_match
        assert "CHANGED" in diff
    
    def test_json_exact_match(self):
        """Test exact JSON matching."""
        from tests.utils.golden_file_utils import compare_json_content
        
        json1 = {"key1": "value1", "key2": 42}
        json2 = {"key2": 42, "key1": "value1"}  # Different order
        
        is_match, diff = compare_json_content(json1, json2)
        assert is_match
        assert diff == ""
    
    def test_json_mismatch(self):
        """Test JSON mismatch detection."""
        from tests.utils.golden_file_utils import compare_json_content
        
        json1 = {"key1": "value1", "key2": 42}
        json2 = {"key1": "value1", "key2": 43}
        
        is_match, diff = compare_json_content(json1, json2)
        assert not is_match
        assert "42" in diff and "43" in diff


@pytest.mark.integration
class TestGoldenFileIntegration:
    """Integration tests for golden file validation."""
    
    def test_end_to_end_bcftools_workflow(self):
        """Test complete bcftools workflow against golden files."""
        test_file = "sample_test_data/small_valid.vcf"
        
        # Test view
        rc, stdout, stderr = bcftools_view(["-h", test_file])
        assert rc == 0
        assert_golden_match(stdout, "bcftools/view/view_header_small_valid.txt")
        
        # Test stats
        rc, stdout, stderr = bcftools_stats([test_file])
        assert rc == 0
        assert_golden_match(stdout, "bcftools/stats/stats_small_valid.txt")
        
        # Test query
        rc, stdout, stderr = bcftools_query(["-f", "%CHROM\t%POS\t%REF\t%ALT\n", test_file])
        assert rc == 0
        assert_golden_match(stdout, "bcftools/query/query_basic_small_valid.txt")
    
    def test_validation_workflow(self):
        """Test validation workflow against golden files."""
        test_files = [
            "sample_test_data/small_valid.vcf",
            "sample_test_data/sample1.vcf"
        ]
        
        for test_file in test_files:
            file_base = Path(test_file).stem
            validation_result = validate_vcf_file(test_file)
            result_json = json.dumps(validation_result, indent=2, sort_keys=True)
            
            assert_golden_match(
                result_json,
                f"validation/valid_files/validation_{file_base}.json",
                content_type="json"
            ) 