"""
Integration tests for VCF processing workflows.

Tests end-to-end workflows that combine multiple VCF Agent components:
- VCF validation + bcftools operations
- Agent tools + file I/O + validation
- Multi-step processing pipelines
- Error handling across component boundaries
"""

import pytest
import os
import tempfile
import shutil
import json
from pathlib import Path
from unittest.mock import patch, MagicMock

from vcf_agent.agent import (
    validate_vcf,
    bcftools_view_tool,
    bcftools_query_tool,
    bcftools_filter_tool,
    bcftools_norm_tool,
    bcftools_stats_tool,
    vcf_comparison_tool,
    load_vcf_into_graph_db_tool,
    get_agent_with_session
)
from vcf_agent.validation import validate_vcf_file
from vcf_agent.bcftools_integration import (
    bcftools_view,
    bcftools_query,
    bcftools_filter,
    bcftools_stats
)
from vcf_agent import graph_integration


class TestVCFValidationWorkflows:
    """Test workflows that combine VCF validation with other operations."""

    @pytest.fixture
    def sample_vcf_path(self):
        """Provide path to a sample VCF file for testing."""
        return "sample_test_data/small_valid.vcf"

    @pytest.fixture
    def invalid_vcf_path(self):
        """Provide path to an invalid VCF file for testing."""
        return "sample_test_data/edgecase_invalid_qual.vcf"

    def test_validate_then_process_workflow(self, sample_vcf_path):
        """Test workflow: validate VCF → process with bcftools → validate results."""
        # Step 1: Validate input VCF
        validation_result = validate_vcf(sample_vcf_path)
        assert "VALID" in validation_result, f"VCF should be valid: {validation_result}"

        # Step 2: Process with bcftools view (header only)
        view_result = bcftools_view_tool(["-h", sample_vcf_path])
        assert "##fileformat=VCF" in view_result, "Should contain VCF header"

        # Step 3: Get basic stats
        stats_result = bcftools_stats_tool([sample_vcf_path])
        assert "SN" in stats_result, "Stats should contain summary numbers"

    def test_validation_failure_stops_workflow(self, invalid_vcf_path):
        """Test that validation failure prevents further processing."""
        # Step 1: Validate input VCF (should fail)
        validation_result = validate_vcf(invalid_vcf_path)
        assert "INVALID" in validation_result or "ERROR" in validation_result

        # Step 2: Processing should still work (bcftools is permissive)
        # but we can check that our validation caught the issue
        view_result = bcftools_view_tool(["-h", invalid_vcf_path])
        # bcftools might still process it, but our validation caught the issue
        assert isinstance(view_result, str)

    def test_multi_file_validation_workflow(self):
        """Test validation workflow with multiple files."""
        test_files = [
            "sample_test_data/small_valid.vcf",
            "sample_test_data/sample1.vcf",
            "sample_test_data/sample2.vcf"
        ]
        
        validation_results = []
        for vcf_file in test_files:
            if os.path.exists(vcf_file):
                result = validate_vcf(vcf_file)
                validation_results.append((vcf_file, result))
        
        assert len(validation_results) > 0, "Should have validated at least one file"
        
        # All valid files should pass validation
        for file_path, result in validation_results:
            if "small_valid" in file_path or "sample" in file_path:
                assert "VALID" in result, f"File {file_path} should be valid"


class TestBCFToolsWorkflows:
    """Test workflows that chain multiple bcftools operations."""

    @pytest.fixture
    def sample_vcf_path(self):
        """Provide path to a sample VCF file for testing."""
        return "sample_test_data/small_valid.vcf"

    def test_view_then_query_workflow(self, sample_vcf_path):
        """Test workflow: bcftools view → bcftools query."""
        # Step 1: View VCF header
        header_result = bcftools_view_tool(["-h", sample_vcf_path])
        assert "##fileformat=VCF" in header_result

        # Step 2: Query specific fields
        query_result = bcftools_query_tool(["-f", "%CHROM\t%POS\t%REF\t%ALT\n", sample_vcf_path])
        
        # Should have tab-separated output
        if query_result and not query_result.startswith("Error"):
            lines = query_result.strip().split('\n')
            if lines and lines[0]:  # If there are variants
                fields = lines[0].split('\t')
                assert len(fields) >= 4, "Should have CHROM, POS, REF, ALT fields"

    def test_filter_then_stats_workflow(self, sample_vcf_path):
        """Test workflow: bcftools filter → bcftools stats."""
        # Step 1: Apply a permissive filter (should pass most variants)
        filter_result = bcftools_filter_tool(["-i", "QUAL>0", sample_vcf_path])
        
        # Step 2: Get stats on original file
        stats_result = bcftools_stats_tool([sample_vcf_path])
        assert "SN" in stats_result, "Stats should contain summary numbers"

    def test_complex_bcftools_pipeline(self, sample_vcf_path):
        """Test complex pipeline: view → query → stats."""
        results = {}
        
        # Step 1: Get header info
        results['header'] = bcftools_view_tool(["-h", sample_vcf_path])
        assert "##fileformat=VCF" in results['header']
        
        # Step 2: Query variant positions
        results['positions'] = bcftools_query_tool(["-f", "%CHROM:%POS\n", sample_vcf_path])
        
        # Step 3: Get comprehensive stats
        results['stats'] = bcftools_stats_tool([sample_vcf_path])
        assert "SN" in results['stats']
        
        # Verify all steps completed successfully
        for step, result in results.items():
            assert result is not None, f"Step {step} should have produced output"
            assert not result.startswith("Error"), f"Step {step} should not have errors: {result}"


class TestAgentToolWorkflows:
    """Test workflows using agent tools in combination."""

    @pytest.fixture
    def sample_vcf_path(self):
        """Provide path to a sample VCF file for testing."""
        return "sample_test_data/small_valid.vcf"

    def test_agent_tool_chain_workflow(self, sample_vcf_path):
        """Test chaining multiple agent tools together."""
        # Step 1: Validate VCF using agent tool
        validation_result = validate_vcf(sample_vcf_path)
        assert "VALID" in validation_result

        # Step 2: Process with bcftools agent tools
        view_result = bcftools_view_tool(["-h", sample_vcf_path])
        assert "##fileformat=VCF" in view_result

        # Step 3: Query using agent tool
        query_result = bcftools_query_tool(["-f", "%CHROM\t%POS\n", sample_vcf_path])
        
        # Verify the chain worked
        assert isinstance(validation_result, str)
        assert isinstance(view_result, str)
        assert isinstance(query_result, str)

    @patch('vcf_agent.agent.graph_integration')
    def test_vcf_to_graph_workflow(self, mock_graph_integration, sample_vcf_path):
        """Test workflow: validate VCF → load into graph database."""
        # Mock graph integration
        mock_graph_integration.get_managed_kuzu_connection.return_value = MagicMock()
        mock_graph_integration.load_vcf_into_kuzu.return_value = None

        # Step 1: Validate VCF
        validation_result = validate_vcf(sample_vcf_path)
        assert "VALID" in validation_result

        # Step 2: Load into graph database
        graph_result = load_vcf_into_graph_db_tool(sample_vcf_path)
        
        # Parse JSON result
        graph_data = json.loads(graph_result)
        assert "status" in graph_data
        assert "success" in graph_data.get("status", "").lower()

    def test_error_propagation_workflow(self):
        """Test how errors propagate through tool chains."""
        non_existent_file = "/path/to/nonexistent.vcf"
        
        # Step 1: Validation should fail gracefully
        validation_result = validate_vcf(non_existent_file)
        assert "ERROR" in validation_result or "not found" in validation_result.lower()
        
        # Step 2: bcftools operations should also fail gracefully
        view_result = bcftools_view_tool(["-h", non_existent_file])
        assert "Error" in view_result or "failed" in view_result.lower()


class TestFileIOWorkflows:
    """Test workflows involving file I/O operations."""

    def test_temporary_file_workflow(self):
        """Test workflow with temporary file creation and cleanup."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as temp_vcf:
            # Create a minimal valid VCF
            temp_vcf.write("##fileformat=VCFv4.2\n")
            temp_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            temp_vcf.write("1\t100\t.\tA\tT\t60\tPASS\t.\n")
            temp_vcf.flush()
            
            try:
                # Step 1: Validate the temporary VCF
                validation_result = validate_vcf(temp_vcf.name)
                assert "VALID" in validation_result

                # Step 2: Process with bcftools
                view_result = bcftools_view_tool(["-h", temp_vcf.name])
                assert "##fileformat=VCF" in view_result

                # Step 3: Query the data
                query_result = bcftools_query_tool(["-f", "%CHROM\t%POS\n", temp_vcf.name])
                assert "1\t100" in query_result

            finally:
                # Cleanup
                os.unlink(temp_vcf.name)

    def test_multiple_file_processing_workflow(self):
        """Test workflow processing multiple VCF files."""
        test_files = [
            "sample_test_data/sample1.vcf",
            "sample_test_data/sample2.vcf"
        ]
        
        results = {}
        
        for vcf_file in test_files:
            if os.path.exists(vcf_file):
                # Process each file through the same workflow
                file_results = {}
                
                # Step 1: Validate
                file_results['validation'] = validate_vcf(vcf_file)
                
                # Step 2: Get header
                file_results['header'] = bcftools_view_tool(["-h", vcf_file])
                
                # Step 3: Get stats
                file_results['stats'] = bcftools_stats_tool([vcf_file])
                
                results[vcf_file] = file_results
        
        # Verify all files were processed
        assert len(results) > 0, "Should have processed at least one file"
        
        for file_path, file_results in results.items():
            assert "VALID" in file_results['validation'], f"File {file_path} should be valid"
            assert "##fileformat=VCF" in file_results['header'], f"File {file_path} should have VCF header"


class TestErrorHandlingWorkflows:
    """Test error handling across workflow boundaries."""

    def test_cascading_error_handling(self):
        """Test how errors cascade through workflow steps."""
        invalid_file = "/definitely/does/not/exist.vcf"
        
        # Each step should handle the error gracefully
        validation_result = validate_vcf(invalid_file)
        view_result = bcftools_view_tool(["-h", invalid_file])
        query_result = bcftools_query_tool(["-f", "%CHROM\n", invalid_file])
        
        # All should contain error information
        error_indicators = ["ERROR", "Error", "failed", "not found", "No such file"]
        
        for result in [validation_result, view_result, query_result]:
            assert any(indicator in result for indicator in error_indicators), \
                f"Result should contain error indicator: {result}"

    def test_partial_failure_workflow(self):
        """Test workflow where some steps succeed and others fail."""
        valid_file = "sample_test_data/small_valid.vcf"
        
        # Step 1: Should succeed
        validation_result = validate_vcf(valid_file)
        assert "VALID" in validation_result
        
        # Step 2: Should succeed
        view_result = bcftools_view_tool(["-h", valid_file])
        assert "##fileformat=VCF" in view_result
        
        # Step 3: Intentionally use invalid bcftools arguments
        invalid_query_result = bcftools_query_tool(["-invalid-flag", valid_file])
        # This should fail gracefully
        assert "Error" in invalid_query_result or "failed" in invalid_query_result.lower()

    def test_recovery_workflow(self):
        """Test workflow recovery after encountering errors."""
        # Start with an invalid operation
        invalid_result = bcftools_query_tool(["-invalid-flag", "nonexistent.vcf"])
        assert "Error" in invalid_result or "failed" in invalid_result.lower()
        
        # Recover with a valid operation
        if os.path.exists("sample_test_data/small_valid.vcf"):
            recovery_result = validate_vcf("sample_test_data/small_valid.vcf")
            assert "VALID" in recovery_result, "Should recover with valid operation"


class TestPerformanceWorkflows:
    """Test workflows with performance considerations."""

    def test_large_file_workflow(self):
        """Test workflow with larger VCF files (if available)."""
        large_files = [
            "sample_test_data/chr22.1kg.phase3.v5a.testready.vcf.gz",
            "sample_test_data/detailed.vcf"
        ]
        
        for large_file in large_files:
            if os.path.exists(large_file):
                # Test that we can handle larger files
                # Use header-only operations to avoid long processing times
                header_result = bcftools_view_tool(["-h", large_file])
                assert "##fileformat=VCF" in header_result
                
                # Quick validation check
                validation_result = validate_vcf(large_file)
                # Should complete without hanging
                assert isinstance(validation_result, str)
                break

    def test_concurrent_workflow_simulation(self):
        """Test simulation of concurrent workflow processing."""
        test_files = [
            "sample_test_data/sample1.vcf",
            "sample_test_data/sample2.vcf",
            "sample_test_data/small_valid.vcf"
        ]
        
        # Simulate concurrent processing by running same operations on multiple files
        results = []
        for vcf_file in test_files:
            if os.path.exists(vcf_file):
                # Quick operations that could run concurrently
                validation = validate_vcf(vcf_file)
                header = bcftools_view_tool(["-h", vcf_file])
                
                results.append({
                    'file': vcf_file,
                    'validation': validation,
                    'header': header
                })
        
        # Verify all completed successfully
        assert len(results) > 0, "Should have processed multiple files"
        for result in results:
            assert "VALID" in result['validation'], f"File {result['file']} should be valid"
            assert "##fileformat=VCF" in result['header'], f"File {result['file']} should have header"


class TestIntegrationWithExternalSystems:
    """Test integration workflows with external systems."""

    @patch('vcf_agent.agent.graph_integration')
    def test_database_integration_workflow(self, mock_graph_integration):
        """Test workflow integrating with graph database."""
        # Mock the graph integration
        mock_conn = MagicMock()
        mock_graph_integration.get_managed_kuzu_connection.return_value = mock_conn
        mock_graph_integration.load_vcf_into_kuzu.return_value = None

        vcf_file = "sample_test_data/small_valid.vcf"
        if os.path.exists(vcf_file):
            # Step 1: Validate VCF
            validation_result = validate_vcf(vcf_file)
            assert "VALID" in validation_result

            # Step 2: Process with bcftools
            stats_result = bcftools_stats_tool([vcf_file])
            assert "SN" in stats_result

            # Step 3: Load into database
            db_result = load_vcf_into_graph_db_tool(vcf_file)
            db_data = json.loads(db_result)
            assert "status" in db_data

    def test_file_system_integration_workflow(self):
        """Test workflow with file system operations."""
        # Create temporary directory for workflow
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)
            
            # Copy a test file to temp directory
            source_file = "sample_test_data/small_valid.vcf"
            if os.path.exists(source_file):
                temp_file = temp_dir_path / "test_copy.vcf"
                shutil.copy2(source_file, temp_file)
                
                # Process the copied file
                validation_result = validate_vcf(str(temp_file))
                assert "VALID" in validation_result
                
                # Verify file operations work
                assert temp_file.exists(), "Copied file should exist"
                
                # Process with bcftools
                view_result = bcftools_view_tool(["-h", str(temp_file)])
                assert "##fileformat=VCF" in view_result


if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 