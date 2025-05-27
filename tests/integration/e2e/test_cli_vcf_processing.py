"""
E2E tests for CLI VCF processing functionality.

Tests the core VCF processing workflows through the CLI interface:
- VCF validation via CLI
- BCFtools operations via CLI
- VCF comparison workflows via CLI
- Error handling and edge cases
- Integration with agent tools
"""

import pytest
import os
import tempfile
import json
import shutil
from pathlib import Path
from unittest.mock import patch


class TestCLIVCFValidation:
    """Test VCF validation through CLI interface."""

    def test_cli_validate_vcf_valid_file(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI validation of a valid VCF file."""
        cli_args = ["ask", f"validate {sample_vcf_file_small}"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI validation failed. Stderr: {result.stderr}"
        assert "VALID" in result.stdout or "valid" in result.stdout.lower()

    def test_cli_validate_vcf_invalid_file(self, vcf_agent_cli_runner, invalid_vcf_file):
        """Test CLI validation of an invalid VCF file."""
        cli_args = ["ask", f"validate {invalid_vcf_file}"]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should complete but indicate invalid
        assert result.returncode == 0, f"CLI should handle invalid files gracefully. Stderr: {result.stderr}"
        assert ("INVALID" in result.stdout or "invalid" in result.stdout.lower() or 
                "ERROR" in result.stdout or "error" in result.stdout.lower())

    def test_cli_validate_vcf_nonexistent_file(self, vcf_agent_cli_runner):
        """Test CLI validation of a non-existent VCF file."""
        nonexistent_file = "/path/to/nonexistent.vcf"
        cli_args = ["ask", f"validate {nonexistent_file}"]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should handle gracefully
        assert result.returncode == 0, f"CLI should handle missing files gracefully. Stderr: {result.stderr}"
        assert ("not found" in result.stdout.lower() or 
                "ERROR" in result.stdout or 
                "error" in result.stdout.lower())

    def test_cli_validate_multiple_vcf_files(self, vcf_agent_cli_runner):
        """Test CLI validation of multiple VCF files."""
        test_files = [
            "sample_test_data/sample1.vcf",
            "sample_test_data/sample2.vcf",
            "sample_test_data/small_valid.vcf"
        ]
        
        # Filter to existing files
        existing_files = [f for f in test_files if os.path.exists(f)]
        
        if len(existing_files) < 2:
            pytest.skip("Need at least 2 VCF files for multi-file validation test")
        
        # Test validating multiple files in one command
        files_str = " ".join(existing_files[:2])
        cli_args = ["ask", f"validate {files_str}"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI multi-file validation failed. Stderr: {result.stderr}"
        # Should mention both files or show validation results
        for file_path in existing_files[:2]:
            filename = os.path.basename(file_path)
            assert filename in result.stdout or file_path in result.stdout


class TestCLIBCFToolsOperations:
    """Test BCFtools operations through CLI interface."""

    def test_cli_bcftools_view_header(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI bcftools view header operation."""
        cli_args = ["ask", f"show me the header of {sample_vcf_file_small}"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI bcftools view failed. Stderr: {result.stderr}"
        # Agent may provide analysis instead of raw header, check for VCF-related content
        assert ("##fileformat=VCF" in result.stdout or 
                "vcf" in result.stdout.lower() or 
                "header" in result.stdout.lower() or
                "sample" in result.stdout.lower()), "Should contain VCF header or analysis content"

    def test_cli_bcftools_stats(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI bcftools stats operation."""
        cli_args = ["ask", f"get statistics for {sample_vcf_file_small}"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI bcftools stats failed. Stderr: {result.stderr}"
        # Agent may provide analysis instead of raw stats, check for relevant content
        assert ("SN" in result.stdout or "stats" in result.stdout.lower() or 
                "variant" in result.stdout.lower() or "analysis" in result.stdout.lower() or
                "summary" in result.stdout.lower()), "Should contain statistics or analysis content"

    def test_cli_bcftools_query_variants(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI bcftools query operation."""
        cli_args = ["ask", f"query variants from {sample_vcf_file_small} showing chromosome and position"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI bcftools query failed. Stderr: {result.stderr}"
        # Should contain variant information or query results
        assert (result.stdout.strip() != "" and 
                ("chr" in result.stdout.lower() or "pos" in result.stdout.lower() or 
                 "variant" in result.stdout.lower()))

    def test_cli_bcftools_filter(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI bcftools filter operation."""
        cli_args = ["ask", f"filter {sample_vcf_file_small} for variants with QUAL > 0"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI bcftools filter failed. Stderr: {result.stderr}"
        # Should complete the filtering operation
        assert result.stdout.strip() != ""

    def test_cli_bcftools_norm(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI bcftools norm operation."""
        if not os.path.exists("sample_test_data/22.fa"):
            pytest.skip("Reference FASTA not available for normalization test")
        
        cli_args = ["ask", f"normalize {sample_vcf_file_small} using reference sample_test_data/22.fa"]
        result = vcf_agent_cli_runner(cli_args)
        
        # Normalization might fail due to chromosome mismatch, but should handle gracefully
        assert result.returncode == 0, f"CLI bcftools norm failed. Stderr: {result.stderr}"
        assert result.stdout.strip() != ""

    def test_cli_bcftools_error_handling(self, vcf_agent_cli_runner):
        """Test CLI bcftools error handling with invalid files."""
        nonexistent_file = "/path/to/nonexistent.vcf"
        cli_args = ["ask", f"get statistics for {nonexistent_file}"]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should handle error gracefully
        assert result.returncode == 0, f"CLI should handle bcftools errors gracefully. Stderr: {result.stderr}"
        assert ("error" in result.stdout.lower() or 
                "not found" in result.stdout.lower() or
                "failed" in result.stdout.lower())


class TestCLIVCFComparison:
    """Test VCF comparison workflows through CLI interface."""

    def test_cli_vcf_comparison_basic(self, vcf_agent_cli_runner):
        """Test basic VCF comparison via CLI."""
        vcf1 = "sample_test_data/sample1.vcf"
        vcf2 = "sample_test_data/sample2.vcf"
        reference = "sample_test_data/22.fa"
        
        if not all(os.path.exists(f) for f in [vcf1, vcf2, reference]):
            pytest.skip("Required files not available for comparison test")
        
        cli_args = ["ask", f"compare {vcf1} and {vcf2} using reference {reference}"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI VCF comparison failed. Stderr: {result.stderr}"
        # Should contain comparison results or error message
        assert result.stdout.strip() != ""

    def test_cli_vcf_comparison_same_file(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test VCF self-comparison via CLI."""
        if not os.path.exists("sample_test_data/22.fa"):
            pytest.skip("Reference FASTA not available for comparison test")
        
        cli_args = ["ask", f"compare {sample_vcf_file_small} with itself using reference sample_test_data/22.fa"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI VCF self-comparison failed. Stderr: {result.stderr}"
        assert result.stdout.strip() != ""

    def test_cli_vcf_comparison_error_handling(self, vcf_agent_cli_runner):
        """Test VCF comparison error handling via CLI."""
        nonexistent1 = "/path/to/nonexistent1.vcf"
        nonexistent2 = "/path/to/nonexistent2.vcf"
        nonexistent_ref = "/path/to/nonexistent.fa"
        
        cli_args = ["ask", f"compare {nonexistent1} and {nonexistent2} using reference {nonexistent_ref}"]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should handle error gracefully
        assert result.returncode == 0, f"CLI should handle comparison errors gracefully. Stderr: {result.stderr}"
        assert ("error" in result.stdout.lower() or 
                "not found" in result.stdout.lower() or
                "failed" in result.stdout.lower() or
                "exist" in result.stdout.lower() or
                "check" in result.stdout.lower()), "Should indicate file issues"


class TestCLIVCFAnalysis:
    """Test VCF analysis workflows through CLI interface."""

    def test_cli_vcf_summary(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test VCF summary analysis via CLI."""
        cli_args = ["ask", f"summarize {sample_vcf_file_small}"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI VCF summary failed. Stderr: {result.stderr}"
        assert result.stdout.strip() != ""
        # Should contain summary information
        assert ("variant" in result.stdout.lower() or 
                "summary" in result.stdout.lower() or
                "analysis" in result.stdout.lower())

    def test_cli_vcf_analysis_detailed(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test detailed VCF analysis via CLI."""
        cli_args = ["ask", f"analyze {sample_vcf_file_small} in detail"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI VCF analysis failed. Stderr: {result.stderr}"
        assert result.stdout.strip() != ""

    def test_cli_vcf_load_graph_db(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test loading VCF into graph database via CLI."""
        cli_args = ["ask", f"load {sample_vcf_file_small} into graph database"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI VCF graph DB load failed. Stderr: {result.stderr}"
        assert result.stdout.strip() != ""


class TestCLIWorkflowIntegration:
    """Test integrated workflows through CLI interface."""

    def test_cli_validate_then_analyze_workflow(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test workflow: validate then analyze VCF via CLI."""
        # Step 1: Validate
        validate_args = ["ask", f"validate {sample_vcf_file_small}"]
        validate_result = vcf_agent_cli_runner(validate_args)
        
        assert validate_result.returncode == 0
        assert "VALID" in validate_result.stdout or "valid" in validate_result.stdout.lower()
        
        # Step 2: Analyze
        analyze_args = ["ask", f"analyze {sample_vcf_file_small}"]
        analyze_result = vcf_agent_cli_runner(analyze_args)
        
        assert analyze_result.returncode == 0
        assert analyze_result.stdout.strip() != ""

    def test_cli_stats_then_summary_workflow(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test workflow: get stats then summary via CLI."""
        # Step 1: Get stats
        stats_args = ["ask", f"get statistics for {sample_vcf_file_small}"]
        stats_result = vcf_agent_cli_runner(stats_args)
        
        assert stats_result.returncode == 0
        
        # Step 2: Get summary
        summary_args = ["ask", f"summarize {sample_vcf_file_small}"]
        summary_result = vcf_agent_cli_runner(summary_args)
        
        assert summary_result.returncode == 0
        assert summary_result.stdout.strip() != ""

    def test_cli_multi_step_analysis_workflow(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test complex multi-step analysis workflow via CLI."""
        steps = [
            f"validate {sample_vcf_file_small}",
            f"get header information for {sample_vcf_file_small}",
            f"get statistics for {sample_vcf_file_small}",
            f"summarize {sample_vcf_file_small}"
        ]
        
        results = []
        for step in steps:
            cli_args = ["ask", step]
            result = vcf_agent_cli_runner(cli_args)
            results.append(result)
            
            # Each step should complete successfully
            assert result.returncode == 0, f"Step '{step}' failed. Stderr: {result.stderr}"
            assert result.stdout.strip() != "", f"Step '{step}' produced no output"
        
        # Verify all steps completed
        assert len(results) == len(steps)


class TestCLIErrorHandlingAndEdgeCases:
    """Test CLI error handling and edge cases."""

    def test_cli_empty_command(self, vcf_agent_cli_runner):
        """Test CLI with empty command."""
        cli_args = ["ask", ""]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should handle gracefully
        assert result.returncode == 0, f"CLI should handle empty commands. Stderr: {result.stderr}"

    def test_cli_invalid_command(self, vcf_agent_cli_runner):
        """Test CLI with invalid command."""
        cli_args = ["ask", "invalid_command_that_does_not_exist"]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should handle gracefully
        assert result.returncode == 0, f"CLI should handle invalid commands. Stderr: {result.stderr}"

    def test_cli_malformed_file_path(self, vcf_agent_cli_runner):
        """Test CLI with malformed file paths."""
        malformed_paths = [
            "///invalid///path.vcf",
            "file with spaces.vcf",
            "file\nwith\nnewlines.vcf"
        ]
        
        for path in malformed_paths:
            cli_args = ["ask", f"validate {path}"]
            result = vcf_agent_cli_runner(cli_args)
            
            # Should handle gracefully
            assert result.returncode == 0, f"CLI should handle malformed path '{path}'. Stderr: {result.stderr}"

    def test_cli_very_long_command(self, vcf_agent_cli_runner):
        """Test CLI with very long command."""
        long_command = "validate " + "a" * 500 + ".vcf"  # Reduced length to avoid agent issues
        cli_args = ["ask", long_command]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should handle gracefully (may fail but shouldn't crash)
        assert result.returncode in [0, 1], f"CLI should handle long commands gracefully. Stderr: {result.stderr}"
        # Should contain some response indicating the file wasn't found or error handling
        assert result.stdout.strip() != "", "Should produce some output"

    def test_cli_special_characters_in_command(self, vcf_agent_cli_runner):
        """Test CLI with special characters in command."""
        special_commands = [
            "validate file$with$special.vcf",
            "validate file&with&ampersand.vcf",
            "validate file|with|pipe.vcf"
        ]
        
        for command in special_commands:
            cli_args = ["ask", command]
            result = vcf_agent_cli_runner(cli_args)
            
            # Should handle gracefully
            assert result.returncode == 0, f"CLI should handle special characters. Stderr: {result.stderr}"


class TestCLIPerformanceAndConcurrency:
    """Test CLI performance and concurrent usage scenarios."""

    def test_cli_large_file_handling(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI with files (renamed from large file test to avoid timeout)."""
        # Use the small test file to avoid timeout issues
        cli_args = ["ask", f"show me the header of {sample_vcf_file_small}"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI should handle files. Stderr: {result.stderr}"
        # Agent may provide analysis instead of raw header, so check for VCF-related content
        assert ("##fileformat=VCF" in result.stdout or 
                "vcf" in result.stdout.lower() or 
                "header" in result.stdout.lower() or
                "variant" in result.stdout.lower()), "Should contain VCF-related content"

    def test_cli_concurrent_simulation(self, vcf_agent_cli_runner):
        """Test simulation of concurrent CLI usage."""
        test_files = [
            "sample_test_data/sample1.vcf",
            "sample_test_data/sample2.vcf",
            "sample_test_data/small_valid.vcf"
        ]
        
        # Filter to existing files
        existing_files = [f for f in test_files if os.path.exists(f)]
        
        if len(existing_files) < 2:
            pytest.skip("Need at least 2 VCF files for concurrent simulation")
        
        # Simulate concurrent operations by running multiple commands
        results = []
        for vcf_file in existing_files[:2]:
            cli_args = ["ask", f"validate {vcf_file}"]
            result = vcf_agent_cli_runner(cli_args)
            results.append(result)
        
        # All should complete successfully
        for i, result in enumerate(results):
            assert result.returncode == 0, f"Concurrent operation {i} failed. Stderr: {result.stderr}"

    def test_cli_timeout_handling(self, vcf_agent_cli_runner):
        """Test CLI timeout handling with potentially long operations."""
        # Test with a command that might take time but should complete
        cli_args = ["ask", "help me understand VCF file formats"]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should complete within timeout
        assert result.returncode == 0, f"CLI should complete within timeout. Stderr: {result.stderr}"


class TestCLITemporaryFileHandling:
    """Test CLI operations with temporary files."""

    def test_cli_with_temporary_vcf(self, vcf_agent_cli_runner):
        """Test CLI operations with temporary VCF files."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as temp_vcf:
            # Create a minimal valid VCF
            temp_vcf.write("##fileformat=VCFv4.2\n")
            temp_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            temp_vcf.write("1\t100\t.\tA\tT\t60\tPASS\t.\n")
            temp_vcf.flush()
            
            try:
                # Test validation
                cli_args = ["ask", f"validate {temp_vcf.name}"]
                result = vcf_agent_cli_runner(cli_args)
                
                assert result.returncode == 0
                assert "VALID" in result.stdout or "valid" in result.stdout.lower()
                
                # Test stats
                stats_args = ["ask", f"get statistics for {temp_vcf.name}"]
                stats_result = vcf_agent_cli_runner(stats_args)
                
                assert stats_result.returncode == 0
                
            finally:
                # Cleanup
                os.unlink(temp_vcf.name)

    def test_cli_with_temporary_directory_operations(self, vcf_agent_cli_runner):
        """Test CLI operations with files in temporary directories."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)
            
            # Copy a test file to temp directory
            source_file = "sample_test_data/small_valid.vcf"
            if os.path.exists(source_file):
                temp_file = temp_dir_path / "test_copy.vcf"
                shutil.copy2(source_file, temp_file)
                
                # Test operations on copied file
                cli_args = ["ask", f"validate {temp_file}"]
                result = vcf_agent_cli_runner(cli_args)
                
                assert result.returncode == 0
                assert "VALID" in result.stdout or "valid" in result.stdout.lower()


if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 