"""
Tests for AI-powered VCF analysis functionality.

This module tests the AI integration features including:
- vcf_analysis_summary_tool
- Enhanced vcf_summarization_tool  
- ai_vcf_comparison_tool
- LLM integration and fallback mechanisms
"""

import pytest
import json
import os
import tempfile
from unittest.mock import patch, MagicMock

from src.vcf_agent.agent import (
    vcf_analysis_summary_tool,
    vcf_summarization_tool,
    ai_vcf_comparison_tool,
    run_llm_analysis_task,
    _generate_basic_vcf_summary,
    _create_basic_summary_result
)


class TestVCFAnalysisSummaryTool:
    """Test the AI-powered VCF analysis summary tool."""
    
    def test_vcf_analysis_summary_tool_success(self):
        """Test successful AI analysis with LLM."""
        test_file = "sample_data/minimal.vcf.gz"
        
        # Mock successful LLM response
        mock_llm_response = {
            "variant_count": 5,
            "variant_types": {"SNP": 3, "INDEL": 2},
            "sample_statistics": {"sample1": {"mean_depth": 25.5, "het_ratio": 0.6}},
            "notable_patterns": ["High heterozygosity in sample1", "Quality distribution normal"]
        }
        
        with patch('src.vcf_agent.agent.run_llm_analysis_task') as mock_llm, \
             patch('src.vcf_agent.bcftools_integration.bcftools_stats') as mock_stats, \
             patch('src.vcf_agent.bcftools_integration.bcftools_query') as mock_query:
            
            # Setup mocks
            mock_llm.return_value = json.dumps(mock_llm_response)
            mock_stats.return_value = (0, "SN\t0\tnumber of records:\t5\nSN\t0\tnumber of SNPs:\t3", "")
            mock_query.return_value = (0, "sample1", "")
            
            # Execute
            result = vcf_analysis_summary_tool(test_file)
            
            # Verify
            result_data = json.loads(result)
            assert result_data["variant_count"] == 5
            assert result_data["variant_types"]["SNP"] == 3
            assert "sample1" in result_data["sample_statistics"]
            assert len(result_data["notable_patterns"]) > 0
            
            # Verify LLM was called with correct parameters
            mock_llm.assert_called_once()
            call_args = mock_llm.call_args
            assert call_args[1]["task"] == "vcf_analysis_summary_v1"
            assert test_file in call_args[1]["file_paths"]
    
    def test_vcf_analysis_summary_tool_llm_failure_fallback(self):
        """Test fallback to basic analysis when LLM fails."""
        test_file = "sample_data/minimal.vcf.gz"
        
        with patch('src.vcf_agent.agent.run_llm_analysis_task') as mock_llm, \
             patch('src.vcf_agent.bcftools_integration.bcftools_stats') as mock_stats, \
             patch('src.vcf_agent.bcftools_integration.bcftools_query') as mock_query:
            
            # Setup mocks - LLM returns error
            mock_llm.return_value = json.dumps({"error": "LLM service unavailable"})
            mock_stats.return_value = (0, "SN\t0\tnumber of records:\t3\nSN\t0\tnumber of SNPs:\t2", "")
            mock_query.return_value = (0, "sample1", "")
            
            # Execute
            result = vcf_analysis_summary_tool(test_file)
            
            # Verify fallback worked
            result_data = json.loads(result)
            assert result_data["variant_count"] == 3
            assert result_data["variant_types"]["SNP"] == 2
            assert "analysis_method" in result_data
            assert result_data["analysis_method"] == "basic_fallback"
    
    def test_vcf_analysis_summary_tool_file_not_found(self):
        """Test handling of non-existent files."""
        result = vcf_analysis_summary_tool("nonexistent.vcf")
        
        result_data = json.loads(result)
        assert "error" in result_data
        assert "Analysis failed" in result_data["error"]
    
    def test_generate_basic_vcf_summary(self):
        """Test the basic VCF summary generation function."""
        stats_output = """
SN	0	number of records:	10
SN	0	number of SNPs:	7
SN	0	number of indels:	3
"""
        samples = ["sample1", "sample2"]
        
        result = _generate_basic_vcf_summary("test.vcf", stats_output, samples)
        result_data = json.loads(result)
        
        assert result_data["variant_count"] == 10
        assert result_data["variant_types"]["SNP"] == 7
        assert result_data["variant_types"]["INDEL"] == 3
        assert len(result_data["sample_statistics"]) == 2
        assert "sample1" in result_data["sample_statistics"]
        assert "sample2" in result_data["sample_statistics"]


class TestVCFSummarizationTool:
    """Test the enhanced VCF summarization tool."""
    
    def test_vcf_summarization_tool_llm_success(self):
        """Test successful LLM-powered summarization."""
        test_file = "sample_data/minimal.vcf.gz"
        
        mock_llm_response = {
            "variant_count": 8,
            "variant_types": {"SNP": 5, "INDEL": 3},
            "sample_statistics": {"sample1": {"mean_depth": 30.2, "het_ratio": 0.55}},
            "notable_patterns": ["Consistent quality scores", "Even distribution across chromosome"]
        }
        
        with patch('src.vcf_agent.agent.run_llm_analysis_task') as mock_llm, \
             patch('src.vcf_agent.bcftools_integration.bcftools_stats') as mock_stats, \
             patch('src.vcf_agent.bcftools_integration.bcftools_query') as mock_query:
            
            # Setup mocks
            mock_llm.return_value = json.dumps(mock_llm_response)
            mock_stats.return_value = (0, "stats output", "")
            mock_query.return_value = (0, "sample1", "")
            
            # Execute
            result = vcf_summarization_tool(test_file)
            
            # Verify
            result_data = json.loads(result)
            assert result_data["variant_count"] == 8
            assert result_data["variant_types"]["SNP"] == 5
            assert "sample1" in result_data["sample_statistics"]
            
            # Verify LLM was called correctly
            mock_llm.assert_called_once()
            call_args = mock_llm.call_args
            assert call_args[1]["task"] == "vcf_summarization_v1"
    
    def test_vcf_summarization_tool_basic_fallback(self):
        """Test fallback to basic analysis when LLM fails."""
        test_file = "sample_data/minimal.vcf.gz"
        
        with patch('src.vcf_agent.agent.run_llm_analysis_task') as mock_llm, \
             patch('src.vcf_agent.bcftools_integration.bcftools_stats') as mock_stats, \
             patch('src.vcf_agent.bcftools_integration.bcftools_query') as mock_query, \
             patch('src.vcf_agent.vcf_utils.extract_variant_summary') as mock_extract:
            
            # Setup mocks - LLM fails, fall back to basic
            mock_llm.side_effect = Exception("LLM service error")
            mock_stats.return_value = (1, "", "bcftools failed")  # bcftools fails too
            mock_query.return_value = (0, "sample1", "")
            mock_extract.return_value = {
                "variant_count": 4,
                "variant_types": {"SNP": 4},
                "samples": ["sample1"]
            }
            
            # Execute
            result = vcf_summarization_tool(test_file)
            
            # Verify fallback worked
            result_data = json.loads(result)
            assert result_data["variant_count"] == 4
            assert result_data["variant_types"]["SNP"] == 4
            assert "analysis_method" in result_data
            assert result_data["analysis_method"] == "basic"
    
    def test_create_basic_summary_result(self):
        """Test the basic summary result creation function."""
        summary = {
            "variant_count": 6,
            "variant_types": {"SNP": 4, "INDEL": 2},
            "samples": ["sample1", "sample2"]
        }
        
        result = _create_basic_summary_result(summary, error_note="Test error")
        result_data = json.loads(result)
        
        assert result_data["variant_count"] == 6
        assert result_data["variant_types"]["SNP"] == 4
        assert len(result_data["sample_statistics"]) == 2
        assert any("Test error" in pattern for pattern in result_data["notable_patterns"])
        assert result_data["analysis_method"] == "basic"


class TestAIVCFComparisonTool:
    """Test the AI-powered VCF comparison tool."""
    
    def test_ai_vcf_comparison_tool_success(self):
        """Test successful AI-powered comparison."""
        file1 = "sample_data/file1.vcf.gz"
        file2 = "sample_data/file2.vcf.gz"
        reference = "sample_data/reference.fa"
        
        # Mock basic comparison result
        basic_comparison = {
            "concordant_variant_count": 10,
            "discordant_variant_count": 5,
            "unique_to_file_1": ["chr1:100:A:T"],
            "unique_to_file_2": ["chr1:200:G:C"],
            "quality_metrics": {"mean_qual": 30.5}
        }
        
        # Mock LLM enhancement
        llm_enhancement = {
            "ai_insights": ["High concordance suggests good quality", "Discordant variants cluster in repetitive regions"],
            "quality_assessment": "Good",
            "recommendations": ["Consider filtering low-quality variants"]
        }
        
        with patch('src.vcf_agent.agent.vcf_comparison_tool') as mock_basic, \
             patch('src.vcf_agent.agent.run_llm_analysis_task') as mock_llm, \
             patch('src.vcf_agent.bcftools_integration.bcftools_stats') as mock_stats:
            
            # Setup mocks
            mock_basic.return_value = json.dumps(basic_comparison)
            mock_llm.return_value = json.dumps(llm_enhancement)
            mock_stats.return_value = (0, "stats output", "")
            
            # Execute
            result = ai_vcf_comparison_tool(file1, file2, reference)
            
            # Verify
            result_data = json.loads(result)
            assert result_data["concordant_variant_count"] == 10
            assert result_data["discordant_variant_count"] == 5
            assert "ai_insights" in result_data
            assert result_data["analysis_method"] == "ai_powered"
            assert len(result_data["ai_insights"]) > 0
            
            # Verify LLM was called correctly
            mock_llm.assert_called_once()
            call_args = mock_llm.call_args
            assert call_args[1]["task"] == "vcf_comparison_v1"
            assert file1 in call_args[1]["file_paths"]
            assert file2 in call_args[1]["file_paths"]
    
    def test_ai_vcf_comparison_tool_llm_failure_fallback(self):
        """Test fallback to basic comparison when LLM fails."""
        file1 = "sample_data/file1.vcf.gz"
        file2 = "sample_data/file2.vcf.gz"
        reference = "sample_data/reference.fa"
        
        basic_comparison = {
            "concordant_variant_count": 8,
            "discordant_variant_count": 3,
            "unique_to_file_1": [],
            "unique_to_file_2": [],
            "quality_metrics": {}
        }
        
        with patch('src.vcf_agent.agent.vcf_comparison_tool') as mock_basic, \
             patch('src.vcf_agent.agent.run_llm_analysis_task') as mock_llm, \
             patch('src.vcf_agent.bcftools_integration.bcftools_stats') as mock_stats:
            
            # Setup mocks - LLM fails
            mock_basic.return_value = json.dumps(basic_comparison)
            mock_llm.return_value = json.dumps({"error": "LLM analysis failed"})
            mock_stats.return_value = (0, "stats output", "")
            
            # Execute
            result = ai_vcf_comparison_tool(file1, file2, reference)
            
            # Verify fallback worked
            result_data = json.loads(result)
            assert result_data["concordant_variant_count"] == 8
            assert result_data["analysis_method"] == "basic_fallback"
            assert "note" in result_data
            assert "AI analysis failed" in result_data["note"]
    
    def test_ai_vcf_comparison_tool_basic_comparison_failure(self):
        """Test handling when basic comparison fails."""
        file1 = "nonexistent1.vcf"
        file2 = "nonexistent2.vcf"
        reference = "nonexistent.fa"
        
        with patch('src.vcf_agent.agent.vcf_comparison_tool') as mock_basic:
            # Setup mock - basic comparison fails
            mock_basic.return_value = json.dumps({"error": "Files not found"})
            
            # Execute
            result = ai_vcf_comparison_tool(file1, file2, reference)
            
            # Verify error handling
            result_data = json.loads(result)
            assert "error" in result_data
            assert "Comparison failed" in result_data["error"]


class TestRunLLMAnalysisTask:
    """Test the core LLM analysis task runner."""
    
    def test_run_llm_analysis_task_success(self):
        """Test successful LLM analysis task execution."""
        with patch('src.vcf_agent.agent.get_agent_with_session') as mock_agent_factory, \
             patch('src.vcf_agent.agent.get_prompt_for_task') as mock_prompt:
            
            # Setup mocks
            mock_agent = MagicMock()
            mock_agent.return_value = '{"result": "success", "data": "test"}'
            mock_agent_factory.return_value = mock_agent
            mock_prompt.return_value = "Test prompt"
            
            # Execute
            result = run_llm_analysis_task(
                task="test_task",
                file_paths=["test.vcf"],
                extra_context={"test": "context"}
            )
            
            # Verify
            result_data = json.loads(result)
            assert result_data["result"] == "success"
            assert result_data["data"] == "test"
            
            # Verify agent was called
            mock_agent.assert_called_once_with("Test prompt")
    
    def test_run_llm_analysis_task_invalid_json_response(self):
        """Test handling of invalid JSON response from LLM."""
        with patch('src.vcf_agent.agent.get_agent_with_session') as mock_agent_factory, \
             patch('src.vcf_agent.agent.get_prompt_for_task') as mock_prompt:
            
            # Setup mocks - agent returns invalid JSON
            mock_agent = MagicMock()
            mock_agent.return_value = "This is not valid JSON"
            mock_agent_factory.return_value = mock_agent
            mock_prompt.return_value = "Test prompt"
            
            # Execute
            result = run_llm_analysis_task(
                task="test_task",
                file_paths=["test.vcf"]
            )
            
            # Verify error handling
            result_data = json.loads(result)
            assert "error" in result_data
            assert "not valid JSON" in result_data["error"]
    
    def test_run_llm_analysis_task_unknown_task(self):
        """Test handling of unknown task."""
        with patch('src.vcf_agent.agent.get_prompt_for_task') as mock_prompt:
            
            # Setup mock - prompt generation fails
            mock_prompt.return_value = None
            
            # Execute
            result = run_llm_analysis_task(
                task="unknown_task",
                file_paths=["test.vcf"]
            )
            
            # Verify error handling
            result_data = json.loads(result)
            assert "error" in result_data
            assert "Unknown task" in result_data["error"]


class TestIntegration:
    """Integration tests for AI analysis functionality."""
    
    def test_ai_tools_integration_with_real_file(self):
        """Test AI tools with a real VCF file (if available)."""
        test_file = "sample_data/minimal.vcf.gz"
        
        if not os.path.exists(test_file):
            pytest.skip(f"Test file {test_file} not available")
        
        # Test with mocked LLM to avoid external dependencies
        with patch('src.vcf_agent.agent.run_llm_analysis_task') as mock_llm:
            mock_llm.return_value = json.dumps({
                "variant_count": 1,
                "variant_types": {"SNP": 1},
                "sample_statistics": {"test_sample": {"mean_depth": 30.0, "het_ratio": 0.5}},
                "notable_patterns": ["Test pattern"]
            })
            
            # Test analysis summary tool
            result = vcf_analysis_summary_tool(test_file)
            result_data = json.loads(result)
            assert "variant_count" in result_data
            assert "variant_types" in result_data
            
            # Test summarization tool
            result = vcf_summarization_tool(test_file)
            result_data = json.loads(result)
            assert "variant_count" in result_data
            assert "sample_statistics" in result_data
    
    def test_error_resilience(self):
        """Test that AI tools are resilient to various error conditions."""
        # Test with non-existent file
        result = vcf_analysis_summary_tool("nonexistent.vcf")
        result_data = json.loads(result)
        assert "error" in result_data
        
        # Test with invalid file path
        result = vcf_summarization_tool("")
        result_data = json.loads(result)
        assert "error" in result_data or "variant_count" in result_data  # Should handle gracefully


if __name__ == "__main__":
    pytest.main([__file__]) 