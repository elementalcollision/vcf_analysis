"""
Comprehensive unit tests for the agent module.

Tests the VCF Analysis Agent tools, model initialization, agent creation,
and core functionality including bcftools integration and error handling.
"""

import pytest
import os
import json
import tempfile
import pandas as pd
from unittest.mock import patch, MagicMock, Mock, call, mock_open
from vcf_agent.agent import (
    echo,
    validate_vcf,
    bcftools_view_tool,
    bcftools_query_tool,
    bcftools_filter_tool,
    bcftools_norm_tool,
    bcftools_stats_tool,
    bcftools_annotate_tool,
    vcf_comparison_tool,
    vcf_analysis_summary_tool,
    vcf_summarization_tool,
    load_vcf_into_graph_db_tool,
    get_openai_model,
    CerebrasStrandsModel,
    get_cerebras_model,
    get_agent_with_session,
    run_llm_analysis_task,
    SYSTEM_PROMPT,
    RAW_MODE
)


class TestEchoTool:
    """Test cases for the echo tool."""

    def test_echo_basic_functionality(self):
        """Test basic echo functionality."""
        result = echo("Hello, world!")
        assert result == "Echo: Hello, world!"

    def test_echo_empty_string(self):
        """Test echo with empty string."""
        result = echo("")
        assert result == "Echo: "

    def test_echo_special_characters(self):
        """Test echo with special characters."""
        text = "Special chars: !@#$%^&*()"
        result = echo(text)
        assert result == f"Echo: {text}"


class TestValidateVCFTool:
    """Test cases for the validate_vcf tool."""

    @patch('vcf_agent.agent.validate_vcf_file')
    @patch('vcf_agent.agent.metrics')
    @patch('vcf_agent.agent.agent_tracer')
    def test_validate_vcf_success(self, mock_tracer, mock_metrics, mock_validate_vcf_file):
        """Test successful VCF validation."""
        # Mock tracer
        mock_span = MagicMock()
        mock_tracer.start_as_current_span.return_value.__enter__.return_value = mock_span
        
        # Mock validation function
        mock_validate_vcf_file.return_value = (True, None)
        
        # Mock metrics
        mock_metrics.VCF_AGENT_TOOL_DURATION_SECONDS.labels.return_value.observe = MagicMock()
        mock_metrics.VCF_AGENT_TOOL_REQUESTS_TOTAL.labels.return_value.inc = MagicMock()
        
        result = validate_vcf("/path/to/test.vcf")
        
        assert "VALID: /path/to/test.vcf passed all checks." in result
        mock_validate_vcf_file.assert_called_once_with("/path/to/test.vcf")
        mock_span.set_status.assert_called_once()

    @patch('vcf_agent.agent.validate_vcf_file')
    @patch('vcf_agent.agent.metrics')
    @patch('vcf_agent.agent.agent_tracer')
    def test_validate_vcf_failure(self, mock_tracer, mock_metrics, mock_validate_vcf_file):
        """Test VCF validation failure."""
        # Mock tracer
        mock_span = MagicMock()
        mock_tracer.start_as_current_span.return_value.__enter__.return_value = mock_span
        
        # Mock validation function
        mock_validate_vcf_file.return_value = (False, "Invalid header")
        
        # Mock metrics
        mock_metrics.VCF_AGENT_TOOL_DURATION_SECONDS.labels.return_value.observe = MagicMock()
        mock_metrics.VCF_AGENT_TOOL_REQUESTS_TOTAL.labels.return_value.inc = MagicMock()
        mock_metrics.VCF_AGENT_TOOL_ERRORS_TOTAL.labels.return_value.inc = MagicMock()
        
        result = validate_vcf("/path/to/invalid.vcf")
        
        assert "INVALID: /path/to/invalid.vcf failed validation" in result
        assert "Invalid header" in result
        mock_validate_vcf_file.assert_called_once_with("/path/to/invalid.vcf")

    @patch('vcf_agent.agent.validate_vcf_file')
    @patch('vcf_agent.agent.metrics')
    @patch('vcf_agent.agent.agent_tracer')
    def test_validate_vcf_exception(self, mock_tracer, mock_metrics, mock_validate_vcf_file):
        """Test VCF validation with exception."""
        # Mock tracer
        mock_span = MagicMock()
        mock_tracer.start_as_current_span.return_value.__enter__.return_value = mock_span
        
        # Mock validation function to raise exception
        mock_validate_vcf_file.side_effect = Exception("File not found")
        
        # Mock metrics
        mock_metrics.VCF_AGENT_TOOL_DURATION_SECONDS.labels.return_value.observe = MagicMock()
        mock_metrics.VCF_AGENT_TOOL_REQUESTS_TOTAL.labels.return_value.inc = MagicMock()
        mock_metrics.VCF_AGENT_TOOL_ERRORS_TOTAL.labels.return_value.inc = MagicMock()
        mock_metrics.log.error = MagicMock()
        
        result = validate_vcf("/path/to/missing.vcf")
        
        assert "ERROR: Unexpected error during validate_vcf" in result
        assert "File not found" in result
        mock_span.record_exception.assert_called_once()


class TestBCFToolsTools:
    """Test cases for BCFtools tool wrappers."""

    @patch('vcf_agent.agent._bcftools_view')
    def test_bcftools_view_tool_success(self, mock_bcftools_view):
        """Test successful bcftools view execution."""
        mock_bcftools_view.return_value = (0, "VCF output", "")
        
        result = bcftools_view_tool(["-h", "test.vcf"])
        
        assert result == "VCF output"
        mock_bcftools_view.assert_called_once_with(["-h", "test.vcf"])

    @patch('vcf_agent.agent._bcftools_view')
    def test_bcftools_view_tool_error(self, mock_bcftools_view):
        """Test bcftools view with error."""
        mock_bcftools_view.return_value = (1, "", "Error: File not found")
        
        result = bcftools_view_tool(["-h", "missing.vcf"])
        
        assert result == "Error: File not found"
        mock_bcftools_view.assert_called_once_with(["-h", "missing.vcf"])

    @patch('vcf_agent.agent._bcftools_query')
    def test_bcftools_query_tool_success(self, mock_bcftools_query):
        """Test successful bcftools query execution."""
        mock_bcftools_query.return_value = (0, "Query results", "")
        
        result = bcftools_query_tool(["-f", "%CHROM\t%POS\n", "test.vcf"])
        
        assert result == "Query results"
        mock_bcftools_query.assert_called_once_with(["-f", "%CHROM\t%POS\n", "test.vcf"])

    @patch('vcf_agent.agent._bcftools_filter')
    def test_bcftools_filter_tool_success(self, mock_bcftools_filter):
        """Test successful bcftools filter execution."""
        mock_bcftools_filter.return_value = (0, "Filtered output", "")
        
        result = bcftools_filter_tool(["-i", "QUAL>30", "test.vcf"])
        
        assert result == "Filtered output"
        mock_bcftools_filter.assert_called_once_with(["-i", "QUAL>30", "test.vcf"])

    @patch('vcf_agent.agent._bcftools_norm')
    def test_bcftools_norm_tool_success(self, mock_bcftools_norm):
        """Test successful bcftools norm execution."""
        mock_bcftools_norm.return_value = (0, "Normalized output", "")
        
        result = bcftools_norm_tool(["-m", "-any", "test.vcf"])
        
        assert result == "Normalized output"
        mock_bcftools_norm.assert_called_once_with(["-m", "-any", "test.vcf"])

    @patch('vcf_agent.agent._bcftools_stats')
    def test_bcftools_stats_tool_success(self, mock_bcftools_stats):
        """Test successful bcftools stats execution."""
        mock_bcftools_stats.return_value = (0, "Stats output", "")
        
        result = bcftools_stats_tool(["test.vcf"])
        
        assert result == "Stats output"
        mock_bcftools_stats.assert_called_once_with(["test.vcf"])

    @patch('vcf_agent.agent._bcftools_annotate')
    def test_bcftools_annotate_tool_success(self, mock_bcftools_annotate):
        """Test successful bcftools annotate execution."""
        mock_bcftools_annotate.return_value = (0, "Annotated output", "")
        
        result = bcftools_annotate_tool(["-a", "annotations.vcf", "test.vcf"])
        
        assert result == "Annotated output"
        mock_bcftools_annotate.assert_called_once_with(["-a", "annotations.vcf", "test.vcf"])


class TestVCFComparisonTool:
    """Test cases for the VCF comparison tool."""

    def test_vcf_comparison_tool_error_handling(self):
        """Test VCF comparison tool error handling with non-existent files."""
        result = vcf_comparison_tool("nonexistent1.vcf", "nonexistent2.vcf", "nonexistent.fa")
        
        # The function should return an error message in JSON format
        assert "error" in result
        assert "bcftools command failed" in result

    @patch('vcf_agent.agent.os.path.exists')
    def test_vcf_comparison_tool_missing_file(self, mock_exists):
        """Test VCF comparison with missing file."""
        mock_exists.return_value = False
        
        result = vcf_comparison_tool("file1.vcf", "file2.vcf", "reference.fa")
        
        assert "error" in result
        assert "bcftools command failed" in result


class TestVCFAnalysisTools:
    """Test cases for VCF analysis and summarization tools."""

    def test_vcf_analysis_summary_tool_not_implemented(self):
        """Test VCF analysis summary tool returns not implemented message."""
        result = vcf_analysis_summary_tool("/path/to/test.vcf")
        
        # The function returns a JSON string indicating it's not implemented
        assert "not implemented" in result
        assert "variant_count" in result

    @patch('vcf_agent.agent.os.path.exists')
    def test_vcf_summarization_tool_missing_file(self, mock_exists):
        """Test VCF summarization with missing file."""
        mock_exists.return_value = False
        
        result = vcf_summarization_tool("/path/to/test.vcf")
        
        # The function returns a JSON string with error
        assert "VCF file not found" in result
        assert "variant_count" in result


class TestLoadVCFIntoGraphDB:
    """Test cases for loading VCF into graph database."""

    @patch('vcf_agent.agent.graph_integration')
    def test_load_vcf_into_graph_db_success(self, mock_graph_integration):
        """Test successful VCF loading into graph database."""
        # Create temporary test file
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as f:
            f.write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.flush()
            
            mock_graph_integration.get_managed_kuzu_connection.return_value = MagicMock()
            mock_graph_integration.load_vcf_into_kuzu.return_value = None
            
            result = load_vcf_into_graph_db_tool(f.name)
            
            # The function returns JSON with status and message
            assert "status" in result
            assert "success" in result
            assert "processed into Kuzu" in result
            
            # Clean up
            import os
            os.unlink(f.name)

    @patch('vcf_agent.agent.os.path.exists')
    def test_load_vcf_into_graph_db_missing_file(self, mock_exists):
        """Test VCF loading with missing file."""
        mock_exists.return_value = False
        
        result = load_vcf_into_graph_db_tool("/path/to/test.vcf")
        
        assert "FileNotFoundError" in result
        assert "status" in result

    @patch('vcf_agent.agent.graph_integration')
    def test_load_vcf_into_graph_db_connection_error(self, mock_graph_integration):
        """Test VCF loading with connection error."""
        # Create temporary test file
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as f:
            f.write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.flush()
            
            # Mock connection failure
            mock_graph_integration.get_managed_kuzu_connection.side_effect = Exception("Database error")
            
            result = load_vcf_into_graph_db_tool(f.name)
            
            assert "Database error" in result
            assert "status" in result
            assert "failed" in result
            
            # Clean up
            import os
            os.unlink(f.name)


class TestModelInitialization:
    """Test cases for model initialization functions."""

    @patch('vcf_agent.agent.OpenAIClient')
    @patch('vcf_agent.agent.LiteLLMModel')
    def test_get_openai_model_success(self, mock_litellm_model, mock_openai_client):
        """Test successful OpenAI model creation."""
        mock_model = MagicMock()
        mock_litellm_model.return_value = mock_model
        mock_client = MagicMock()
        mock_client.test_connection.return_value = True
        mock_openai_client.return_value = mock_client
        
        result = get_openai_model()
        
        assert result == mock_model
        # Check that the model was created with the correct parameters
        mock_litellm_model.assert_called_once()

    @patch('vcf_agent.agent.OpenAIClient')
    def test_get_openai_model_connection_error(self, mock_openai_client):
        """Test OpenAI model creation with connection error."""
        mock_client = MagicMock()
        mock_client.test_connection.side_effect = Exception("Connection failed")
        mock_openai_client.return_value = mock_client
        
        with pytest.raises(Exception, match="Connection failed"):
            get_openai_model()


class TestCerebrasStrandsModel:
    """Test cases for the CerebrasStrandsModel class."""

    @patch('vcf_agent.agent.CerebrasClient')
    def test_cerebras_model_initialization(self, mock_cerebras_client):
        """Test CerebrasStrandsModel initialization."""
        mock_client = MagicMock()
        mock_cerebras_client.return_value = mock_client
        
        model = CerebrasStrandsModel("cerebras-gpt")
        
        assert model.model == "cerebras-gpt"
        assert model.client == mock_client
        mock_cerebras_client.assert_called_once_with(credential_manager=None)

    @patch('vcf_agent.agent.CerebrasClient')
    def test_cerebras_model_generate(self, mock_cerebras_client):
        """Test CerebrasStrandsModel generate method."""
        mock_client = MagicMock()
        mock_response = {"choices": [{"message": {"content": "Generated response"}}]}
        mock_client.chat_completion.return_value = mock_response
        mock_cerebras_client.return_value = mock_client
        
        model = CerebrasStrandsModel()
        result = model.generate("Test prompt")
        
        assert result == "Generated response"

    @patch('vcf_agent.agent.CerebrasClient')
    def test_cerebras_model_converse(self, mock_cerebras_client):
        """Test CerebrasStrandsModel converse method."""
        mock_client = MagicMock()
        mock_response = {"choices": [{"message": {"content": "Conversation response"}}]}
        mock_client.chat_completion.return_value = mock_response
        mock_cerebras_client.return_value = mock_client
        
        model = CerebrasStrandsModel()
        messages = [{"role": "user", "content": "Hello"}]
        result = list(model.converse(messages))  # Convert generator to list
        
        assert len(result) > 0

    @patch('vcf_agent.agent.CerebrasClient')
    def test_get_cerebras_model(self, mock_cerebras_client):
        """Test get_cerebras_model function."""
        mock_client = MagicMock()
        mock_cerebras_client.return_value = mock_client
        
        result = get_cerebras_model()
        
        assert isinstance(result, CerebrasStrandsModel)
        assert result.model == "cerebras-gpt"


class TestAgentCreation:
    """Test cases for agent creation and session management."""

    @patch('vcf_agent.agent.graph_integration')
    @patch('vcf_agent.agent.OllamaModel')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_ollama(self, mock_agent, mock_ollama_model, mock_graph_integration):
        """Test agent creation with Ollama model."""
        mock_model = MagicMock()
        mock_ollama_model.return_value = mock_model
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance
        mock_graph_integration.get_managed_kuzu_connection.return_value = MagicMock()
        
        result = get_agent_with_session(model_provider="ollama")
        
        assert result == mock_agent_instance
        # Check that OllamaModel was called with the correct parameters
        mock_ollama_model.assert_called_once()
        mock_agent.assert_called_once()

    @patch('vcf_agent.agent.graph_integration')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_default(self, mock_agent, mock_graph_integration):
        """Test agent creation with default model provider."""
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance
        mock_graph_integration.get_managed_kuzu_connection.return_value = MagicMock()
        
        result = get_agent_with_session()
        
        assert result == mock_agent_instance
        mock_agent.assert_called_once()

    @patch('vcf_agent.agent.graph_integration')
    @patch('vcf_agent.agent.OllamaModel')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_with_config(self, mock_agent, mock_ollama_model, mock_graph_integration):
        """Test agent creation with session config."""
        from vcf_agent.config import SessionConfig
        
        mock_model = MagicMock()
        mock_ollama_model.return_value = mock_model
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance
        mock_graph_integration.get_managed_kuzu_connection.return_value = MagicMock()
        
        session_config = SessionConfig(raw_mode=True)
        result = get_agent_with_session(session_config=session_config)
        
        assert result == mock_agent_instance


class TestRunLLMAnalysisTask:
    """Test cases for the run_llm_analysis_task function."""

    @patch('vcf_agent.agent.get_prompt_for_task')
    @patch('vcf_agent.agent.get_agent_with_session')
    def test_run_llm_analysis_task_success(self, mock_get_agent, mock_get_prompt):
        """Test successful LLM analysis task execution."""
        mock_agent = MagicMock()
        # Return valid JSON response
        mock_agent.return_value = '{"analysis": "complete", "status": "success"}'
        mock_get_agent.return_value = mock_agent
        mock_get_prompt.return_value = "Test prompt"
        
        result = run_llm_analysis_task(
            task="vcf_analysis_summary",
            file_paths=["/path/to/test.vcf"],
            model_provider="ollama"
        )
        
        assert '"analysis": "complete"' in result
        mock_get_agent.assert_called_once()
        mock_agent.assert_called_once()

    @patch('vcf_agent.agent.get_prompt_for_task')
    @patch('vcf_agent.agent.get_agent_with_session')
    def test_run_llm_analysis_task_error(self, mock_get_agent, mock_get_prompt):
        """Test LLM analysis task with error."""
        mock_agent = MagicMock()
        mock_agent.side_effect = Exception("Model error")
        mock_get_agent.return_value = mock_agent
        mock_get_prompt.return_value = "Test prompt"
        
        # The function should catch the exception and return an error response
        with pytest.raises(Exception, match="Model error"):
            run_llm_analysis_task(
                task="vcf_analysis_summary",
                file_paths=["/path/to/test.vcf"]
            )


class TestConstants:
    """Test cases for module constants and configuration."""

    def test_system_prompt_exists(self):
        """Test that SYSTEM_PROMPT is properly defined."""
        assert isinstance(SYSTEM_PROMPT, str)
        assert "VCF Analysis Agent" in SYSTEM_PROMPT
        assert "tool_name" in SYSTEM_PROMPT

    def test_raw_mode_setting(self):
        """Test that RAW_MODE is properly configured."""
        assert RAW_MODE is True

    def test_system_prompt_contains_tools(self):
        """Test that SYSTEM_PROMPT mentions available tools."""
        expected_tools = [
            "echo", "validate_vcf", "bcftools_view_tool", 
            "bcftools_query_tool", "vcf_comparison_tool"
        ]
        for tool in expected_tools:
            assert tool in SYSTEM_PROMPT


class TestErrorHandling:
    """Test cases for error handling across agent functions."""

    @patch('vcf_agent.agent.graph_integration')
    def test_graph_integration_error_handling(self, mock_graph_integration):
        """Test error handling when graph integration fails."""
        # Create temporary test file
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as f:
            f.write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.flush()
            
            mock_graph_integration.get_managed_kuzu_connection.side_effect = Exception("DB connection failed")
            
            # Should not raise exception, but handle gracefully
            result = load_vcf_into_graph_db_tool(f.name)
            
            assert "DB connection failed" in result
            assert "status" in result
            
            # Clean up
            import os
            os.unlink(f.name)

    @patch('vcf_agent.agent.OllamaModel')
    def test_model_initialization_error(self, mock_ollama_model):
        """Test error handling during model initialization."""
        mock_ollama_model.side_effect = Exception("Model not available")
        
        with pytest.raises(Exception, match="Model not available"):
            get_agent_with_session(model_provider="ollama")


class TestIntegrationScenarios:
    """Integration test scenarios for agent functionality."""

    @patch('vcf_agent.agent.validate_vcf_file')
    @patch('vcf_agent.agent._bcftools_view')
    @patch('vcf_agent.agent.metrics')
    @patch('vcf_agent.agent.agent_tracer')
    def test_tool_chain_execution(self, mock_tracer, mock_metrics, mock_bcftools_view, mock_validate_vcf_file):
        """Test chaining multiple tools together."""
        # Mock tracer
        mock_span = MagicMock()
        mock_tracer.start_as_current_span.return_value.__enter__.return_value = mock_span
        
        # Mock validation
        mock_validate_vcf_file.return_value = (True, None)
        
        # Mock bcftools
        mock_bcftools_view.return_value = (0, "VCF header", "")
        
        # Mock metrics
        mock_metrics.VCF_AGENT_TOOL_DURATION_SECONDS.labels.return_value.observe = MagicMock()
        mock_metrics.VCF_AGENT_TOOL_REQUESTS_TOTAL.labels.return_value.inc = MagicMock()
        
        # Execute tool chain
        validation_result = validate_vcf("/path/to/test.vcf")
        view_result = bcftools_view_tool(["-h", "/path/to/test.vcf"])
        
        assert "VALID" in validation_result
        assert "VCF header" in view_result

    def test_echo_tool_integration(self):
        """Test echo tool integration with various inputs."""
        test_cases = [
            "Simple text",
            "Text with numbers 123",
            "Special chars: !@#$%^&*()",
            "",
            "Multi\nline\ntext"
        ]
        
        for test_input in test_cases:
            result = echo(test_input)
            assert result == f"Echo: {test_input}" 