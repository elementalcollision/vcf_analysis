"""
Unit tests for agent model initialization and factory functions.

Tests the model creation functions, error handling, and agent factory methods
in the agent.py module.
"""

import pytest
import os
from unittest.mock import patch, MagicMock, Mock
from vcf_agent.agent import (
    get_openai_model, 
    get_cerebras_model, 
    CerebrasStrandsModel,
    get_agent_with_session
)
from vcf_agent.config import SessionConfig
from vcf_agent.api_clients import APIClientError


class TestGetOpenAIModel:
    """Test cases for the get_openai_model function."""

    @patch('vcf_agent.agent.OpenAIClient')
    @patch('vcf_agent.agent.LiteLLMModel')
    def test_get_openai_model_success(self, mock_litellm_model, mock_openai_client):
        """Test successful OpenAI model creation."""
        # Setup mocks
        mock_client_instance = MagicMock()
        mock_client_instance.test_connection.return_value = True
        mock_client_instance.credential_manager.get_credential.return_value = "test-api-key"
        mock_openai_client.return_value = mock_client_instance
        
        mock_model_instance = MagicMock()
        mock_litellm_model.return_value = mock_model_instance

        # Execute
        result = get_openai_model()

        # Verify
        assert result == mock_model_instance
        mock_openai_client.assert_called_once_with(credential_manager=None)
        mock_client_instance.test_connection.assert_called_once()
        mock_client_instance.credential_manager.get_credential.assert_called_once_with("openai")
        mock_litellm_model.assert_called_once_with(
            client_args={"api_key": "test-api-key"},
            model_id="openai/gpt-4o",
            params={
                "temperature": 0.7,
                "max_tokens": 4000
            }
        )

    @patch('vcf_agent.agent.OpenAIClient')
    @patch('vcf_agent.agent.LiteLLMModel')
    def test_get_openai_model_with_credential_manager(self, mock_litellm_model, mock_openai_client):
        """Test OpenAI model creation with custom credential manager."""
        # Setup mocks
        mock_credential_manager = MagicMock()
        mock_client_instance = MagicMock()
        mock_client_instance.test_connection.return_value = True
        mock_client_instance.credential_manager.get_credential.return_value = "custom-api-key"
        mock_openai_client.return_value = mock_client_instance
        
        mock_model_instance = MagicMock()
        mock_litellm_model.return_value = mock_model_instance

        # Execute
        result = get_openai_model(credential_manager=mock_credential_manager)

        # Verify
        assert result == mock_model_instance
        mock_openai_client.assert_called_once_with(credential_manager=mock_credential_manager)
        mock_client_instance.test_connection.assert_called_once()
        mock_client_instance.credential_manager.get_credential.assert_called_once_with("openai")
        mock_litellm_model.assert_called_once_with(
            client_args={"api_key": "custom-api-key"},
            model_id="openai/gpt-4o",
            params={
                "temperature": 0.7,
                "max_tokens": 4000
            }
        )

    @patch('vcf_agent.agent.OpenAIClient')
    @patch('vcf_agent.agent.logging')
    def test_get_openai_model_api_client_error(self, mock_logging, mock_openai_client):
        """Test OpenAI model creation with API client error."""
        # Setup mock to raise APIClientError
        mock_openai_client.side_effect = APIClientError("API connection failed")

        # Execute and verify exception
        with pytest.raises(APIClientError, match="API connection failed"):
            get_openai_model()

        # Verify logging
        mock_logging.error.assert_called_once_with("Failed to initialize OpenAI model: API connection failed")

    @patch('vcf_agent.agent.OpenAIClient')
    @patch('vcf_agent.agent.logging')
    def test_get_openai_model_test_connection_error(self, mock_logging, mock_openai_client):
        """Test OpenAI model creation when test_connection fails."""
        # Setup mocks
        mock_client_instance = MagicMock()
        mock_client_instance.test_connection.side_effect = APIClientError("Connection test failed")
        mock_openai_client.return_value = mock_client_instance

        # Execute and verify exception
        with pytest.raises(APIClientError, match="Connection test failed"):
            get_openai_model()

        # Verify logging
        mock_logging.error.assert_called_once_with("Failed to initialize OpenAI model: Connection test failed")


class TestCerebrasStrandsModel:
    """Test cases for the CerebrasStrandsModel class."""

    @patch('vcf_agent.agent.CerebrasClient')
    def test_cerebras_strands_model_init_success(self, mock_cerebras_client):
        """Test successful CerebrasStrandsModel initialization."""
        # Setup mock
        mock_client_instance = MagicMock()
        mock_client_instance.test_connection.return_value = True
        mock_cerebras_client.return_value = mock_client_instance

        # Execute
        model = CerebrasStrandsModel(model="test-model")

        # Verify
        assert model.model == "test-model"
        assert model.client == mock_client_instance
        mock_cerebras_client.assert_called_once_with(credential_manager=None)
        mock_client_instance.test_connection.assert_called_once()

    @patch('vcf_agent.agent.CerebrasClient')
    def test_cerebras_strands_model_init_with_credential_manager(self, mock_cerebras_client):
        """Test CerebrasStrandsModel initialization with custom credential manager."""
        # Setup mock
        mock_credential_manager = MagicMock()
        mock_client_instance = MagicMock()
        mock_client_instance.test_connection.return_value = True
        mock_cerebras_client.return_value = mock_client_instance

        # Execute
        model = CerebrasStrandsModel(credential_manager=mock_credential_manager)

        # Verify
        assert model.model == "cerebras-gpt"  # default value
        assert model.client == mock_client_instance
        mock_cerebras_client.assert_called_once_with(credential_manager=mock_credential_manager)
        mock_client_instance.test_connection.assert_called_once()

    @patch('vcf_agent.agent.CerebrasClient')
    def test_cerebras_strands_model_generate(self, mock_cerebras_client):
        """Test CerebrasStrandsModel generate method."""
        # Setup mock
        mock_client_instance = MagicMock()
        mock_client_instance.test_connection.return_value = True
        mock_client_instance.chat_completion.return_value = {
            "choices": [{"message": {"content": "Generated response"}}]
        }
        mock_cerebras_client.return_value = mock_client_instance

        # Execute
        model = CerebrasStrandsModel()
        result = model.generate("Test prompt", temperature=0.5)

        # Verify
        assert result == "Generated response"
        mock_client_instance.chat_completion.assert_called_once_with(
            messages=[{"role": "user", "content": "Test prompt"}],
            model="cerebras-gpt",
            temperature=0.5
        )

    @patch('vcf_agent.agent.CerebrasClient')
    def test_cerebras_strands_model_converse(self, mock_cerebras_client):
        """Test CerebrasStrandsModel converse method."""
        # Setup mock
        mock_client_instance = MagicMock()
        mock_client_instance.test_connection.return_value = True
        mock_cerebras_client.return_value = mock_client_instance

        # Execute
        model = CerebrasStrandsModel()
        messages = [{"role": "user", "content": "Test message"}]
        result = list(model.converse(messages))

        # Verify
        assert len(result) == 1
        response = result[0]
        assert response["role"] == "assistant"
        assert response["content"] == "This is a mock Cerebras response."
        assert response["tool_calls"] == []
        assert response["usage"]["prompt_tokens"] == 1
        assert response["usage"]["completion_tokens"] == 1
        assert response["usage"]["total_tokens"] == 2


class TestGetCerebrasModel:
    """Test cases for the get_cerebras_model function."""

    @patch('vcf_agent.agent.CerebrasStrandsModel')
    def test_get_cerebras_model_success(self, mock_cerebras_strands_model):
        """Test successful Cerebras model creation."""
        # Setup mock
        mock_model_instance = MagicMock()
        mock_cerebras_strands_model.return_value = mock_model_instance

        # Execute
        result = get_cerebras_model()

        # Verify
        assert result == mock_model_instance
        mock_cerebras_strands_model.assert_called_once_with(credential_manager=None)

    @patch('vcf_agent.agent.CerebrasStrandsModel')
    def test_get_cerebras_model_with_credential_manager(self, mock_cerebras_strands_model):
        """Test Cerebras model creation with custom credential manager."""
        # Setup mock
        mock_credential_manager = MagicMock()
        mock_model_instance = MagicMock()
        mock_cerebras_strands_model.return_value = mock_model_instance

        # Execute
        result = get_cerebras_model(credential_manager=mock_credential_manager)

        # Verify
        assert result == mock_model_instance
        mock_cerebras_strands_model.assert_called_once_with(credential_manager=mock_credential_manager)

    @patch('vcf_agent.agent.CerebrasStrandsModel')
    @patch('vcf_agent.agent.logging')
    def test_get_cerebras_model_api_client_error(self, mock_logging, mock_cerebras_strands_model):
        """Test Cerebras model creation with API client error."""
        # Setup mock to raise APIClientError
        mock_cerebras_strands_model.side_effect = APIClientError("Cerebras API connection failed")

        # Execute and verify exception
        with pytest.raises(APIClientError, match="Cerebras API connection failed"):
            get_cerebras_model()

        # Verify logging
        mock_logging.error.assert_called_once_with("Failed to initialize Cerebras model: Cerebras API connection failed")


class TestGetAgentWithSession:
    """Test cases for the get_agent_with_session function."""

    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    @patch('vcf_agent.agent.OllamaModel')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_ollama_default(self, mock_agent, mock_ollama_model, mock_kuzu_connection):
        """Test agent creation with default Ollama model."""
        # Setup mocks
        mock_kuzu_connection.return_value = MagicMock()
        mock_model_instance = MagicMock()
        mock_ollama_model.return_value = mock_model_instance
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance

        # Execute
        result = get_agent_with_session()

        # Verify
        assert result == mock_agent_instance
        mock_kuzu_connection.assert_called_once()
        mock_ollama_model.assert_called_once_with(
            host="http://127.0.0.1:11434",
            model_id="qwen3:4b"
        )
        mock_agent.assert_called_once()

    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    @patch('vcf_agent.agent.LiteLLMModel')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_openai(self, mock_agent, mock_litellm_model, mock_kuzu_connection):
        """Test agent creation with OpenAI model."""
        # Setup mocks
        mock_kuzu_connection.return_value = MagicMock()
        mock_model_instance = MagicMock()
        mock_litellm_model.return_value = mock_model_instance
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance

        # Execute
        result = get_agent_with_session(model_provider="openai")

        # Verify
        assert result == mock_agent_instance
        mock_kuzu_connection.assert_called_once()
        mock_litellm_model.assert_called_once_with(model_id="gpt-4o")
        mock_agent.assert_called_once()

    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    @patch('vcf_agent.agent.LiteLLMModel')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_cerebras(self, mock_agent, mock_litellm_model, mock_kuzu_connection):
        """Test agent creation with Cerebras model."""
        # Setup mocks
        mock_kuzu_connection.return_value = MagicMock()
        mock_model_instance = MagicMock()
        mock_litellm_model.return_value = mock_model_instance
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance

        # Execute
        result = get_agent_with_session(model_provider="cerebras")

        # Verify
        assert result == mock_agent_instance
        mock_kuzu_connection.assert_called_once()
        mock_litellm_model.assert_called_once_with(model_id="azure/cerebras-gpt-13b")
        mock_agent.assert_called_once()

    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    def test_get_agent_with_session_unsupported_provider(self, mock_kuzu_connection):
        """Test agent creation with unsupported model provider."""
        # Setup mock
        mock_kuzu_connection.return_value = MagicMock()

        # Execute and verify exception
        with pytest.raises(ValueError, match="Unsupported model provider: invalid_provider"):
            get_agent_with_session(model_provider="invalid_provider")

    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    @patch('vcf_agent.agent.OllamaModel')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_custom_session_config(self, mock_agent, mock_ollama_model, mock_kuzu_connection):
        """Test agent creation with custom session config."""
        # Setup mocks
        mock_kuzu_connection.return_value = MagicMock()
        mock_model_instance = MagicMock()
        mock_ollama_model.return_value = mock_model_instance
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance

        # Create custom session config
        session_config = SessionConfig(ollama_model_name="custom-model", raw_mode=True)

        # Execute
        result = get_agent_with_session(session_config=session_config)

        # Verify
        assert result == mock_agent_instance
        mock_kuzu_connection.assert_called_once()
        mock_ollama_model.assert_called_once_with(
            host="http://127.0.0.1:11434",
            model_id="custom-model"
        )
        mock_agent.assert_called_once()

    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    @patch('vcf_agent.agent.OllamaModel')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_kuzu_error(self, mock_agent, mock_ollama_model, mock_kuzu_connection):
        """Test agent creation when Kuzu initialization fails."""
        # Setup mocks
        mock_kuzu_connection.side_effect = RuntimeError("Kuzu connection failed")
        mock_model_instance = MagicMock()
        mock_ollama_model.return_value = mock_model_instance
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance

        # Execute (should not raise exception, just log error)
        result = get_agent_with_session()

        # Verify agent is still created despite Kuzu error
        assert result == mock_agent_instance
        mock_kuzu_connection.assert_called_once()
        mock_ollama_model.assert_called_once()
        mock_agent.assert_called_once()

    @patch.dict(os.environ, {'VCF_AGENT_RAW_MODE': 'true'})
    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    @patch('vcf_agent.agent.OllamaModel')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_raw_mode_env_var(self, mock_agent, mock_ollama_model, mock_kuzu_connection):
        """Test agent creation with raw mode set via environment variable."""
        # Setup mocks
        mock_kuzu_connection.return_value = MagicMock()
        mock_model_instance = MagicMock()
        mock_ollama_model.return_value = mock_model_instance
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance

        # Execute
        result = get_agent_with_session()

        # Verify agent is created (raw mode logic is internal)
        assert result == mock_agent_instance
        mock_kuzu_connection.assert_called_once()
        mock_ollama_model.assert_called_once()
        mock_agent.assert_called_once()

    @patch.dict(os.environ, {'OLLAMA_HOST': 'http://custom-host:11434'})
    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    @patch('vcf_agent.agent.OllamaModel')
    @patch('vcf_agent.agent.Agent')
    def test_get_agent_with_session_custom_ollama_host(self, mock_agent, mock_ollama_model, mock_kuzu_connection):
        """Test agent creation with custom Ollama host from environment."""
        # Setup mocks
        mock_kuzu_connection.return_value = MagicMock()
        mock_model_instance = MagicMock()
        mock_ollama_model.return_value = mock_model_instance
        mock_agent_instance = MagicMock()
        mock_agent.return_value = mock_agent_instance

        # Execute
        result = get_agent_with_session()

        # Verify
        assert result == mock_agent_instance
        mock_kuzu_connection.assert_called_once()
        mock_ollama_model.assert_called_once_with(
            host="http://custom-host:11434",
            model_id="qwen3:4b"
        )
        mock_agent.assert_called_once() 