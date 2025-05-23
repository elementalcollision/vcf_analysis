"""Tests for API clients."""

import os
import json
import tempfile
from unittest import mock
from typing import Dict, Any

import pytest

from vcf_agent.api_clients import (
    CredentialManager, 
    OpenAIClient, 
    CerebrasClient, 
    APIClientError
)


class TestCredentialManager:
    """Tests for CredentialManager."""
    
    def test_load_from_env(self):
        """Test loading credentials from environment variables."""
        with mock.patch.dict(os.environ, {
            "OPENAI_API_KEY": "test-openai-key",
            "CEREBRAS_API_KEY": "test-cerebras-key"
        }):
            manager = CredentialManager()
            
            assert manager.get_credential("openai") == "test-openai-key"
            assert manager.get_credential("cerebras") == "test-cerebras-key"
    
    def test_load_from_file(self):
        """Test loading credentials from file."""
        creds = {
            "openai": {"api_key": "file-openai-key"},
            "cerebras": {"api_key": "file-cerebras-key"}
        }
        
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            json.dump(creds, f)
            creds_file = f.name
            
        try:
            manager = CredentialManager(credentials_file=creds_file)
            
            assert manager.get_credential("openai") == "file-openai-key"
            assert manager.get_credential("cerebras") == "file-cerebras-key"
        finally:
            os.unlink(creds_file)
    
    def test_file_overrides_env(self):
        """Test that file credentials override environment variables."""
        creds = {
            "openai": {"api_key": "file-openai-key"},
            "cerebras": {"api_key": "file-cerebras-key"}
        }
        
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            json.dump(creds, f)
            creds_file = f.name
            
        try:
            with mock.patch.dict(os.environ, {
                "OPENAI_API_KEY": "env-openai-key",
                "CEREBRAS_API_KEY": "env-cerebras-key"
            }):
                manager = CredentialManager(credentials_file=creds_file)
                
                assert manager.get_credential("openai") == "file-openai-key"
                assert manager.get_credential("cerebras") == "file-cerebras-key"
        finally:
            os.unlink(creds_file)
    
    def test_missing_credential(self):
        """Test that missing credential raises error."""
        manager = CredentialManager()
        
        with pytest.raises(APIClientError, match="No credentials found for service"):
            manager.get_credential("nonexistent")


class TestOpenAIClient:
    """Tests for OpenAIClient."""
    
    def test_initialization(self):
        """Test client initialization."""
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.return_value = "test-key"
        
        client = OpenAIClient(credential_manager=mock_cred_manager)
        
        mock_cred_manager.get_credential.assert_called_once_with("openai")
        assert client.client is not None
    
    def test_missing_credential(self):
        """Test that missing credential raises error."""
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.side_effect = APIClientError("Test error")
        
        with pytest.raises(APIClientError, match="Test error"):
            OpenAIClient(credential_manager=mock_cred_manager)
    
    @mock.patch("openai.OpenAI")
    def test_test_connection(self, mock_openai):
        """Test connection testing."""
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.return_value = "test-key"
        
        mock_client = mock.MagicMock()
        mock_openai.return_value = mock_client
        
        client = OpenAIClient(credential_manager=mock_cred_manager)
        result = client.test_connection()
        
        assert result is True
        mock_client.models.list.assert_called_once()
        
    @mock.patch("openai.OpenAI")
    def test_chat_completion(self, mock_openai):
        """Test chat completion."""
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.return_value = "test-key"
        
        mock_client = mock.MagicMock()
        mock_openai.return_value = mock_client
        
        mock_response = {"id": "test-id", "choices": [{"message": {"content": "Test response"}}]}
        mock_client.chat.completions.create.return_value = mock_response
        
        client = OpenAIClient(credential_manager=mock_cred_manager)
        messages = [{"role": "user", "content": "Test message"}]
        response = client.chat_completion(messages, model="test-model")
        
        assert response == mock_response
        mock_client.chat.completions.create.assert_called_once_with(
            model="test-model",
            messages=messages,
            temperature=0.7,
            max_tokens=None
        )


class TestCerebrasClient:
    """Tests for CerebrasClient."""
    
    def test_initialization(self):
        """Test client initialization."""
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.return_value = "test-key"
        
        client = CerebrasClient(credential_manager=mock_cred_manager)
        
        mock_cred_manager.get_credential.assert_called_once_with("cerebras")
        assert client.api_key == "test-key"
    
    def test_missing_credential(self):
        """Test that missing credential raises error."""
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.side_effect = APIClientError("Test error")
        
        with pytest.raises(APIClientError, match="Test error"):
            CerebrasClient(credential_manager=mock_cred_manager)
    
    def test_test_connection(self):
        """Test connection testing."""
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.return_value = "test-key"
        
        client = CerebrasClient(credential_manager=mock_cred_manager)
        result = client.test_connection()
        
        assert result is True
    
    def test_chat_completion(self):
        """Test chat completion (mock implementation)."""
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.return_value = "test-key"
        
        client = CerebrasClient(credential_manager=mock_cred_manager)
        messages = [{"role": "user", "content": "Test message"}]
        response = client.chat_completion(messages, model="test-model")
        
        assert "id" in response
        assert "choices" in response
        assert response["model"] == "test-model"
        assert "This is a mock response from Cerebras test-model" in response["choices"][0]["message"]["content"] 