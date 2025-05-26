"""Tests for API clients."""

import os
import json
import tempfile
from unittest import mock
from typing import Dict, Any
import unittest

import pytest

from vcf_agent.api_clients import (
    CredentialManager, 
    OpenAIClient, 
    CerebrasClient, 
    APIClientError
)


class TestCredentialManager(unittest.TestCase):
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


class TestOpenAIClient(unittest.TestCase):
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
    
    @mock.patch("vcf_agent.api_clients.OpenAI")
    def test_test_connection(self, mock_openai_constructor):
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.return_value = "test-api-key"

        mock_sdk_instance = mock.MagicMock()
        mock_openai_constructor.return_value = mock_sdk_instance
        
        mock_sdk_instance.models.list.return_value = {"data": [{"id": "model-1"}]}

        client = OpenAIClient(credential_manager=mock_cred_manager)
        self.assertTrue(client.test_connection())
        
        mock_openai_constructor.assert_called_once_with(api_key="test-api-key")
        mock_sdk_instance.models.list.assert_called_once()

    @mock.patch("vcf_agent.api_clients.OpenAI")
    def test_chat_completion(self, mock_openai_constructor):
        mock_cred_manager = mock.MagicMock()
        mock_cred_manager.get_credential.return_value = "test-api-key"
    
        mock_sdk_instance = mock.MagicMock()
        mock_openai_constructor.return_value = mock_sdk_instance
    
        mock_response_dict = {"id": "test-id", "choices": [{"message": {"content": "Test response"}}]}
        mock_sdk_instance.chat.completions.create.return_value = mock_response_dict
    
        client = OpenAIClient(credential_manager=mock_cred_manager)
        messages = [{"role": "user", "content": "Test message"}]
        response = client.chat_completion(messages, model="test-model")

        self.assertEqual(response["choices"][0]["message"]["content"], "Test response")
        mock_openai_constructor.assert_called_once_with(api_key="test-api-key")
        mock_sdk_instance.chat.completions.create.assert_called_once_with(
            model="test-model",
            messages=messages,
            temperature=0.7, 
            max_tokens=None,
        )


class TestCerebrasClient(unittest.TestCase):
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

class TestOllamaClient(unittest.TestCase):
    pass 