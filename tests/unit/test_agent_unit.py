"""
Unit tests for the Agent (src/vcf_agent/agent.py),
focusing on Kuzu-related tools and initialization.
"""

import pytest
import json
import os # Added for os.environ patch
from unittest.mock import patch, MagicMock, ANY

# Import the module to be tested
from vcf_agent import agent
# Removed direct tool imports as they are not used in this specific test for initialization

@pytest.fixture
def mock_kuzu_init():
    """Mocks graph_integration.get_managed_kuzu_connection."""
    with patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection') as mock_get_conn:
        mock_get_conn.return_value = MagicMock() # Ensure it returns a mock connection
        yield mock_get_conn

class TestAgentKuzuInitialization:
    """Tests for Kuzu initialization points within the agent."""

    def test_get_agent_with_session_initializes_kuzu(self, mock_kuzu_init):
        """Test that get_agent_with_session calls get_managed_kuzu_connection."""
        
        # Mock environment variables required by get_agent_with_session
        # Ensure default values are reasonable if not all are Kuzu-related.
        with patch.dict(os.environ, {
            'LANCEDB_TABLE_NAME': 'test_table', 
            'LANCEDB_DB_PATH': '/tmp/test_db',
            'OPENAI_API_KEY': 'test_key', # Added as it might be needed by agent instantiation
            'KUZU_DB_PATH': '/tmp/test_kuzu.db' # Ensure Kuzu path is also available if used by get_managed_kuzu_connection internally
        }, clear=True):
            # test_session_id = "test_session_123" # No longer used
            # Call the function that initializes the agent and Kuzu connection
            vcf_agent_instance = agent.get_agent_with_session() # REMOVED session_id
            
            # Assert that Kuzu connection manager was called
            # The specific path argument to get_managed_kuzu_connection depends on its internal logic (e.g. env var or default)
            # For this test, we only care that it *was* called.
            mock_kuzu_init.assert_called_once()
            assert vcf_agent_instance is not None
            # Removed assertion for vcf_agent_instance.session_id as it's not an attribute

    # TODO: Add tests for load_vcf_into_graph_db_tool
    #   - Successful VCF load
    #   - FileNotFoundError for VCF
    #   - Kuzu connection not available (simulated error from get_managed_kuzu_connection)

class TestLoadVcfTool:
    """Tests for the load_vcf_into_graph_db_tool."""

    @pytest.fixture
    def agent_instance(self, mock_kuzu_init): # mock_kuzu_init fixture from TestAgentKuzuInitialization will run
        """Provides a VCFAgent instance with mocked Kuzu initialization for tool testing."""
        with patch.dict(os.environ, {
            'LANCEDB_TABLE_NAME': 'test_table', 
            'LANCEDB_DB_PATH': '/tmp/test_db',
            'OPENAI_API_KEY': 'test_key',
            'KUZU_DB_PATH': '/tmp/test_kuzu.db'
        }, clear=True):
            # The mock_kuzu_init fixture (autouse=False here by not being in class scope directly)
            # ensures get_managed_kuzu_connection is already mocked when get_agent_with_session is called.
            # For the tool test, we need finer control over get_managed_kuzu_connection for some scenarios.
            return agent.get_agent_with_session() # REMOVED session_id

    @patch('vcf_agent.agent.vcf_utils.populate_kuzu_from_vcf')
    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    def test_load_vcf_tool_success(self, mock_get_kuzu_conn, mock_populate_vcf, agent_instance):
        """Test successful VCF load via the agent tool."""
        mock_kuzu_conn_instance = MagicMock()
        mock_get_kuzu_conn.return_value = mock_kuzu_conn_instance
        mock_populate_vcf.return_value = {"variants": 10, "samples": 1, "links": 10}
        
        filepath = "/test/dummy.vcf"
        sample_override = "sample_001"
        
        # Call using the function name
        try:
            result_json = agent_instance.load_vcf_into_graph_db_tool(filepath=filepath, sample_name_override=sample_override)
        except AttributeError:
            # Fallback: call via .tools dict
            result_json = agent_instance.tools['load_vcf_into_graph_db_tool'](filepath=filepath, sample_name_override=sample_override)
        result = json.loads(result_json)

        mock_get_kuzu_conn.assert_called_once()
        mock_populate_vcf.assert_called_once_with(mock_kuzu_conn_instance, filepath, sample_override)
        assert result["status"] == "success"
        assert result["counts"] == {"variants": 10, "samples": 1, "links": 10}

    @patch('vcf_agent.agent.vcf_utils.populate_kuzu_from_vcf')
    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    def test_load_vcf_tool_file_not_found(self, mock_get_kuzu_conn, mock_populate_vcf, agent_instance):
        """Test VCF load tool when VCF file is not found."""
        mock_kuzu_conn_instance = MagicMock()
        mock_get_kuzu_conn.return_value = mock_kuzu_conn_instance
        mock_populate_vcf.side_effect = FileNotFoundError("VCF not found")
        
        filepath = "/test/nonexistent.vcf"
        try:
            result_json = agent_instance.load_vcf_into_graph_db_tool(filepath=filepath) # Call by function name
        except AttributeError:
            result_json = agent_instance.tools['load_vcf_into_graph_db_tool'](filepath=filepath)
        result = json.loads(result_json)

        mock_get_kuzu_conn.assert_called_once()
        mock_populate_vcf.assert_called_once_with(mock_kuzu_conn_instance, filepath, None)
        assert result["status"] == "failed"
        assert "FileNotFoundError" in result["error"]

    @patch('vcf_agent.agent.vcf_utils.populate_kuzu_from_vcf') # Mock to ensure it's not called
    @patch('vcf_agent.agent.graph_integration.get_managed_kuzu_connection')
    def test_load_vcf_tool_kuzu_conn_unavailable(self, mock_get_kuzu_conn, mock_populate_vcf, agent_instance):
        """Test VCF load tool when Kuzu connection is unavailable."""
        mock_get_kuzu_conn.return_value = None # Simulate Kuzu connection failure
        
        filepath = "/test/any.vcf"
        try:
            result_json = agent_instance.load_vcf_into_graph_db_tool(filepath=filepath) # Call by function name
        except AttributeError:
            result_json = agent_instance.tools['load_vcf_into_graph_db_tool'](filepath=filepath)
        result = json.loads(result_json)

        mock_get_kuzu_conn.assert_called_once()
        mock_populate_vcf.assert_not_called()
        assert result["status"] == "failed"
        assert "Kuzu connection not available" in result["error"]
