"""
Unit tests for the Kuzu graph integration module (src/vcf_agent/graph_integration.py).
"""

import pytest
from unittest.mock import patch, MagicMock, call, ANY
import pandas as pd
import kuzu

# Import the module to be tested
from vcf_agent import graph_integration

# Since DEFAULT_KUZU_DB_PATH is in graph_integration, we can use it or mock it too.
# For unit tests, we often want to avoid actual file system operations.

# Fixture to manage the global _kuzu_main_connection state
@pytest.fixture
def managed_kuzu_conn_fixture():
    # Patch the correct singleton variable for the test
    with patch.object(graph_integration, '_kuzu_main_connection', None, create=True):
        yield

class TestGraphIntegrationUnit:
    """Unit tests for graph_integration.py"""

    def test_example_placeholder(self):
        """A placeholder test to ensure the file is picked up by pytest."""
        assert True

    # TODO: Add unit tests for get_managed_kuzu_connection
    # TODO: Add unit tests for create_schema
    # TODO: Add unit tests for add_variant
    # TODO: Add unit tests for add_sample
    # TODO: Add unit tests for link_variant_to_sample
    # TODO: Add unit tests for execute_query
    # TODO: Add unit tests for get_variant_context

    @patch('vcf_agent.graph_integration.kuzu.Connection')
    @patch('vcf_agent.graph_integration.kuzu.Database')
    def test_get_kuzu_db_connection_default_path(self, mock_kuzu_database_constructor, mock_kuzu_connection_constructor):
        """Test get_kuzu_db_connection with the default database path."""
        mock_db_instance = MagicMock(name="MockKuzuDB")
        mock_conn_instance = MagicMock(name="MockKuzuConn")
        
        mock_kuzu_database_constructor.return_value = mock_db_instance
        mock_kuzu_connection_constructor.return_value = mock_conn_instance

        default_path = graph_integration.DEFAULT_KUZU_DB_PATH
        conn = graph_integration.get_kuzu_db_connection()

        mock_kuzu_database_constructor.assert_called_once_with(default_path)
        mock_kuzu_connection_constructor.assert_called_once_with(mock_db_instance)
        assert conn == mock_conn_instance

    @patch('vcf_agent.graph_integration.kuzu.Connection')
    @patch('vcf_agent.graph_integration.kuzu.Database')
    def test_get_kuzu_db_connection_custom_path(self, mock_kuzu_database_constructor, mock_kuzu_connection_constructor):
        """Test get_kuzu_db_connection with a custom database path."""
        mock_db_instance = MagicMock(name="MockKuzuDB")
        mock_conn_instance = MagicMock(name="MockKuzuConn")
        mock_kuzu_database_constructor.return_value = mock_db_instance
        mock_kuzu_connection_constructor.return_value = mock_conn_instance
        
        custom_path = "/tmp/test_custom.kuzu"
        conn = graph_integration.get_kuzu_db_connection(db_path=custom_path)

        mock_kuzu_database_constructor.assert_called_once_with(custom_path)
        mock_kuzu_connection_constructor.assert_called_once_with(mock_db_instance)
        assert conn == mock_conn_instance

    @patch('vcf_agent.graph_integration.kuzu.Connection')
    @patch('vcf_agent.graph_integration.kuzu.Database')
    def test_get_kuzu_db_connection_with_readonly(self, mock_kuzu_database_constructor, mock_kuzu_connection_constructor):
        """Test get_kuzu_db_connection with read_only flag (though unused in current kuzu calls)."""
        mock_db_instance = MagicMock(name="MockKuzuDB")
        mock_conn_instance = MagicMock(name="MockKuzuConn")
        mock_kuzu_database_constructor.return_value = mock_db_instance
        mock_kuzu_connection_constructor.return_value = mock_conn_instance

        default_path = graph_integration.DEFAULT_KUZU_DB_PATH
        conn = graph_integration.get_kuzu_db_connection(read_only=True)

        mock_kuzu_database_constructor.assert_called_once_with(default_path)
        mock_kuzu_connection_constructor.assert_called_once_with(mock_db_instance)
        assert conn == mock_conn_instance

    def test_create_schema(self):
        """Test the create_schema function to ensure it executes the correct Cypher queries."""
        mock_conn = MagicMock()
        
        graph_integration.create_schema(mock_conn)
        
        expected_calls = [
            call("CREATE NODE TABLE IF NOT EXISTS Sample(sample_id STRING, PRIMARY KEY (sample_id))"),
            call("CREATE NODE TABLE IF NOT EXISTS Variant(variant_id STRING, chrom STRING, pos INT64, ref STRING, alt STRING, rs_id STRING, PRIMARY KEY (variant_id))"),
            call("CREATE REL TABLE IF NOT EXISTS ObservedIn(FROM Variant TO Sample, zygosity STRING)")
        ]
        
        mock_conn.execute.assert_has_calls(expected_calls, any_order=True)
        assert mock_conn.execute.call_count == len(expected_calls)

    @patch('vcf_agent.graph_integration.create_schema')
    @patch('vcf_agent.graph_integration.get_kuzu_db_connection')
    def test_get_managed_kuzu_connection_first_call(self, mock_get_kuzu_db_conn, mock_create_schema, managed_kuzu_conn_fixture):
        """Test get_managed_kuzu_connection when no instance exists."""
        mock_new_conn = MagicMock()
        mock_get_kuzu_db_conn.return_value = mock_new_conn

        returned_conn = graph_integration.get_managed_kuzu_connection()

        mock_get_kuzu_db_conn.assert_called_once_with(db_path=graph_integration.DEFAULT_KUZU_DB_PATH)
        mock_create_schema.assert_called_once_with(mock_new_conn)
        assert returned_conn == mock_new_conn
        assert graph_integration._kuzu_main_connection == mock_new_conn

    @patch('vcf_agent.graph_integration.create_schema')
    @patch('vcf_agent.graph_integration.get_kuzu_db_connection')
    def test_get_managed_kuzu_connection_subsequent_calls(self, mock_get_kuzu_db_conn, mock_create_schema, managed_kuzu_conn_fixture):
        """Test get_managed_kuzu_connection when an instance already exists."""
        # First call to establish the connection
        mock_first_conn = MagicMock()
        mock_get_kuzu_db_conn.return_value = mock_first_conn

        first_returned_conn = graph_integration.get_managed_kuzu_connection()
        assert first_returned_conn == mock_first_conn
        mock_get_kuzu_db_conn.assert_called_once_with(db_path=graph_integration.DEFAULT_KUZU_DB_PATH)
        mock_create_schema.assert_called_once_with(mock_first_conn)

        # Reset mocks for the second call, but _kuzu_conn_instance should now be mock_first_conn
        mock_get_kuzu_db_conn.reset_mock()
        mock_create_schema.reset_mock()

        # Second call
        second_returned_conn = graph_integration.get_managed_kuzu_connection()

        mock_get_kuzu_db_conn.assert_not_called() # Should not create a new DB connection
        mock_create_schema.assert_not_called()    # Should not re-create schema
        assert second_returned_conn == first_returned_conn # Should return the existing connection
        assert graph_integration._kuzu_main_connection == first_returned_conn

    def test_add_variant(self):
        """Test the add_variant function."""
        mock_conn = MagicMock()
        variant_data = {
            "variant_id": "test_variant_1",
            "chrom": "chr1",
            "pos": 100,
            "ref": "A",
            "alt": "T",
            "rs_id": "rsTest1"
        }

        graph_integration.add_variant(mock_conn, variant_data)

        expected_query = """
        CREATE (v:Variant {
            variant_id: $variant_id,
            chrom: $chrom,
            pos: $pos,
            ref: $ref,
            alt: $alt,
            rs_id: $rs_id
        })
        """
        mock_conn.execute.assert_called_once_with(expected_query, parameters=variant_data)

    def test_add_sample(self):
        """Test the add_sample function."""
        mock_conn = MagicMock()
        sample_data = {"sample_id": "test_sample_1"}

        graph_integration.add_sample(mock_conn, sample_data)

        expected_query = "CREATE (s:Sample {sample_id: $sample_id})"
        mock_conn.execute.assert_called_once_with(expected_query, parameters=sample_data)

    def test_link_variant_to_sample(self):
        """Test the link_variant_to_sample function."""
        mock_conn = MagicMock()
        sample_id = "test_sample_1"
        variant_id = "test_variant_1"
        properties = {"zygosity": "HOM"}

        graph_integration.link_variant_to_sample(mock_conn, sample_id, variant_id, properties)

        expected_query = """
        MATCH (v:Variant {variant_id: $v_id}), (s:Sample {sample_id: $s_id})
        CREATE (v)-[r:ObservedIn {zygosity: $zygosity}]->(s)
        """
        expected_params_dict = {
            "v_id": variant_id,
            "s_id": sample_id,
            "zygosity": properties['zygosity']
        }
        mock_conn.execute.assert_called_once_with(expected_query, parameters=expected_params_dict)

    def test_execute_query_with_params(self):
        """Test execute_query with parameters."""
        mock_conn = MagicMock()
        mock_query_result = MagicMock()
        # mock_df should be a real DataFrame instance for isinstance check
        expected_result_list = [{"colA": 1, "colB": "val1"}, {"colA": 2, "colB": "val2"}]
        mock_df_instance = pd.DataFrame(expected_result_list) # Create a real DataFrame

        mock_conn.execute.return_value = mock_query_result
        mock_query_result.get_as_df.return_value = mock_df_instance # Return the real DataFrame
        # mock_df.to_dict.return_value = expected_result_list # No longer needed if mock_df_instance is used directly

        test_query = "MATCH (n) WHERE n.prop = $param RETURN n.prop AS colA, n.other AS colB"
        test_params = {"param": "test_value"}

        result = graph_integration.execute_query(mock_conn, test_query, test_params)

        mock_conn.execute.assert_called_once_with(test_query, parameters=test_params)
        mock_query_result.get_as_df.assert_called_once()
        # The function itself calls .to_dict(orient='records') on the DataFrame
        assert result == expected_result_list

    def test_execute_query_no_params(self):
        """Test execute_query without parameters."""
        mock_conn = MagicMock()
        mock_query_result = MagicMock()
        # mock_df should be a real DataFrame instance
        expected_result_list = [{"colA": 3, "colB": "val3"}]
        mock_df_instance = pd.DataFrame(expected_result_list) # Create a real DataFrame

        mock_conn.execute.return_value = mock_query_result
        mock_query_result.get_as_df.return_value = mock_df_instance # Return the real DataFrame
        # mock_df.to_dict.return_value = expected_result_list

        test_query = "MATCH (n) RETURN n.prop AS colA, n.other AS colB"

        result = graph_integration.execute_query(mock_conn, test_query)

        mock_conn.execute.assert_called_once_with(test_query, parameters={})
        mock_query_result.get_as_df.assert_called_once()
        # The function itself calls .to_dict(orient='records') on the DataFrame
        assert result == expected_result_list

    @patch('vcf_agent.graph_integration.execute_query')
    def test_get_variant_context(self, mock_execute_query):
        """Test the get_variant_context function."""
        mock_conn = MagicMock() 
        variant_ids = ["var1", "var2"]
        
        # This is what execute_query would return (a flat list)
        mock_execute_query_return_value = [
            {"variant_id": "var1", "sample_id": "sampleA", "zygosity": "HOM"},
            {"variant_id": "var1", "sample_id": "sampleB", "zygosity": "HET"}, # var1 observed in two samples
            {"variant_id": "var2", "sample_id": "sampleA", "zygosity": "HET"}  # var2 observed in one sample
        ]
        mock_execute_query.return_value = mock_execute_query_return_value

        # This is the expected final structure from get_variant_context
        expected_context_map = {
            "var1": [
                {"sample_id": "sampleA", "zygosity": "HOM"},
                {"sample_id": "sampleB", "zygosity": "HET"}
            ],
            "var2": [
                {"sample_id": "sampleA", "zygosity": "HET"}
            ]
        }

        result_map = graph_integration.get_variant_context(mock_conn, variant_ids)

        expected_query_text = """
    UNWIND $v_ids AS target_variant_id
    MATCH (v:Variant {variant_id: target_variant_id})-[o:ObservedIn]->(s:Sample)
    RETURN v.variant_id AS variant_id, s.sample_id AS sample_id, o.zygosity AS zygosity
    """
        expected_params = {"v_ids": variant_ids}
        
        mock_execute_query.assert_called_once_with(mock_conn, expected_query_text, params=expected_params)
        assert result_map == expected_context_map
