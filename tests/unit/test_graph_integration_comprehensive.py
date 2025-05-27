"""
Comprehensive unit tests for graph integration module.

Tests the Kuzu database integration functions including connection management,
schema creation, data insertion, and query execution.
"""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock, Mock
from vcf_agent.graph_integration import (
    get_managed_kuzu_connection,
    get_kuzu_db_connection,
    create_schema,
    add_variant,
    add_sample,
    link_variant_to_sample,
    execute_query,
    get_variant_context,
    close_kuzu_connection,
    DEFAULT_KUZU_DB_PATH
)

# Create a mock class that will pass isinstance checks
class MockKuzuQueryResult:
    """Mock class to simulate kuzu.query_result.QueryResult"""
    def get_as_df(self):
        return pd.DataFrame()


class TestGetManagedKuzuConnection:
    """Test cases for the get_managed_kuzu_connection function."""

    @patch('vcf_agent.graph_integration.get_kuzu_db_connection')
    @patch('vcf_agent.graph_integration.create_schema')
    def test_get_managed_kuzu_connection_success(self, mock_create_schema, mock_get_kuzu_db_connection):
        """Test successful managed Kuzu connection creation."""
        # Setup mocks
        mock_connection = MagicMock()
        mock_get_kuzu_db_connection.return_value = mock_connection
        
        # Reset global connection to None for clean test
        import vcf_agent.graph_integration
        vcf_agent.graph_integration._kuzu_main_connection = None

        # Execute
        result = get_managed_kuzu_connection()

        # Verify
        assert result == mock_connection
        mock_get_kuzu_db_connection.assert_called_once_with(db_path=DEFAULT_KUZU_DB_PATH)
        mock_create_schema.assert_called_once_with(mock_connection)

    @patch('vcf_agent.graph_integration.get_kuzu_db_connection')
    def test_get_managed_kuzu_connection_existing(self, mock_get_kuzu_db_connection):
        """Test that existing connection is returned without recreating."""
        # Setup existing connection
        existing_connection = MagicMock()
        import vcf_agent.graph_integration
        vcf_agent.graph_integration._kuzu_main_connection = existing_connection

        # Execute
        result = get_managed_kuzu_connection()

        # Verify
        assert result == existing_connection
        mock_get_kuzu_db_connection.assert_not_called()

    @patch('vcf_agent.graph_integration.get_kuzu_db_connection')
    @patch('vcf_agent.graph_integration.create_schema')
    def test_get_managed_kuzu_connection_error(self, mock_create_schema, mock_get_kuzu_db_connection):
        """Test error handling in managed Kuzu connection creation."""
        # Setup mocks to raise exception
        mock_get_kuzu_db_connection.side_effect = Exception("Database connection failed")
        
        # Reset global connection to None for clean test
        import vcf_agent.graph_integration
        vcf_agent.graph_integration._kuzu_main_connection = None

        # Execute and verify exception
        with pytest.raises(RuntimeError, match="Kuzu setup failed: Database connection failed"):
            get_managed_kuzu_connection()


class TestGetKuzuDbConnection:
    """Test cases for the get_kuzu_db_connection function."""

    @patch('vcf_agent.graph_integration.kuzu.Database')
    @patch('vcf_agent.graph_integration.kuzu.Connection')
    def test_get_kuzu_db_connection_success(self, mock_connection_class, mock_database_class):
        """Test successful Kuzu database connection."""
        # Setup mocks
        mock_database = MagicMock()
        mock_database_class.return_value = mock_database
        mock_connection = MagicMock()
        mock_connection_class.return_value = mock_connection

        # Execute
        result = get_kuzu_db_connection("/test/path")

        # Verify
        assert result == mock_connection
        mock_database_class.assert_called_once_with("/test/path")
        mock_connection_class.assert_called_once_with(mock_database)

    @patch('vcf_agent.graph_integration.kuzu.Database')
    def test_get_kuzu_db_connection_error(self, mock_database_class):
        """Test error handling in Kuzu database connection."""
        # Setup mock to raise exception
        mock_database_class.side_effect = Exception("Database initialization failed")

        # Execute and verify exception
        with pytest.raises(Exception, match="Database initialization failed"):
            get_kuzu_db_connection("/test/path")


class TestCreateSchema:
    """Test cases for the create_schema function."""

    def test_create_schema_success(self):
        """Test successful schema creation."""
        # Setup mock connection
        mock_connection = MagicMock()
        mock_query_result = MagicMock()
        mock_connection.execute.return_value = mock_query_result

        # Execute
        create_schema(mock_connection)

        # Verify all schema queries were executed
        assert mock_connection.execute.call_count == 3
        expected_queries = [
            "CREATE NODE TABLE IF NOT EXISTS Variant(variant_id STRING, chrom STRING, pos INT64, ref STRING, alt STRING, rs_id STRING, PRIMARY KEY (variant_id))",
            "CREATE NODE TABLE IF NOT EXISTS Sample(sample_id STRING, PRIMARY KEY (sample_id))",
            "CREATE REL TABLE IF NOT EXISTS ObservedIn(FROM Variant TO Sample, zygosity STRING)"
        ]
        for i, expected_query in enumerate(expected_queries):
            actual_call = mock_connection.execute.call_args_list[i]
            assert actual_call[0][0] == expected_query

    def test_create_schema_error(self):
        """Test error handling in schema creation."""
        # Setup mock connection to raise exception
        mock_connection = MagicMock()
        mock_connection.execute.side_effect = Exception("Schema creation failed")

        # Execute and verify exception
        with pytest.raises(Exception, match="Schema creation failed"):
            create_schema(mock_connection)


class TestAddVariant:
    """Test cases for the add_variant function."""

    def test_add_variant_success(self):
        """Test successful variant addition."""
        # Setup mock connection
        mock_connection = MagicMock()
        mock_query_result = MagicMock()
        mock_connection.execute.return_value = mock_query_result

        # Test data
        variant_data = {
            'variant_id': 'chr1-123-A-G',
            'chrom': '1',
            'pos': 123,
            'ref': 'A',
            'alt': 'G',
            'rs_id': 'rs12345'
        }

        # Execute
        add_variant(mock_connection, variant_data)

        # Verify
        mock_connection.execute.assert_called_once()
        call_args = mock_connection.execute.call_args
        assert "CREATE (v:Variant" in call_args[0][0]
        assert call_args[1]['parameters']['variant_id'] == 'chr1-123-A-G'
        assert call_args[1]['parameters']['chrom'] == '1'
        assert call_args[1]['parameters']['pos'] == 123

    def test_add_variant_without_rs_id(self):
        """Test variant addition without rs_id."""
        # Setup mock connection
        mock_connection = MagicMock()
        mock_query_result = MagicMock()
        mock_connection.execute.return_value = mock_query_result

        # Test data without rs_id
        variant_data = {
            'variant_id': 'chr1-123-A-G',
            'chrom': '1',
            'pos': 123,
            'ref': 'A',
            'alt': 'G'
        }

        # Execute
        add_variant(mock_connection, variant_data)

        # Verify rs_id defaults to empty string
        call_args = mock_connection.execute.call_args
        assert call_args[1]['parameters']['rs_id'] == ''

    def test_add_variant_error(self):
        """Test error handling in variant addition."""
        # Setup mock connection to raise exception
        mock_connection = MagicMock()
        mock_connection.execute.side_effect = Exception("Variant addition failed")

        variant_data = {
            'variant_id': 'chr1-123-A-G',
            'chrom': '1',
            'pos': 123,
            'ref': 'A',
            'alt': 'G'
        }

        # Execute and verify exception
        with pytest.raises(Exception, match="Variant addition failed"):
            add_variant(mock_connection, variant_data)


class TestAddSample:
    """Test cases for the add_sample function."""

    def test_add_sample_success(self):
        """Test successful sample addition."""
        # Setup mock connection
        mock_connection = MagicMock()
        mock_query_result = MagicMock()
        mock_connection.execute.return_value = mock_query_result

        # Test data
        sample_data = {'sample_id': 'SAMPLE_001'}

        # Execute
        add_sample(mock_connection, sample_data)

        # Verify
        mock_connection.execute.assert_called_once()
        call_args = mock_connection.execute.call_args
        assert "CREATE (s:Sample" in call_args[0][0]
        assert call_args[1]['parameters']['sample_id'] == 'SAMPLE_001'

    def test_add_sample_error(self):
        """Test error handling in sample addition."""
        # Setup mock connection to raise exception
        mock_connection = MagicMock()
        mock_connection.execute.side_effect = Exception("Sample addition failed")

        sample_data = {'sample_id': 'SAMPLE_001'}

        # Execute and verify exception
        with pytest.raises(Exception, match="Sample addition failed"):
            add_sample(mock_connection, sample_data)


class TestLinkVariantToSample:
    """Test cases for the link_variant_to_sample function."""

    def test_link_variant_to_sample_success(self):
        """Test successful variant-sample linking."""
        # Setup mock connection
        mock_connection = MagicMock()
        mock_query_result = MagicMock()
        mock_connection.execute.return_value = mock_query_result

        # Execute
        link_variant_to_sample(mock_connection, 'SAMPLE_001', 'chr1-123-A-G', {'zygosity': 'HET'})

        # Verify
        mock_connection.execute.assert_called_once()
        call_args = mock_connection.execute.call_args
        assert "MATCH (v:Variant" in call_args[0][0]
        assert "CREATE (v)-[r:ObservedIn" in call_args[0][0]
        assert call_args[1]['parameters']['s_id'] == 'SAMPLE_001'
        assert call_args[1]['parameters']['v_id'] == 'chr1-123-A-G'
        assert call_args[1]['parameters']['zygosity'] == 'HET'

    def test_link_variant_to_sample_missing_zygosity(self):
        """Test error when zygosity is missing."""
        mock_connection = MagicMock()

        # Execute and verify exception
        with pytest.raises(ValueError, match="Missing 'zygosity' in properties"):
            link_variant_to_sample(mock_connection, 'SAMPLE_001', 'chr1-123-A-G', {})

    def test_link_variant_to_sample_error(self):
        """Test error handling in variant-sample linking."""
        # Setup mock connection to raise exception
        mock_connection = MagicMock()
        mock_connection.execute.side_effect = Exception("Linking failed")

        # Execute and verify exception
        with pytest.raises(Exception, match="Linking failed"):
            link_variant_to_sample(mock_connection, 'SAMPLE_001', 'chr1-123-A-G', {'zygosity': 'HET'})


class TestExecuteQuery:
    """Test cases for the execute_query function."""

    @patch('vcf_agent.graph_integration.kuzu')
    def test_execute_query_success(self, mock_kuzu):
        """Test successful query execution."""
        # Setup mock connection and query result
        mock_connection = MagicMock()
        mock_query_result = MagicMock()
        mock_df = pd.DataFrame([{'variant_id': 'chr1-123-A-G', 'sample_id': 'SAMPLE_001'}])
        mock_query_result.get_as_df.return_value = mock_df
        
        # Configure the mock to pass isinstance checks
        mock_kuzu.query_result.QueryResult = type('QueryResult', (), {})
        mock_query_result.__class__ = mock_kuzu.query_result.QueryResult
        mock_connection.execute.return_value = mock_query_result

        # Execute
        result = execute_query(mock_connection, "MATCH (v:Variant) RETURN v.variant_id", {'param': 'value'})

        # Verify
        assert len(result) == 1
        assert result[0]['variant_id'] == 'chr1-123-A-G'
        assert result[0]['sample_id'] == 'SAMPLE_001'
        mock_connection.execute.assert_called_once_with("MATCH (v:Variant) RETURN v.variant_id", parameters={'param': 'value'})

    @patch('vcf_agent.graph_integration.kuzu')
    def test_execute_query_no_params(self, mock_kuzu):
        """Test query execution without parameters."""
        # Setup mock connection and query result
        mock_connection = MagicMock()
        mock_query_result = MagicMock()
        mock_df = pd.DataFrame([])
        mock_query_result.get_as_df.return_value = mock_df
        
        # Configure the mock to pass isinstance checks
        mock_kuzu.query_result.QueryResult = type('QueryResult', (), {})
        mock_query_result.__class__ = mock_kuzu.query_result.QueryResult
        mock_connection.execute.return_value = mock_query_result

        # Execute
        result = execute_query(mock_connection, "MATCH (v:Variant) RETURN v.variant_id")

        # Verify
        assert result == []
        mock_connection.execute.assert_called_once_with("MATCH (v:Variant) RETURN v.variant_id", parameters={})

    @patch('vcf_agent.graph_integration.kuzu')
    def test_execute_query_list_result(self, mock_kuzu):
        """Test query execution when result is a list."""
        # Setup mock connection and query result
        mock_connection = MagicMock()
        mock_query_result = MagicMock()
        mock_df = pd.DataFrame([{'variant_id': 'chr1-123-A-G'}])
        mock_query_result.get_as_df.return_value = mock_df
        
        # Configure the mock to pass isinstance checks
        mock_kuzu.query_result.QueryResult = type('QueryResult', (), {})
        mock_query_result.__class__ = mock_kuzu.query_result.QueryResult
        mock_connection.execute.return_value = [mock_query_result]

        # Execute
        result = execute_query(mock_connection, "MATCH (v:Variant) RETURN v.variant_id")

        # Verify
        assert len(result) == 1
        assert result[0]['variant_id'] == 'chr1-123-A-G'

    def test_execute_query_empty_list_result(self):
        """Test query execution when result is an empty list."""
        # Setup mock connection
        mock_connection = MagicMock()
        mock_connection.execute.return_value = []

        # Execute
        result = execute_query(mock_connection, "MATCH (v:Variant) RETURN v.variant_id")

        # Verify
        assert result == []

    def test_execute_query_none_result(self):
        """Test query execution when result is None."""
        # Setup mock connection
        mock_connection = MagicMock()
        mock_connection.execute.return_value = None

        # Execute
        result = execute_query(mock_connection, "MATCH (v:Variant) RETURN v.variant_id")

        # Verify
        assert result == []

    def test_execute_query_error(self):
        """Test error handling in query execution."""
        # Setup mock connection to raise exception
        mock_connection = MagicMock()
        mock_connection.execute.side_effect = Exception("Query execution failed")

        # Execute and verify exception
        with pytest.raises(Exception, match="Query execution failed"):
            execute_query(mock_connection, "INVALID QUERY")


class TestGetVariantContext:
    """Test cases for the get_variant_context function."""

    @patch('vcf_agent.graph_integration.kuzu')
    @patch('vcf_agent.graph_integration.get_managed_kuzu_connection')
    def test_get_variant_context_success(self, mock_get_managed_kuzu_connection, mock_kuzu):
        """Test successful variant context retrieval."""
        # Setup mock connection and query result
        mock_connection = MagicMock()
        mock_get_managed_kuzu_connection.return_value = mock_connection
        mock_query_result = MagicMock()
        mock_df = pd.DataFrame([
            {'variant_id': 'chr1-123-A-G', 'sample_id': 'SAMPLE_001', 'zygosity': 'HET'},
            {'variant_id': 'chr1-456-C-T', 'sample_id': 'SAMPLE_002', 'zygosity': 'HOM'}
        ])
        mock_query_result.get_as_df.return_value = mock_df
        
        # Configure the mock to pass isinstance checks
        mock_kuzu.query_result.QueryResult = type('QueryResult', (), {})
        mock_query_result.__class__ = mock_kuzu.query_result.QueryResult
        mock_connection.execute.return_value = mock_query_result

        # Execute
        result = get_variant_context(None, ['chr1-123-A-G', 'chr1-456-C-T'])

        # Verify
        assert 'chr1-123-A-G' in result
        assert 'chr1-456-C-T' in result
        assert len(result['chr1-123-A-G']) == 1
        assert result['chr1-123-A-G'][0]['sample_id'] == 'SAMPLE_001'
        assert result['chr1-123-A-G'][0]['zygosity'] == 'HET'

    @patch('vcf_agent.graph_integration.kuzu')
    def test_get_variant_context_with_connection(self, mock_kuzu):
        """Test variant context retrieval with provided connection."""
        # Setup mock connection and query result
        mock_connection = MagicMock()
        mock_query_result = MagicMock()
        mock_df = pd.DataFrame([])
        mock_query_result.get_as_df.return_value = mock_df
        
        # Configure the mock to pass isinstance checks
        mock_kuzu.query_result.QueryResult = type('QueryResult', (), {})
        mock_query_result.__class__ = mock_kuzu.query_result.QueryResult
        mock_connection.execute.return_value = mock_query_result

        # Execute
        result = get_variant_context(mock_connection, ['chr1-123-A-G'])

        # Verify
        assert 'chr1-123-A-G' in result
        assert result['chr1-123-A-G'] == []

    def test_get_variant_context_empty_list(self):
        """Test variant context retrieval with empty variant list."""
        # Execute
        result = get_variant_context(None, [])

        # Verify
        assert result == {}

    @patch('vcf_agent.graph_integration.get_managed_kuzu_connection')
    def test_get_variant_context_no_connection(self, mock_get_managed_kuzu_connection):
        """Test error when managed connection is not available."""
        # Setup mock to return None
        mock_get_managed_kuzu_connection.return_value = None

        # Execute and verify exception
        with pytest.raises(RuntimeError, match="Managed Kuzu connection could not be established"):
            get_variant_context(None, ['chr1-123-A-G'])


class TestCloseKuzuConnection:
    """Test cases for the close_kuzu_connection function."""

    @patch('vcf_agent.graph_integration.get_kuzu_db_connection')
    def test_close_kuzu_connection_success(self, mock_get_kuzu_db_connection):
        """Test successful Kuzu connection closure."""
        # Setup mock connection
        mock_connection = MagicMock()
        mock_get_kuzu_db_connection.return_value = mock_connection

        # Execute
        close_kuzu_connection("/test/path")

        # Verify
        mock_get_kuzu_db_connection.assert_called_once_with("/test/path")
        mock_connection.close.assert_called_once()

    @patch('vcf_agent.graph_integration.get_kuzu_db_connection')
    def test_close_kuzu_connection_error(self, mock_get_kuzu_db_connection):
        """Test error handling in connection closure."""
        # Setup mock to raise exception
        mock_get_kuzu_db_connection.side_effect = Exception("Connection closure failed")

        # Execute and verify exception
        with pytest.raises(Exception, match="Connection closure failed"):
            close_kuzu_connection("/test/path") 