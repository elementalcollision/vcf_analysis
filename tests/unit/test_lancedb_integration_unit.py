"""
Unit tests for LanceDB integration module (src/vcf_agent/lancedb_integration.py),
focusing on Kuzu enrichment aspects.
"""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock, call
import re

# Import the module to be tested
from vcf_agent import lancedb_integration
from vcf_agent import graph_integration # For mocking
from vcf_agent.lancedb_integration import mask_sensitive_sql

class TestLanceDBIntegrationUnit:
    """Unit tests for lancedb_integration.py Kuzu enrichment."""

    def test_example_placeholder_lancedb(self):
        """A placeholder test for lancedb_integration Kuzu parts."""
        assert True

    @patch('vcf_agent.lancedb_integration.graph_integration')
    # No need to patch get_lancedb_table, we pass a mock table directly
    def test_search_by_embedding_successful_enrichment(self, mock_graph_integration_module):
        """Test search_by_embedding with successful Kuzu enrichment."""
        mock_table = MagicMock() # This is the LanceDB table object

        # Sample LanceDB results
        lancedb_results_data = {
            'variant_id': ['var1', 'var2'],
            'score': [0.9, 0.85],
        }
        mock_lancedb_df = pd.DataFrame(lancedb_results_data)
        # Mock the .to_pandas() method call
        mock_q_object = mock_table.search.return_value.limit.return_value
        mock_q_object.to_pandas.return_value = mock_lancedb_df

        # Mock Kuzu connection and context
        mock_kuzu_conn = MagicMock()
        mock_graph_integration_module.get_managed_kuzu_connection.return_value = mock_kuzu_conn
        kuzu_context_results = [
            {'variant_id': 'var1', 'sample_id': 'sampleA', 'zygosity': 'HOM'},
            {'variant_id': 'var1', 'sample_id': 'sampleB', 'zygosity': 'HET'},
            {'variant_id': 'var2', 'sample_id': 'sampleA', 'zygosity': 'HET'},
        ]
        mock_graph_integration_module.get_variant_context.return_value = kuzu_context_results

        # Call the function being tested
        result_df = lancedb_integration.search_by_embedding(
            table=mock_table, 
            query_embedding=[0.1]*1024, 
            limit=5
        )

        # Assertions for LanceDB part (interaction with the provided table mock)
        mock_table.search.assert_called_once_with([0.1]*1024, vector_column_name='embedding')
        mock_table.search.return_value.limit.assert_called_once_with(5)

        # Assertions for Kuzu part
        mock_graph_integration_module.get_managed_kuzu_connection.assert_called_once()
        mock_graph_integration_module.get_variant_context.assert_called_once_with(
            mock_kuzu_conn, 
            ['var1', 'var2']
        )

        # Assertions for DataFrame enrichment
        assert 'kuzu_observed_samples' in result_df.columns
        expected_kuzu_col_data = [
            [('sampleA', 'HOM'), ('sampleB', 'HET')], # For var1
            [('sampleA', 'HET')]                      # For var2
        ]
        
        # Compare the contents of the kuzu_observed_samples column
        # Ensuring lists of tuples are sorted for consistent comparison if necessary.
        # The function itself sorts the tuples by sample_id, then zygosity.
        for i, expected_list in enumerate(expected_kuzu_col_data):
            actual_list = result_df['kuzu_observed_samples'].iloc[i]
            assert sorted(actual_list) == sorted(expected_list)
        
        assert result_df['variant_id'].tolist() == ['var1', 'var2']
        assert result_df['score'].tolist() == [0.9, 0.85]

    @patch('vcf_agent.lancedb_integration.graph_integration')
    def test_search_by_embedding_no_lancedb_results(self, mock_graph_integration_module):
        """Test search_by_embedding when LanceDB returns an empty DataFrame."""
        mock_table = MagicMock()

        mock_empty_df = pd.DataFrame() # Explicitly create an empty DataFrame
        mock_q_object = mock_table.search.return_value.limit.return_value
        mock_q_object.to_pandas.return_value = mock_empty_df # Return the actual empty DataFrame

        result_df = lancedb_integration.search_by_embedding(
            table=mock_table,
            query_embedding=[0.2]*1024,
            limit=3
        )

        mock_table.search.assert_called_once_with([0.2]*1024, vector_column_name='embedding')
        mock_table.search.return_value.limit.assert_called_once_with(3)

        mock_graph_integration_module.get_managed_kuzu_connection.assert_not_called()
        mock_graph_integration_module.get_variant_context.assert_not_called()

        assert result_df.empty
        # search_by_embedding adds this column even to an empty DataFrame
        assert 'kuzu_observed_samples' in result_df.columns

    @patch('vcf_agent.lancedb_integration.logger')
    @patch('vcf_agent.lancedb_integration.graph_integration')
    def test_search_by_embedding_lancedb_found_kuzu_empty_context(self, mock_graph_integration_module, mock_logger):
        """Test when LanceDB finds variants, but Kuzu returns no context for them."""
        mock_table = MagicMock()

        lancedb_results_data = {'variant_id': ['var3', 'var4'], 'score': [0.7, 0.6]}
        mock_lancedb_df = pd.DataFrame(lancedb_results_data)
        mock_q_object = mock_table.search.return_value.limit.return_value
        mock_q_object.to_pandas.return_value = mock_lancedb_df

        mock_kuzu_conn = MagicMock()
        mock_graph_integration_module.get_managed_kuzu_connection.return_value = mock_kuzu_conn
        # Kuzu returns an empty list, meaning no context found
        mock_graph_integration_module.get_variant_context.return_value = [] 

        result_df = lancedb_integration.search_by_embedding(
            table=mock_table,
            query_embedding=[0.3]*1024,
            limit=2
        )

        mock_table.search.assert_called_once_with([0.3]*1024, vector_column_name='embedding')
        mock_table.search.return_value.limit.assert_called_once_with(2)

        mock_graph_integration_module.get_managed_kuzu_connection.assert_called_once()
        mock_graph_integration_module.get_variant_context.assert_called_once_with(
            mock_kuzu_conn, 
            ['var3', 'var4']
        )

        mock_logger.warning.assert_not_called()

        assert 'kuzu_observed_samples' in result_df.columns
        # Expect lists of empty lists, or lists of pd.NA, or just empty lists if no variants had context
        # The current implementation in search_by_embedding will result in empty lists for variants with no context
        for item in result_df['kuzu_observed_samples']:
            assert item == [] # Each variant row should have an empty list for kuzu_observed_samples
        
        assert result_df['variant_id'].tolist() == ['var3', 'var4']

    @patch('vcf_agent.lancedb_integration.logger')
    @patch('vcf_agent.lancedb_integration.graph_integration')
    def test_search_by_embedding_kuzu_partial_context(self, mock_graph_integration_module, mock_logger):
        """Test when Kuzu returns context for a subset of LanceDB variants."""
        mock_table = MagicMock()

        lancedb_results_data = {'variant_id': ['var5', 'var6', 'var7'], 'score': [0.5, 0.4, 0.3]}
        mock_lancedb_df = pd.DataFrame(lancedb_results_data)
        mock_q_object = mock_table.search.return_value.limit.return_value
        mock_q_object.to_pandas.return_value = mock_lancedb_df

        mock_kuzu_conn = MagicMock()
        mock_graph_integration_module.get_managed_kuzu_connection.return_value = mock_kuzu_conn
        # Kuzu returns context only for 'var5' and 'var7'
        kuzu_partial_context = [
            {'variant_id': 'var5', 'sample_id': 'sampleC', 'zygosity': 'HOM'},
            {'variant_id': 'var7', 'sample_id': 'sampleD', 'zygosity': 'HET'},
        ]
        mock_graph_integration_module.get_variant_context.return_value = kuzu_partial_context

        result_df = lancedb_integration.search_by_embedding(
            table=mock_table,
            query_embedding=[0.4]*1024,
            limit=3
        )

        mock_table.search.assert_called_once_with([0.4]*1024, vector_column_name='embedding')
        mock_table.search.return_value.limit.assert_called_once_with(3)

        mock_graph_integration_module.get_variant_context.assert_called_once_with(
            mock_kuzu_conn, 
            ['var5', 'var6', 'var7']
        )

        mock_logger.warning.assert_not_called()

        assert 'kuzu_observed_samples' in result_df.columns
        expected_kuzu_data = [
            [('sampleC', 'HOM')], # Context for var5
            [],                   # No context for var6
            [('sampleD', 'HET')]  # Context for var7
        ]
        
        actual_kuzu_data = result_df['kuzu_observed_samples'].tolist()
        for i, expected_list in enumerate(expected_kuzu_data):
            assert sorted(actual_kuzu_data[i]) == sorted(expected_list)
            
        assert result_df['variant_id'].tolist() == ['var5', 'var6', 'var7']

    @patch('vcf_agent.lancedb_integration.logger') # To check log messages
    @patch('vcf_agent.lancedb_integration.graph_integration')
    def test_search_by_embedding_kuzu_error(self, mock_graph_integration_module, mock_logger):
        """Test when a Kuzu call (e.g., get_variant_context) raises an exception."""
        mock_table = MagicMock()

        lancedb_results_data = {'variant_id': ['var8'], 'score': [0.2]}
        mock_lancedb_df = pd.DataFrame(lancedb_results_data)
        mock_q_object = mock_table.search.return_value.limit.return_value
        mock_q_object.to_pandas.return_value = mock_lancedb_df

        mock_kuzu_conn = MagicMock()
        mock_graph_integration_module.get_managed_kuzu_connection.return_value = mock_kuzu_conn
        # Configure get_variant_context to raise an exception
        mock_graph_integration_module.get_variant_context.side_effect = Exception("Kuzu Cconnection Error")

        result_df = lancedb_integration.search_by_embedding(
            table=mock_table,
            query_embedding=[0.5]*1024,
            limit=1
        )

        mock_table.search.assert_called_once_with([0.5]*1024, vector_column_name='embedding')
        mock_table.search.return_value.limit.assert_called_once_with(1)

        mock_graph_integration_module.get_managed_kuzu_connection.assert_called_once()
        mock_graph_integration_module.get_variant_context.assert_called_once_with(mock_kuzu_conn, ['var8'])

        # Check that an error was logged
        mock_logger.error.assert_called_once()
        assert "Error during Kuzu context enrichment" in mock_logger.error.call_args[0][0]
        assert "Kuzu Cconnection Error" in mock_logger.error.call_args[0][0]

        # DataFrame should still contain LanceDB results
        assert not result_df.empty
        assert result_df['variant_id'].tolist() == ['var8']
        # The 'kuzu_observed_samples' column should exist and be filled with pd.NA or empty lists
        # The function initializes it with pd.NA, and on error, it should remain as such or be empty lists.
        # Current implementation leaves it as pd.NA from initialization if enrichment fails.
        assert 'kuzu_observed_samples' in result_df.columns
        # Check that all values in kuzu_observed_samples are NA or empty lists
        # Based on current code: it's initialized with pd.NA
        # And if context_map is empty or error, it applies this (which result in pd.NA for non-matches)
        # So, expect pd.NA here. Pandas .isna() is good for this.
        assert result_df['kuzu_observed_samples'].isna().all() 

    # TODO: Add tests for agent.py (Kuzu-related tools and initialization)
    #   - Successful enrichment: Kuzu returns context for found variants
    #   - No variants found in LanceDB (Kuzu part should not be called)
    #   - Variants found in LanceDB, but Kuzu returns no context for them
    #   - Kuzu returns context for a subset of variants
    #   - Error during Kuzu connection/query (e.g., graph_integration.get_variant_context raises an exception) 

class TestMaskSensitiveSQL:
    """Unit tests for mask_sensitive_sql utility function."""

    def test_mask_email_and_ssn(self):
        s = "email = 'john.doe@example.com' AND ssn = '123-45-6789'"
        masked = mask_sensitive_sql(s)
        assert "[MASKED]" in masked
        assert masked.count("[MASKED]") >= 2

    def test_mask_long_number(self):
        s = "user_id = 12345678901"
        masked = mask_sensitive_sql(s)
        assert "[MASKED]" in masked or masked == '[MASKED_SQL]'

    def test_mask_name_heuristic(self):
        s = "name = 'John Smith'"
        masked = mask_sensitive_sql(s)
        assert "[MASKED]" in masked or masked == '[MASKED_SQL]'

    def test_mask_sensitive_column_value(self):
        s = "patient_id = 'abc123' AND genotype = 'A/B'"
        masked = mask_sensitive_sql(s)
        assert masked.count("[MASKED]") >= 2

    def test_unparseable_returns_masked_sql(self):
        s = None
        masked = mask_sensitive_sql(s)
        assert masked == '[MASKED_SQL]'
        s = 12345
        masked = mask_sensitive_sql(s)
        assert masked == '[MASKED_SQL]'

    def test_no_masking_for_non_sensitive(self):
        s = "chrom = '1' AND pos > 1000"
        masked = mask_sensitive_sql(s)
        assert masked == s

    def test_custom_sensitive_columns_and_patterns(self):
        s = "foo = 'bar' AND secret = 'baz'"
        masked = mask_sensitive_sql(s, sensitive_columns=['secret'], patterns=[re.compile(r'baz')])
        assert '[MASKED]' in masked or masked == '[MASKED_SQL]' 