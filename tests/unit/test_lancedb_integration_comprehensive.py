"""
Comprehensive unit tests for LanceDB integration module.

Tests the LanceDB integration functions including database connection,
table management, variant operations, search functionality, and security features.
"""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock, Mock, call
from vcf_agent.lancedb_integration import (
    mask_sensitive_sql,
    Variant,
    get_db,
    get_or_create_table,
    add_variants,
    search_by_embedding,
    update_variant,
    delete_variants,
    create_scalar_index,
    filter_variants_by_metadata,
    SENSITIVE_COLUMNS,
    SENSITIVE_PATTERNS,
    lancedb_write_lock
)


class TestMaskSensitiveSQL:
    """Test cases for the mask_sensitive_sql function."""

    def test_mask_sensitive_sql_empty_input(self):
        """Test masking with empty or None input."""
        assert mask_sensitive_sql("") == '[MASKED_SQL]'
        assert mask_sensitive_sql(None) == '[MASKED_SQL]'
        assert mask_sensitive_sql(123) == '[MASKED_SQL]'

    def test_mask_sensitive_sql_email_pattern(self):
        """Test masking email addresses."""
        sql = "email = 'john.doe@example.com' AND name = 'John'"
        result = mask_sensitive_sql(sql)
        assert 'john.doe@example.com' not in result
        assert '[MASKED]' in result

    def test_mask_sensitive_sql_ssn_pattern(self):
        """Test masking SSN-like patterns."""
        sql = "ssn = '123-45-6789' AND id = 123456789"
        result = mask_sensitive_sql(sql)
        assert '123-45-6789' not in result
        assert '[MASKED]' in result

    def test_mask_sensitive_sql_sensitive_columns(self):
        """Test masking sensitive column values."""
        sql = "patient_id = 'P12345' AND chrom = '1'"
        result = mask_sensitive_sql(sql)
        assert 'P12345' not in result
        assert '[MASKED]' in result
        assert "chrom = '1'" in result  # Non-sensitive column should remain

    def test_mask_sensitive_sql_custom_columns(self):
        """Test masking with custom sensitive columns."""
        sql = "custom_field = 'secret' AND public_field = 'visible'"
        result = mask_sensitive_sql(sql, sensitive_columns=['custom_field'])
        assert 'secret' not in result
        assert '[MASKED]' in result
        assert 'visible' in result

    def test_mask_sensitive_sql_custom_patterns(self):
        """Test masking with custom patterns."""
        import re
        custom_patterns = [re.compile(r'SECRET_\d+')]
        sql = "field = 'SECRET_123' AND other = 'normal'"
        result = mask_sensitive_sql(sql, patterns=custom_patterns)
        assert 'SECRET_123' not in result
        assert '[MASKED]' in result

    @patch('vcf_agent.lancedb_integration.sqlparse')
    def test_mask_sensitive_sql_with_sqlparse(self, mock_sqlparse):
        """Test masking with sqlparse available."""
        # Mock sqlparse parsing
        mock_token = MagicMock()
        mock_token.ttype = 'Literal.String.Single'
        mock_token.value = "'sensitive_value'"
        
        mock_parsed = MagicMock()
        mock_parsed.flatten.return_value = [mock_token]
        
        mock_sqlparse.parse.return_value = [mock_parsed]
        mock_sqlparse.tokens.Literal.String.Single = 'Literal.String.Single'
        mock_sqlparse.tokens.Literal.Number = 'Literal.Number'
        mock_sqlparse.tokens.Name = 'Name'
        
        sql = "field = 'sensitive_value'"
        result = mask_sensitive_sql(sql)
        assert "'[MASKED]'" in result

    @patch('vcf_agent.lancedb_integration.sqlparse')
    def test_mask_sensitive_sql_sqlparse_exception(self, mock_sqlparse):
        """Test fallback when sqlparse raises exception."""
        mock_sqlparse.parse.side_effect = Exception("Parse error")
        
        sql = "email = 'test@example.com'"
        result = mask_sensitive_sql(sql)
        # Should fall back to regex masking
        assert 'test@example.com' not in result


class TestVariantModel:
    """Test cases for the Variant Pydantic model."""

    def test_variant_model_creation(self):
        """Test creating a Variant model instance."""
        variant_data = {
            "variant_id": "chr1-123-A-G",
            "embedding": [0.1] * 1024,
            "chrom": "1",
            "pos": 123,
            "ref": "A",
            "alt": "G",
            "clinical_significance": "Benign"
        }
        variant = Variant(**variant_data)
        assert variant.variant_id == "chr1-123-A-G"
        assert variant.chrom == "1"
        assert variant.pos == 123
        assert variant.ref == "A"
        assert variant.alt == "G"
        assert variant.clinical_significance == "Benign"
        assert len(variant.embedding) == 1024

    def test_variant_model_optional_fields(self):
        """Test Variant model with optional fields."""
        variant_data = {
            "variant_id": "chr1-123-A-G",
            "embedding": [0.1] * 1024,
            "chrom": "1",
            "pos": 123,
            "ref": "A",
            "alt": "G"
        }
        variant = Variant(**variant_data)
        assert variant.clinical_significance is None


class TestGetDB:
    """Test cases for the get_db function."""

    @patch('vcf_agent.lancedb_integration.lancedb.connect')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_get_db_success(self, mock_logger, mock_connect):
        """Test successful database connection."""
        mock_db = MagicMock()
        mock_connect.return_value = mock_db
        
        result = get_db("./test_db")
        
        assert result == mock_db
        mock_connect.assert_called_once_with("./test_db")
        mock_logger.info.assert_any_call("Attempting to connect to LanceDB at ./test_db...")
        mock_logger.info.assert_any_call("Successfully connected to LanceDB at ./test_db.")

    @patch('vcf_agent.lancedb_integration.lancedb.connect')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_get_db_failure(self, mock_logger, mock_connect):
        """Test database connection failure."""
        mock_connect.side_effect = Exception("Connection failed")
        
        with pytest.raises(Exception, match="Connection failed"):
            get_db("./test_db")
        
        mock_logger.error.assert_called_once()

    @patch('vcf_agent.lancedb_integration.lancedb.connect')
    def test_get_db_default_path(self, mock_connect):
        """Test get_db with default path."""
        mock_db = MagicMock()
        mock_connect.return_value = mock_db
        
        result = get_db()
        
        mock_connect.assert_called_once_with("./lancedb")


class TestGetOrCreateTable:
    """Test cases for the get_or_create_table function."""

    @patch('vcf_agent.lancedb_integration.logger')
    def test_get_or_create_table_existing(self, mock_logger):
        """Test opening existing table."""
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_db.table_names.return_value = ["variants", "other_table"]
        mock_db.open_table.return_value = mock_table
        
        result = get_or_create_table(mock_db, "variants")
        
        assert result == mock_table
        mock_db.open_table.assert_called_once_with("variants")
        mock_logger.info.assert_any_call("Opened existing table: 'variants'.")

    @patch('vcf_agent.lancedb_integration.logger')
    def test_get_or_create_table_new(self, mock_logger):
        """Test creating new table."""
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_db.table_names.return_value = ["other_table"]
        mock_db.create_table.return_value = mock_table
        
        result = get_or_create_table(mock_db, "variants")
        
        assert result == mock_table
        mock_db.create_table.assert_called_once_with("variants", schema=Variant, mode="overwrite")
        mock_logger.info.assert_any_call("Created new table: 'variants' with schema Variant.")

    @patch('vcf_agent.lancedb_integration.logger')
    def test_get_or_create_table_failure(self, mock_logger):
        """Test table creation/opening failure."""
        mock_db = MagicMock()
        mock_db.table_names.side_effect = Exception("Database error")
        
        with pytest.raises(Exception, match="Database error"):
            get_or_create_table(mock_db, "variants")
        
        mock_logger.error.assert_called_once()

    @patch('vcf_agent.lancedb_integration.logger')
    def test_get_or_create_table_default_name(self, mock_logger):
        """Test get_or_create_table with default table name."""
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_db.table_names.return_value = []
        mock_db.create_table.return_value = mock_table
        
        result = get_or_create_table(mock_db)
        
        mock_db.create_table.assert_called_once_with("variants", schema=Variant, mode="overwrite")


class TestAddVariants:
    """Test cases for the add_variants function."""

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_add_variants_success(self, mock_logger, mock_lock):
        """Test successful variant addition."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        variants = [
            {"variant_id": "chr1-123-A-G", "chrom": "1", "pos": 123},
            {"variant_id": "chr2-456-C-T", "chrom": "2", "pos": 456}
        ]
        
        add_variants(mock_table, variants)
        
        mock_table.add.assert_called_once_with(variants)
        mock_logger.info.assert_any_call("Attempting to add 2 variants to table 'test_table'...")
        mock_logger.info.assert_any_call("Successfully added 2 variants to table 'test_table'.")

    @patch('vcf_agent.lancedb_integration.logger')
    def test_add_variants_empty_list(self, mock_logger):
        """Test adding empty variant list."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        add_variants(mock_table, [])
        
        mock_table.add.assert_not_called()
        mock_logger.warning.assert_called_once()

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_add_variants_failure(self, mock_logger, mock_lock):
        """Test variant addition failure."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        mock_table.add.side_effect = Exception("Add failed")
        variants = [{"variant_id": "chr1-123-A-G"}]
        
        with pytest.raises(Exception, match="Add failed"):
            add_variants(mock_table, variants)
        
        mock_logger.error.assert_called_once()


class TestSearchByEmbedding:
    """Test cases for the search_by_embedding function."""

    @patch('vcf_agent.lancedb_integration.graph_integration')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_search_by_embedding_success(self, mock_logger, mock_graph_integration):
        """Test successful embedding search."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        # Mock search chain
        mock_query = MagicMock()
        mock_query.limit.return_value = mock_query
        mock_query.where.return_value = mock_query
        mock_table.search.return_value = mock_query
        
        # Mock results
        mock_df = pd.DataFrame([
            {"variant_id": "chr1-123-A-G", "chrom": "1", "pos": 123},
            {"variant_id": "chr2-456-C-T", "chrom": "2", "pos": 456}
        ])
        mock_query.to_pandas.return_value = mock_df
        
        # Mock Kuzu integration
        mock_kuzu_conn = MagicMock()
        mock_graph_integration.get_managed_kuzu_connection.return_value = mock_kuzu_conn
        mock_graph_integration.get_variant_context.return_value = {
            "chr1-123-A-G": [("SAMPLE_001", "HET")],
            "chr2-456-C-T": [("SAMPLE_002", "HOM")]
        }
        
        embedding = [0.1] * 1024
        result = search_by_embedding(mock_table, embedding, limit=5, filter_sql="chrom = '1'")
        
        assert len(result) == 2
        assert "kuzu_observed_samples" in result.columns
        mock_table.search.assert_called_once_with(embedding, vector_column_name='embedding')
        mock_query.limit.assert_called_once_with(5)
        mock_query.where.assert_called_once_with("chrom = '1'")

    @patch('vcf_agent.lancedb_integration.graph_integration')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_search_by_embedding_no_filter(self, mock_logger, mock_graph_integration):
        """Test embedding search without filter."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        mock_query = MagicMock()
        mock_query.limit.return_value = mock_query
        mock_table.search.return_value = mock_query
        
        mock_df = pd.DataFrame([{"variant_id": "chr1-123-A-G"}])
        mock_query.to_pandas.return_value = mock_df
        
        mock_graph_integration.get_managed_kuzu_connection.return_value = None
        
        embedding = [0.1] * 1024
        result = search_by_embedding(mock_table, embedding, limit=10)
        
        mock_query.where.assert_not_called()
        assert "kuzu_observed_samples" in result.columns

    @patch('vcf_agent.lancedb_integration.graph_integration')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_search_by_embedding_empty_results(self, mock_logger, mock_graph_integration):
        """Test embedding search with empty results."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        mock_query = MagicMock()
        mock_query.limit.return_value = mock_query
        mock_table.search.return_value = mock_query
        
        mock_df = pd.DataFrame()  # Empty DataFrame
        mock_query.to_pandas.return_value = mock_df
        
        embedding = [0.1] * 1024
        result = search_by_embedding(mock_table, embedding)
        
        assert len(result) == 0
        assert "kuzu_observed_samples" in result.columns

    @patch('vcf_agent.lancedb_integration.graph_integration')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_search_by_embedding_kuzu_error(self, mock_logger, mock_graph_integration):
        """Test embedding search with Kuzu enrichment error."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        mock_query = MagicMock()
        mock_query.limit.return_value = mock_query
        mock_table.search.return_value = mock_query
        
        mock_df = pd.DataFrame([{"variant_id": "chr1-123-A-G"}])
        mock_query.to_pandas.return_value = mock_df
        
        mock_kuzu_conn = MagicMock()
        mock_graph_integration.get_managed_kuzu_connection.return_value = mock_kuzu_conn
        mock_graph_integration.get_variant_context.side_effect = Exception("Kuzu error")
        
        embedding = [0.1] * 1024
        result = search_by_embedding(mock_table, embedding)
        
        assert len(result) == 1
        assert "kuzu_observed_samples" in result.columns
        mock_logger.error.assert_any_call("Error during Kuzu context enrichment: Kuzu error")

    @patch('vcf_agent.lancedb_integration.logger')
    def test_search_by_embedding_failure(self, mock_logger):
        """Test embedding search failure."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        mock_table.search.side_effect = Exception("Search failed")
        
        embedding = [0.1] * 1024
        with pytest.raises(Exception, match="Search failed"):
            search_by_embedding(mock_table, embedding)
        
        mock_logger.error.assert_called_once()


class TestUpdateVariant:
    """Test cases for the update_variant function."""

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_update_variant_success(self, mock_logger, mock_lock):
        """Test successful variant update."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        updates = {"clinical_significance": "Pathogenic", "pos": 12346}
        
        update_variant(mock_table, "chr1-123-A-G", updates)
        
        mock_table.update.assert_called_once_with(
            where="variant_id = 'chr1-123-A-G'",
            values=updates
        )
        mock_logger.info.assert_any_call(
            "Attempting to update variant 'chr1-123-A-G' in table 'test_table' with updates for keys: ['clinical_significance', 'pos']."
        )

    @patch('vcf_agent.lancedb_integration.logger')
    def test_update_variant_empty_updates(self, mock_logger):
        """Test update with empty updates."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        update_variant(mock_table, "chr1-123-A-G", {})
        
        mock_table.update.assert_not_called()
        mock_logger.warning.assert_called_once()

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_update_variant_failure(self, mock_logger, mock_lock):
        """Test variant update failure."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        mock_table.update.side_effect = Exception("Update failed")
        updates = {"clinical_significance": "Benign"}
        
        with pytest.raises(Exception, match="Update failed"):
            update_variant(mock_table, "chr1-123-A-G", updates)
        
        mock_logger.error.assert_called_once()


class TestDeleteVariants:
    """Test cases for the delete_variants function."""

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_delete_variants_success(self, mock_logger, mock_lock):
        """Test successful variant deletion."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        filter_sql = "chrom = '1' AND pos < 1000"
        
        delete_variants(mock_table, filter_sql)
        
        mock_table.delete.assert_called_once_with(filter_sql)
        mock_logger.info.assert_any_call(
            "Attempting to delete variants from table 'test_table' matching filter: 'chrom = '1' AND pos < 1000'."
        )

    @patch('vcf_agent.lancedb_integration.logger')
    def test_delete_variants_empty_filter(self, mock_logger):
        """Test delete with empty filter."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        with pytest.raises(ValueError, match="filter_sql cannot be empty"):
            delete_variants(mock_table, "")
        
        mock_table.delete.assert_not_called()

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_delete_variants_failure(self, mock_logger, mock_lock):
        """Test variant deletion failure."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        mock_table.delete.side_effect = Exception("Delete failed")
        filter_sql = "variant_id = 'chr1-123-A-G'"
        
        with pytest.raises(Exception, match="Delete failed"):
            delete_variants(mock_table, filter_sql)
        
        mock_logger.error.assert_called_once()


class TestCreateScalarIndex:
    """Test cases for the create_scalar_index function."""

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_create_scalar_index_success(self, mock_logger, mock_lock):
        """Test successful scalar index creation."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        create_scalar_index(mock_table, "chrom", index_type="BTREE", replace=True)
        
        mock_table.create_scalar_index.assert_called_once_with(
            "chrom", index_type="BTREE", replace=True
        )
        mock_logger.info.assert_any_call(
            "Attempting to create scalar index on column 'chrom' in table 'test_table' with type 'BTREE', replace=True, kwargs_keys=[]."
        )

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_create_scalar_index_default_type(self, mock_logger, mock_lock):
        """Test scalar index creation with default type."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        create_scalar_index(mock_table, "pos")
        
        mock_table.create_scalar_index.assert_called_once_with(
            "pos", index_type="BTREE", replace=False
        )

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_create_scalar_index_with_kwargs(self, mock_logger, mock_lock):
        """Test scalar index creation with additional kwargs."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        create_scalar_index(mock_table, "chrom", index_type="BITMAP", num_partitions=4)
        
        mock_table.create_scalar_index.assert_called_once_with(
            "chrom", index_type="BITMAP", replace=False, num_partitions=4
        )

    @patch('vcf_agent.lancedb_integration.lancedb_write_lock')
    @patch('vcf_agent.lancedb_integration.logger')
    def test_create_scalar_index_failure(self, mock_logger, mock_lock):
        """Test scalar index creation failure."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        mock_table.create_scalar_index.side_effect = Exception("Index creation failed")
        
        with pytest.raises(Exception, match="Index creation failed"):
            create_scalar_index(mock_table, "chrom")
        
        mock_logger.error.assert_called_once()


class TestFilterVariantsByMetadata:
    """Test cases for the filter_variants_by_metadata function."""

    @patch('vcf_agent.lancedb_integration.logger')
    def test_filter_variants_by_metadata_success(self, mock_logger):
        """Test successful metadata filtering."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        # Mock search chain
        mock_query = MagicMock()
        mock_query.select.return_value = mock_query
        mock_query.where.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_table.search.return_value = mock_query
        
        mock_df = pd.DataFrame([{"variant_id": "chr1-123-A-G", "chrom": "1"}])
        mock_query.to_pandas.return_value = mock_df
        
        result = filter_variants_by_metadata(
            mock_table, 
            "chrom = '1'", 
            select_columns=["variant_id", "chrom"], 
            limit=100
        )
        
        assert len(result) == 1
        mock_query.select.assert_called_once_with(["variant_id", "chrom"])
        mock_query.where.assert_called_once_with("chrom = '1'")
        mock_query.limit.assert_called_once_with(100)

    @patch('vcf_agent.lancedb_integration.logger')
    def test_filter_variants_by_metadata_no_select_columns(self, mock_logger):
        """Test metadata filtering without select columns."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        mock_query = MagicMock()
        mock_query.where.return_value = mock_query
        mock_table.search.return_value = mock_query
        
        mock_df = pd.DataFrame([{"variant_id": "chr1-123-A-G"}])
        mock_query.to_pandas.return_value = mock_df
        
        result = filter_variants_by_metadata(mock_table, "chrom = '1'")
        
        mock_query.select.assert_not_called()
        mock_query.limit.assert_not_called()

    @patch('vcf_agent.lancedb_integration.logger')
    def test_filter_variants_by_metadata_empty_filter(self, mock_logger):
        """Test metadata filtering with empty filter."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        with pytest.raises(ValueError, match="filter_sql cannot be empty"):
            filter_variants_by_metadata(mock_table, "")

    @patch('vcf_agent.lancedb_integration.logger')
    def test_filter_variants_by_metadata_unsafe_sql(self, mock_logger):
        """Test metadata filtering with unsafe SQL."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        
        unsafe_filters = [
            "chrom = '1' OR 1=1",
            "chrom = '1'; DROP TABLE variants",
            "chrom = '1' -- comment",
            "chrom = '1' AND DELETE FROM variants"
        ]
        
        for unsafe_filter in unsafe_filters:
            with pytest.raises(ValueError, match="Unsafe SQL filter detected"):
                filter_variants_by_metadata(mock_table, unsafe_filter)

    @patch('vcf_agent.lancedb_integration.logger')
    def test_filter_variants_by_metadata_failure(self, mock_logger):
        """Test metadata filtering failure."""
        mock_table = MagicMock()
        mock_table.name = "test_table"
        mock_table.search.side_effect = Exception("Filter failed")
        
        with pytest.raises(Exception, match="Filter failed"):
            filter_variants_by_metadata(mock_table, "chrom = '1'")
        
        mock_logger.error.assert_called_once()


class TestThreadSafety:
    """Test cases for thread safety features."""

    def test_lancedb_write_lock_exists(self):
        """Test that the write lock exists and is a threading.Lock."""
        import threading
        assert isinstance(lancedb_write_lock, threading.Lock)

    def test_sensitive_columns_constant(self):
        """Test that SENSITIVE_COLUMNS constant is properly defined."""
        assert isinstance(SENSITIVE_COLUMNS, list)
        assert 'patient_id' in SENSITIVE_COLUMNS
        assert 'email' in SENSITIVE_COLUMNS

    def test_sensitive_patterns_constant(self):
        """Test that SENSITIVE_PATTERNS constant is properly defined."""
        import re
        assert isinstance(SENSITIVE_PATTERNS, list)
        assert all(isinstance(pattern, re.Pattern) for pattern in SENSITIVE_PATTERNS) 