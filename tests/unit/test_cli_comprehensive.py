"""
Comprehensive unit tests for CLI module.

Tests the command-line interface including argument parsing, subcommands,
LanceDB operations, Kuzu operations, and agent interactions.
"""

import pytest
import sys
import os
import json
from unittest.mock import patch, MagicMock, Mock, call
from io import StringIO
from vcf_agent.cli import main


class TestCLIMockResponse:
    """Test cases for CLI mock response functionality."""

    @patch.dict(os.environ, {'VCF_AGENT_CLI_MOCK_RESPONSE': 'Mock response for testing'})
    @patch('builtins.print')
    def test_cli_mock_response(self, mock_print):
        """Test CLI with mock response environment variable."""
        # Execute
        main()

        # Verify
        mock_print.assert_called_once_with('Mock response for testing')

    @patch.dict(os.environ, {}, clear=True)
    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.cli.argparse.ArgumentParser.print_help')
    def test_cli_no_mock_response_no_command(self, mock_print_help, mock_parse_args):
        """Test CLI without mock response and no command."""
        # Setup mock to return args with no command
        mock_args = MagicMock()
        mock_args.command = None
        mock_parse_args.return_value = mock_args

        # Execute and verify SystemExit
        with pytest.raises(SystemExit) as exc_info:
            main()
        
        assert exc_info.value.code == 0
        mock_print_help.assert_called_once()


class TestCLILanceDBCommands:
    """Test cases for LanceDB CLI commands."""

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('builtins.print')
    def test_init_lancedb_command(self, mock_print, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test init-lancedb command."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "init-lancedb"
        mock_args.db_path = "./test_lancedb"
        mock_args.table_name = "test_variants"
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table

        # Execute
        main()

        # Verify
        mock_get_db.assert_called_once_with("./test_lancedb")
        mock_get_or_create_table.assert_called_once_with(mock_db, "test_variants")
        mock_print.assert_called_with("Initialized LanceDB table 'test_variants' at './test_lancedb'")

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.add_variants')
    @patch('builtins.print')
    def test_add_variant_command(self, mock_print, mock_add_variants, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test add-variant command."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "add-variant"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.variant_id = "rs123"
        mock_args.chrom = "1"
        mock_args.pos = 12345
        mock_args.ref = "A"
        mock_args.alt = "G"
        mock_args.embedding = "0.1,0.2,0.3"
        mock_args.clinical_significance = "Benign"
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table

        # Execute
        main()

        # Verify
        mock_get_db.assert_called_once_with("./lancedb")
        mock_get_or_create_table.assert_called_once_with(mock_db, "variants")
        expected_variant = {
            "variant_id": "rs123",
            "chrom": "1",
            "pos": 12345,
            "ref": "A",
            "alt": "G",
            "embedding": [0.1, 0.2, 0.3],
            "clinical_significance": "Benign",
        }
        mock_add_variants.assert_called_once_with(mock_table, [expected_variant])
        mock_print.assert_called_with("Added variant rs123 to LanceDB.")

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.search_by_embedding')
    @patch('builtins.print')
    def test_search_embedding_command(self, mock_print, mock_search_by_embedding, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test search-embedding command."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "search-embedding"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.embedding = "0.1,0.2,0.3"
        mock_args.limit = 5
        mock_args.filter_sql = "chrom = '1'"
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table
        mock_search_results = [{"variant_id": "rs123", "distance": 0.1}]
        mock_search_by_embedding.return_value = mock_search_results

        # Execute
        main()

        # Verify
        mock_get_db.assert_called_once_with("./lancedb")
        mock_get_or_create_table.assert_called_once_with(mock_db, "variants")
        mock_search_by_embedding.assert_called_once_with(mock_table, [0.1, 0.2, 0.3], limit=5, filter_sql="chrom = '1'")
        mock_print.assert_called_with(mock_search_results)

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.update_variant')
    @patch('builtins.print')
    def test_update_variant_command(self, mock_print, mock_update_variant, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test update-variant command."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "update-variant"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.variant_id = "rs123"
        mock_args.updates = '{"clinical_significance": "Pathogenic", "pos": 12346}'
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table

        # Execute
        main()

        # Verify
        mock_get_db.assert_called_once_with("./lancedb")
        mock_get_or_create_table.assert_called_once_with(mock_db, "variants")
        expected_updates = {"clinical_significance": "Pathogenic", "pos": 12346}
        mock_update_variant.assert_called_once_with(mock_table, "rs123", expected_updates)
        mock_print.assert_called_with("Update command processed for variant rs123.")

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('sys.stderr', new_callable=StringIO)
    def test_update_variant_invalid_json(self, mock_stderr, mock_parse_args):
        """Test update-variant command with invalid JSON."""
        # Setup mock args with invalid JSON
        mock_args = MagicMock()
        mock_args.command = "update-variant"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.variant_id = "rs123"
        mock_args.updates = '{"invalid": json}'
        mock_parse_args.return_value = mock_args

        # Execute and verify SystemExit
        with pytest.raises(SystemExit) as exc_info:
            main()
        
        assert exc_info.value.code == 1
        assert "Error: Invalid JSON string for updates:" in mock_stderr.getvalue()

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.delete_variants')
    @patch('builtins.print')
    def test_delete_variants_command(self, mock_print, mock_delete_variants, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test delete-variants command."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "delete-variants"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.filter_sql = "chrom = '1' AND pos < 1000"
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table

        # Execute
        main()

        # Verify
        mock_get_db.assert_called_once_with("./lancedb")
        mock_get_or_create_table.assert_called_once_with(mock_db, "variants")
        mock_delete_variants.assert_called_once_with(mock_table, "chrom = '1' AND pos < 1000")
        mock_print.assert_called_with("Delete command processed with filter: chrom = '1' AND pos < 1000")

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.create_scalar_index')
    @patch('builtins.print')
    def test_create_lancedb_index_command(self, mock_print, mock_create_scalar_index, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test create-lancedb-index command."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "create-lancedb-index"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.column = "chrom"
        mock_args.index_type = "BTREE"
        mock_args.replace = True
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table

        # Execute
        main()

        # Verify
        mock_get_db.assert_called_once_with("./lancedb")
        mock_get_or_create_table.assert_called_once_with(mock_db, "variants")
        mock_create_scalar_index.assert_called_once_with(mock_table, "chrom", index_type="BTREE", replace=True)
        mock_print.assert_called_with("Index creation command processed for column 'chrom' on table 'variants'.")

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.filter_variants_by_metadata')
    @patch('builtins.print')
    def test_filter_lancedb_command(self, mock_print, mock_filter_variants, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test filter-lancedb command."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "filter-lancedb"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.filter_sql = "chrom = '1'"
        mock_args.select_columns = "variant_id,chrom,pos"
        mock_args.limit = 100
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table
        
        # Mock DataFrame result
        import pandas as pd
        mock_df = pd.DataFrame([{"variant_id": "rs123", "chrom": "1", "pos": 12345}])
        mock_filter_variants.return_value = mock_df

        # Execute
        main()

        # Verify
        mock_get_db.assert_called_once_with("./lancedb")
        mock_get_or_create_table.assert_called_once_with(mock_db, "variants")
        mock_filter_variants.assert_called_once_with(
            mock_table, 
            "chrom = '1'", 
            select_columns=["variant_id", "chrom", "pos"], 
            limit=100
        )
        # Verify print calls
        assert mock_print.call_count == 2
        mock_print.assert_any_call("Found 1 variants matching filter:")

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.filter_variants_by_metadata')
    @patch('sys.stderr', new_callable=StringIO)
    def test_filter_lancedb_unsafe_sql_error(self, mock_stderr, mock_filter_variants, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test filter-lancedb command with unsafe SQL error."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "filter-lancedb"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.filter_sql = "DROP TABLE variants"
        mock_args.select_columns = None
        mock_args.limit = None
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table
        mock_filter_variants.side_effect = ValueError("Unsafe SQL filter detected: DROP")

        # Execute and verify SystemExit
        with pytest.raises(SystemExit) as exc_info:
            main()
        
        assert exc_info.value.code == 1
        assert "Error: Unsafe SQL filter detected: DROP" in mock_stderr.getvalue()


class TestCLIKuzuCommands:
    """Test cases for Kuzu CLI commands."""

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.graph_integration.get_kuzu_db_connection')
    @patch('vcf_agent.graph_integration.create_schema')
    @patch('vcf_agent.vcf_utils.populate_kuzu_from_vcf')
    @patch('vcf_agent.cli.metrics')
    @patch('vcf_agent.cli.cli_tracer')
    @patch('vcf_agent.cli.time')
    def test_populate_kuzu_from_vcf_success(self, mock_time, mock_tracer, mock_metrics, mock_populate, mock_create_schema, mock_get_kuzu_connection, mock_parse_args):
        """Test populate-kuzu-from-vcf command success."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "populate-kuzu-from-vcf"
        mock_args.vcf_file_path = "/path/to/test.vcf"
        mock_args.kuzu_db_path = "./kuzu_db"
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_connection = MagicMock()
        mock_get_kuzu_connection.return_value = mock_connection
        mock_populate.return_value = {"variants": 100, "samples": 5, "links": 150}
        
        # Mock time
        mock_time.time.side_effect = [1000.0, 1010.0]  # 10 second duration
        
        # Mock tracer
        mock_span = MagicMock()
        mock_tracer.start_as_current_span.return_value.__enter__.return_value = mock_span
        
        # Mock metrics
        mock_duration_metric = MagicMock()
        mock_requests_metric = MagicMock()
        mock_metrics.CLI_COMMAND_DURATION_SECONDS = mock_duration_metric
        mock_metrics.CLI_COMMAND_REQUESTS_TOTAL = mock_requests_metric

        # Execute
        main()

        # Verify
        mock_get_kuzu_connection.assert_called_once_with(db_path="./kuzu_db")
        mock_create_schema.assert_called_once_with(mock_connection)
        mock_populate.assert_called_once_with(kuzu_conn=mock_connection, vcf_path="/path/to/test.vcf")
        
        # Verify span attributes
        mock_span.set_attribute.assert_any_call("vcf.file_path", "/path/to/test.vcf")
        mock_span.set_attribute.assert_any_call("kuzu.db_path", "./kuzu_db")
        mock_span.set_attribute.assert_any_call("kuzu.variants_added", 100)
        mock_span.set_attribute.assert_any_call("kuzu.samples_added_found", 5)
        mock_span.set_attribute.assert_any_call("kuzu.links_created", 150)

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.graph_integration.get_kuzu_db_connection')
    @patch('vcf_agent.graph_integration.create_schema')
    @patch('vcf_agent.vcf_utils.populate_kuzu_from_vcf')
    @patch('vcf_agent.cli.metrics')
    @patch('vcf_agent.cli.cli_tracer')
    @patch('vcf_agent.cli.time')
    @patch('sys.stderr', new_callable=StringIO)
    def test_populate_kuzu_from_vcf_file_not_found(self, mock_stderr, mock_time, mock_tracer, mock_metrics, mock_populate, mock_create_schema, mock_get_kuzu_connection, mock_parse_args):
        """Test populate-kuzu-from-vcf command with file not found error."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "populate-kuzu-from-vcf"
        mock_args.vcf_file_path = "/path/to/nonexistent.vcf"
        mock_args.kuzu_db_path = "./kuzu_db"
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_connection = MagicMock()
        mock_get_kuzu_connection.return_value = mock_connection
        mock_populate.side_effect = FileNotFoundError("VCF file not found")
        
        # Mock time
        mock_time.time.side_effect = [1000.0, 1010.0]
        
        # Mock tracer
        mock_span = MagicMock()
        mock_tracer.start_as_current_span.return_value.__enter__.return_value = mock_span
        
        # Mock metrics
        mock_duration_metric = MagicMock()
        mock_requests_metric = MagicMock()
        mock_errors_metric = MagicMock()
        mock_metrics.CLI_COMMAND_DURATION_SECONDS = mock_duration_metric
        mock_metrics.CLI_COMMAND_REQUESTS_TOTAL = mock_requests_metric
        mock_metrics.CLI_COMMAND_ERRORS_TOTAL = mock_errors_metric

        # Execute and verify SystemExit
        with pytest.raises(SystemExit) as exc_info:
            main()
        
        assert exc_info.value.code == 1
        
        # Verify error handling
        mock_span.record_exception.assert_called()
        mock_span.set_status.assert_called()


class TestCLIAgentCommands:
    """Test cases for agent CLI commands."""

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.agent.get_agent_with_session')
    @patch('vcf_agent.config.SessionConfig')
    @patch('vcf_agent.cli.cli_tracer')
    @patch('vcf_agent.cli.metrics')
    @patch('vcf_agent.cli.time')
    @patch('builtins.print')
    def test_ask_command_success(self, mock_print, mock_time, mock_metrics, mock_tracer, mock_session_config, mock_get_agent, mock_parse_args):
        """Test ask command success."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "ask"
        mock_args.prompt_text = "Validate this VCF file"
        mock_args.raw = False
        mock_args.model = "ollama"
        mock_args.credentials = None
        mock_args.reference = None
        mock_args.ollama_model = None
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_agent_instance = MagicMock()
        mock_agent_instance.return_value = "VCF file is valid"
        mock_get_agent.return_value = mock_agent_instance
        
        # Mock time
        mock_time.time.side_effect = [1000.0, 1005.0, 1065.0]  # 5 second AI response, 60 second sleep
        
        # Mock tracer
        mock_span = MagicMock()
        mock_tracer.start_as_current_span.return_value.__enter__.return_value = mock_span
        
        # Mock metrics
        mock_ai_response_metric = MagicMock()
        mock_ai_requests_metric = MagicMock()
        mock_metrics.VCF_AGENT_AI_RESPONSE_SECONDS = mock_ai_response_metric
        mock_metrics.VCF_AGENT_AI_REQUESTS_TOTAL = mock_ai_requests_metric

        # Execute
        main()

        # Verify
        mock_get_agent.assert_called_once()
        mock_agent_instance.assert_called_once_with("Validate this VCF file")
        mock_print.assert_any_call("VCF file is valid")
        
        # Verify span attributes
        mock_span.set_attribute.assert_any_call("agent.prompt", "Validate this VCF file")
        mock_span.set_attribute.assert_any_call("agent.model_provider", "ollama")
        mock_span.set_attribute.assert_any_call("agent.raw_mode", False)

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.agent.get_agent_with_session')
    @patch('vcf_agent.config.SessionConfig')
    @patch('vcf_agent.cli.cli_tracer')
    @patch('vcf_agent.cli.metrics')
    @patch('vcf_agent.cli.time')
    @patch.dict(os.environ, {'VCF_AGENT_RAW_MODE': '1'})
    def test_ask_command_raw_mode(self, mock_time, mock_metrics, mock_tracer, mock_session_config, mock_get_agent, mock_parse_args):
        """Test ask command with raw mode."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "ask"
        mock_args.prompt_text = "Quick validation"
        mock_args.raw = True
        mock_args.model = "openai"
        mock_args.credentials = "/path/to/creds.json"
        mock_args.reference = "/path/to/ref.fasta"
        mock_args.ollama_model = "custom-model"
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_agent_instance = MagicMock()
        mock_response_obj = MagicMock()
        mock_response_obj.response = "Quick validation result"
        mock_agent_instance.return_value = mock_response_obj
        mock_get_agent.return_value = mock_agent_instance
        
        # Mock time
        mock_time.time.side_effect = [1000.0, 1003.0, 1063.0]
        
        # Mock tracer
        mock_span = MagicMock()
        mock_tracer.start_as_current_span.return_value.__enter__.return_value = mock_span

        # Execute
        main()

        # Verify environment variable was set
        assert os.environ.get("VCF_AGENT_RAW_MODE") == "1"
        
        # Verify SessionConfig was called with correct parameters
        mock_session_config.assert_called_once_with(
            raw_mode=True,
            model_provider="openai",
            credentials_file="/path/to/creds.json",
            reference_fasta="/path/to/ref.fasta",
            ollama_model_name="custom-model"
        )

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.agent.get_agent_with_session')
    @patch('vcf_agent.config.SessionConfig')
    @patch('vcf_agent.cli.cli_tracer')
    @patch('vcf_agent.cli.metrics')
    @patch('vcf_agent.cli.time')
    @patch('sys.stderr', new_callable=StringIO)
    def test_ask_command_agent_error(self, mock_stderr, mock_time, mock_metrics, mock_tracer, mock_session_config, mock_get_agent, mock_parse_args):
        """Test ask command with agent error."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "ask"
        mock_args.prompt_text = "Invalid request"
        mock_args.raw = False
        mock_args.model = "ollama"
        mock_args.credentials = None
        mock_args.reference = None
        mock_args.ollama_model = None
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_agent_instance = MagicMock()
        mock_agent_instance.side_effect = Exception("Agent processing error")
        mock_get_agent.return_value = mock_agent_instance
        
        # Mock time
        mock_time.time.side_effect = [1000.0, 1005.0]
        
        # Mock tracer
        mock_span = MagicMock()
        mock_tracer.start_as_current_span.return_value.__enter__.return_value = mock_span
        
        # Mock metrics
        mock_ai_response_metric = MagicMock()
        mock_ai_requests_metric = MagicMock()
        mock_ai_errors_metric = MagicMock()
        mock_metrics.VCF_AGENT_AI_RESPONSE_SECONDS = mock_ai_response_metric
        mock_metrics.VCF_AGENT_AI_REQUESTS_TOTAL = mock_ai_requests_metric
        mock_metrics.VCF_AGENT_AI_ERRORS_TOTAL = mock_ai_errors_metric

        # Execute and verify exception is raised
        with pytest.raises(Exception, match="Agent processing error"):
            main()
        
        # Verify error handling
        mock_span.record_exception.assert_called()
        mock_span.set_status.assert_called()


class TestCLITracingAndMetrics:
    """Test cases for CLI tracing and metrics functionality."""

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('opentelemetry.trace.get_tracer_provider')
    @patch('opentelemetry.sdk.trace.TracerProvider')
    @patch('vcf_agent.metrics.start_metrics_http_server')
    @patch('vcf_agent.cli.time.sleep')
    @patch('builtins.print')
    def test_tracing_flush_success(self, mock_print, mock_sleep, mock_start_metrics, mock_sdk_tracer, mock_get_tracer_provider, mock_parse_args):
        """Test successful tracing flush."""
        # Setup mock args for ask command
        mock_args = MagicMock()
        mock_args.command = "ask"
        mock_args.prompt_text = "test"
        mock_args.raw = False
        mock_args.model = "ollama"
        mock_args.credentials = None
        mock_args.reference = None
        mock_args.ollama_model = None
        mock_parse_args.return_value = mock_args

        # Setup tracer provider mock
        mock_provider = MagicMock()
        mock_get_tracer_provider.return_value = mock_provider
        
        # Make isinstance check return True
        with patch('vcf_agent.cli.isinstance', return_value=True):
            # Mock other dependencies
            with patch('vcf_agent.agent.get_agent_with_session'), \
                 patch('vcf_agent.config.SessionConfig'), \
                 patch('vcf_agent.cli.cli_tracer'), \
                 patch('vcf_agent.cli.metrics'), \
                 patch('vcf_agent.cli.time'):
                
                # Execute
                main()

        # Verify flush was called
        mock_provider.force_flush.assert_called_once_with(timeout_millis=5000)
        mock_print.assert_any_call("\n[CLI] Forcing flush of OTel spans...")
        mock_print.assert_any_call("[CLI] OTel spans flushed.")

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('opentelemetry.trace.get_tracer_provider')
    @patch('builtins.print')
    def test_tracing_flush_no_sdk_provider(self, mock_print, mock_get_tracer_provider, mock_parse_args):
        """Test tracing flush with no SDK provider."""
        # Setup mock args for ask command
        mock_args = MagicMock()
        mock_args.command = "ask"
        mock_args.prompt_text = "test"
        mock_args.raw = False
        mock_args.model = "ollama"
        mock_args.credentials = None
        mock_args.reference = None
        mock_args.ollama_model = None
        mock_parse_args.return_value = mock_args

        # Setup tracer provider mock (not SDK provider)
        mock_provider = MagicMock()
        mock_get_tracer_provider.return_value = mock_provider
        
        # Make isinstance check return False
        with patch('vcf_agent.cli.isinstance', return_value=False):
            # Mock other dependencies
            with patch('vcf_agent.agent.get_agent_with_session'), \
                 patch('vcf_agent.config.SessionConfig'), \
                 patch('vcf_agent.cli.cli_tracer'), \
                 patch('vcf_agent.cli.metrics'), \
                 patch('vcf_agent.cli.time'):
                
                # Execute
                main()

        # Verify no flush was attempted
        mock_provider.force_flush.assert_not_called()
        mock_print.assert_any_call("\n[CLI] No SDK TracerProvider found, skipping flush.")


class TestCLIEdgeCases:
    """Test edge cases and error conditions for CLI."""

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.update_variant')
    @patch('sys.stderr', new_callable=StringIO)
    def test_update_variant_non_dict_updates(self, mock_stderr, mock_update_variant, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test update-variant command with non-dictionary updates."""
        # Setup mock args with array instead of object
        mock_args = MagicMock()
        mock_args.command = "update-variant"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.variant_id = "rs123"
        mock_args.updates = '["not", "a", "dict"]'
        mock_parse_args.return_value = mock_args

        # Execute and verify SystemExit
        with pytest.raises(SystemExit) as exc_info:
            main()
        
        assert exc_info.value.code == 1
        assert "Error: Updates must be a valid JSON object (dictionary)." in mock_stderr.getvalue()

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.update_variant')
    @patch('sys.stderr', new_callable=StringIO)
    def test_update_variant_update_error(self, mock_stderr, mock_update_variant, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test update-variant command with update operation error."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "update-variant"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.variant_id = "rs123"
        mock_args.updates = '{"pos": 12346}'
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table
        mock_update_variant.side_effect = Exception("Update operation failed")

        # Execute and verify SystemExit
        with pytest.raises(SystemExit) as exc_info:
            main()
        
        assert exc_info.value.code == 1
        assert "Error during variant update: Update operation failed" in mock_stderr.getvalue()

    @patch('vcf_agent.cli.argparse.ArgumentParser.parse_args')
    @patch('vcf_agent.lancedb_integration.get_db')
    @patch('vcf_agent.lancedb_integration.get_or_create_table')
    @patch('vcf_agent.lancedb_integration.filter_variants_by_metadata')
    @patch('sys.stderr', new_callable=StringIO)
    def test_filter_lancedb_unexpected_error(self, mock_stderr, mock_filter_variants, mock_get_or_create_table, mock_get_db, mock_parse_args):
        """Test filter-lancedb command with unexpected error."""
        # Setup mock args
        mock_args = MagicMock()
        mock_args.command = "filter-lancedb"
        mock_args.db_path = "./lancedb"
        mock_args.table_name = "variants"
        mock_args.filter_sql = "chrom = '1'"
        mock_args.select_columns = None
        mock_args.limit = None
        mock_parse_args.return_value = mock_args

        # Setup mocks
        mock_db = MagicMock()
        mock_table = MagicMock()
        mock_get_db.return_value = mock_db
        mock_get_or_create_table.return_value = mock_table
        mock_filter_variants.side_effect = Exception("Unexpected database error")

        # Execute and verify SystemExit
        with pytest.raises(SystemExit) as exc_info:
            main()
        
        assert exc_info.value.code == 1
        assert "An unexpected error occurred: Unexpected database error" in mock_stderr.getvalue() 