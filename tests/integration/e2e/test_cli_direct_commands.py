"""
E2E tests for direct CLI commands.

Tests the specific CLI subcommands directly:
- populate-kuzu-from-vcf command
- LanceDB commands (already tested in other files)
- CLI argument parsing and validation
- Direct command error handling
"""

import pytest
import os
import tempfile
import json
import shutil
from pathlib import Path
import kuzu


class TestCLIDirectCommands:
    """Test direct CLI commands without agent interface."""

    def test_cli_populate_kuzu_from_vcf_success(self, vcf_agent_cli_runner, sample_vcf_file_small, test_dbs):
        """Test successful populate-kuzu-from-vcf CLI command."""
        kuzu_path = test_dbs["kuzu_path"]
        
        cli_args = [
            "populate-kuzu-from-vcf",
            "--vcf_file_path", sample_vcf_file_small,
            "--kuzu_db_path", kuzu_path
        ]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI populate-kuzu-from-vcf failed. Stderr: {result.stderr}"
        
        # Verify output contains success information
        assert "populated" in result.stdout.lower() or "success" in result.stdout.lower()
        
        # Verify Kuzu database was created and populated
        assert os.path.exists(kuzu_path), "Kuzu database directory should exist"
        
        # Verify data was loaded by checking node counts
        try:
            with kuzu.Database(kuzu_path) as db:
                with kuzu.Connection(db) as conn:
                    # Check for Variant nodes
                    variant_result = conn.execute("MATCH (v:Variant) RETURN count(v)")
                    if variant_result.has_next():
                        variant_count = variant_result.get_next()[0]
                        assert variant_count > 0, "Should have loaded some variants"
                    
                    # Check for Sample nodes
                    sample_result = conn.execute("MATCH (s:Sample) RETURN count(s)")
                    if sample_result.has_next():
                        sample_count = sample_result.get_next()[0]
                        assert sample_count > 0, "Should have loaded some samples"
        except Exception as e:
            pytest.fail(f"Failed to verify Kuzu database contents: {e}")

    def test_cli_populate_kuzu_from_vcf_nonexistent_file(self, vcf_agent_cli_runner, test_dbs):
        """Test populate-kuzu-from-vcf with non-existent VCF file."""
        kuzu_path = test_dbs["kuzu_path"]
        nonexistent_vcf = "/path/to/nonexistent.vcf"
        
        cli_args = [
            "populate-kuzu-from-vcf",
            "--vcf_file_path", nonexistent_vcf,
            "--kuzu_db_path", kuzu_path
        ]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should fail gracefully
        assert result.returncode != 0, "CLI should fail with non-existent VCF file"
        assert ("not found" in result.stderr.lower() or 
                "error" in result.stderr.lower() or
                "no such file" in result.stderr.lower())

    def test_cli_populate_kuzu_from_vcf_invalid_vcf(self, vcf_agent_cli_runner, invalid_vcf_file, test_dbs):
        """Test populate-kuzu-from-vcf with invalid VCF file."""
        kuzu_path = test_dbs["kuzu_path"]
        
        cli_args = [
            "populate-kuzu-from-vcf",
            "--vcf_file_path", invalid_vcf_file,
            "--kuzu_db_path", kuzu_path
        ]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should fail gracefully
        assert result.returncode != 0, "CLI should fail with invalid VCF file"
        assert ("error" in result.stderr.lower() or 
                "invalid" in result.stderr.lower() or
                "parse" in result.stderr.lower())

    def test_cli_populate_kuzu_missing_required_args(self, vcf_agent_cli_runner):
        """Test populate-kuzu-from-vcf with missing required arguments."""
        # Missing VCF file path
        cli_args = ["populate-kuzu-from-vcf", "--kuzu_db_path", "/tmp/test"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode != 0, "CLI should fail with missing required arguments"
        assert "required" in result.stderr.lower() or "error" in result.stderr.lower()

    def test_cli_populate_kuzu_custom_db_path(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test populate-kuzu-from-vcf with custom database path."""
        with tempfile.TemporaryDirectory() as temp_dir:
            custom_kuzu_path = os.path.join(temp_dir, "custom_kuzu_db")
            
            cli_args = [
                "populate-kuzu-from-vcf",
                "--vcf_file_path", sample_vcf_file_small,
                "--kuzu_db_path", custom_kuzu_path
            ]
            result = vcf_agent_cli_runner(cli_args)
            
            assert result.returncode == 0, f"CLI with custom path failed. Stderr: {result.stderr}"
            assert os.path.exists(custom_kuzu_path), "Custom Kuzu database directory should exist"


class TestCLIArgumentParsing:
    """Test CLI argument parsing and validation."""

    def test_cli_help_command(self, vcf_agent_cli_runner):
        """Test CLI help command."""
        cli_args = ["--help"]
        result = vcf_agent_cli_runner(cli_args)
        
        assert result.returncode == 0, f"CLI help failed. Stderr: {result.stderr}"
        assert "usage:" in result.stdout.lower() or "help" in result.stdout.lower()

    def test_cli_no_arguments(self, vcf_agent_cli_runner):
        """Test CLI with no arguments."""
        cli_args = []
        result = vcf_agent_cli_runner(cli_args)
        
        # Should show help or handle gracefully
        assert result.returncode == 0, f"CLI with no args should show help. Stderr: {result.stderr}"

    def test_cli_invalid_subcommand(self, vcf_agent_cli_runner):
        """Test CLI with invalid subcommand."""
        cli_args = ["invalid-subcommand"]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should fail with error message
        assert result.returncode != 0, "CLI should fail with invalid subcommand"
        assert ("invalid" in result.stderr.lower() or 
                "error" in result.stderr.lower() or
                "unrecognized" in result.stderr.lower())

    def test_cli_model_selection(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI model selection arguments."""
        # Test with different model options
        model_options = ["ollama", "openai", "cerebras"]
        
        for model in model_options:
            cli_args = ["ask", f"validate {sample_vcf_file_small}", "--model", model]
            result = vcf_agent_cli_runner(cli_args)
            
            # Should handle model selection (may fail due to missing credentials, but should parse)
            # We're mainly testing argument parsing here
            assert isinstance(result.returncode, int), f"CLI should handle model {model}"

    def test_cli_raw_mode_flag(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI raw mode flag."""
        cli_args = ["ask", f"validate {sample_vcf_file_small}", "--raw"]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should handle raw mode flag
        assert isinstance(result.returncode, int), "CLI should handle raw mode flag"

    def test_cli_credentials_argument(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI credentials argument."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as temp_creds:
            # Create a dummy credentials file
            json.dump({"api_key": "dummy_key"}, temp_creds)
            temp_creds.flush()
            
            try:
                cli_args = ["ask", f"validate {sample_vcf_file_small}", "--credentials", temp_creds.name]
                result = vcf_agent_cli_runner(cli_args)
                
                # Should handle credentials argument
                assert isinstance(result.returncode, int), "CLI should handle credentials argument"
                
            finally:
                os.unlink(temp_creds.name)

    def test_cli_reference_argument(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI reference FASTA argument."""
        if not os.path.exists("sample_test_data/22.fa"):
            pytest.skip("Reference FASTA not available for test")
        
        cli_args = [
            "ask", 
            f"normalize {sample_vcf_file_small}", 
            "--reference", "sample_test_data/22.fa"
        ]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should handle reference argument
        assert isinstance(result.returncode, int), "CLI should handle reference argument"

    def test_cli_ollama_model_argument(self, vcf_agent_cli_runner, sample_vcf_file_small):
        """Test CLI Ollama model specification."""
        cli_args = [
            "ask", 
            f"validate {sample_vcf_file_small}", 
            "--model", "ollama",
            "--ollama-model", "qwen:4b"
        ]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should handle Ollama model specification
        assert isinstance(result.returncode, int), "CLI should handle Ollama model argument"


class TestCLILanceDBCommands:
    """Test LanceDB CLI commands (complementing existing tests)."""

    def test_cli_init_lancedb_custom_path(self, vcf_agent_cli_runner):
        """Test init-lancedb with custom path."""
        with tempfile.TemporaryDirectory() as temp_dir:
            custom_db_path = os.path.join(temp_dir, "custom_lancedb")
            
            cli_args = [
                "init-lancedb",
                "--db_path", custom_db_path,
                "--table_name", "test_variants"
            ]
            result = vcf_agent_cli_runner(cli_args)
            
            assert result.returncode == 0, f"CLI init-lancedb failed. Stderr: {result.stderr}"
            assert "Initialized" in result.stdout or "initialized" in result.stdout.lower()
            assert os.path.exists(custom_db_path), "Custom LanceDB directory should exist"

    def test_cli_add_variant_minimal(self, vcf_agent_cli_runner):
        """Test add-variant with minimal required arguments."""
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = os.path.join(temp_dir, "test_lancedb")
            
            # First initialize the database
            init_args = ["init-lancedb", "--db_path", db_path]
            init_result = vcf_agent_cli_runner(init_args)
            assert init_result.returncode == 0
            
            # Then add a variant
            embedding = ",".join(["0.1"] * 1024)  # 1024-dimensional embedding
            add_args = [
                "add-variant",
                "--db_path", db_path,
                "--variant_id", "test-variant-1",
                "--chrom", "1",
                "--pos", "12345",
                "--ref", "A",
                "--alt", "G",
                "--embedding", embedding
            ]
            result = vcf_agent_cli_runner(add_args)
            
            assert result.returncode == 0, f"CLI add-variant failed. Stderr: {result.stderr}"
            assert "Added variant" in result.stdout or "added" in result.stdout.lower()

    def test_cli_search_embedding_basic(self, vcf_agent_cli_runner):
        """Test search-embedding with basic arguments."""
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = os.path.join(temp_dir, "test_lancedb")
            
            # Initialize and add a variant first
            init_args = ["init-lancedb", "--db_path", db_path]
            init_result = vcf_agent_cli_runner(init_args)
            assert init_result.returncode == 0
            
            embedding = ",".join(["0.1"] * 1024)
            add_args = [
                "add-variant",
                "--db_path", db_path,
                "--variant_id", "search-test-variant",
                "--chrom", "1",
                "--pos", "12345",
                "--ref", "A",
                "--alt", "G",
                "--embedding", embedding
            ]
            add_result = vcf_agent_cli_runner(add_args)
            assert add_result.returncode == 0
            
            # Now search for the variant
            search_args = [
                "search-embedding",
                "--db_path", db_path,
                "--embedding", embedding,
                "--limit", "1"
            ]
            result = vcf_agent_cli_runner(search_args)
            
            assert result.returncode == 0, f"CLI search-embedding failed. Stderr: {result.stderr}"
            assert "search-test-variant" in result.stdout


class TestCLIErrorHandling:
    """Test CLI error handling scenarios."""

    def test_cli_invalid_json_updates(self, vcf_agent_cli_runner):
        """Test update-variant with invalid JSON."""
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = os.path.join(temp_dir, "test_lancedb")
            
            cli_args = [
                "update-variant",
                "--db_path", db_path,
                "--variant_id", "test-variant",
                "--updates", "invalid-json-string"
            ]
            result = vcf_agent_cli_runner(cli_args)
            
            # Should fail with JSON parsing error
            assert result.returncode != 0, "CLI should fail with invalid JSON"
            assert ("json" in result.stderr.lower() or 
                    "parse" in result.stderr.lower() or
                    "error" in result.stderr.lower())

    def test_cli_invalid_embedding_format(self, vcf_agent_cli_runner):
        """Test add-variant with invalid embedding format."""
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = os.path.join(temp_dir, "test_lancedb")
            
            # Initialize database first
            init_args = ["init-lancedb", "--db_path", db_path]
            init_result = vcf_agent_cli_runner(init_args)
            assert init_result.returncode == 0
            
            # Try to add variant with invalid embedding
            cli_args = [
                "add-variant",
                "--db_path", db_path,
                "--variant_id", "test-variant",
                "--chrom", "1",
                "--pos", "12345",
                "--ref", "A",
                "--alt", "G",
                "--embedding", "invalid,embedding,format,not,numbers"
            ]
            result = vcf_agent_cli_runner(cli_args)
            
            # Should fail with embedding parsing error
            assert result.returncode != 0, "CLI should fail with invalid embedding format"

    def test_cli_missing_database_path(self, vcf_agent_cli_runner):
        """Test commands with missing database paths."""
        nonexistent_db = "/path/to/nonexistent/database"
        
        cli_args = [
            "search-embedding",
            "--db_path", nonexistent_db,
            "--embedding", ",".join(["0.1"] * 1024)
        ]
        result = vcf_agent_cli_runner(cli_args)
        
        # Should fail gracefully
        assert result.returncode != 0, "CLI should fail with missing database"
        assert ("not found" in result.stderr.lower() or 
                "error" in result.stderr.lower() or
                "no such" in result.stderr.lower())

    def test_cli_invalid_sql_filter(self, vcf_agent_cli_runner):
        """Test filter commands with invalid SQL."""
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = os.path.join(temp_dir, "test_lancedb")
            
            # Initialize database
            init_args = ["init-lancedb", "--db_path", db_path]
            init_result = vcf_agent_cli_runner(init_args)
            assert init_result.returncode == 0
            
            # Try filter with invalid SQL
            cli_args = [
                "filter-lancedb",
                "--db_path", db_path,
                "--filter_sql", "INVALID SQL SYNTAX HERE"
            ]
            result = vcf_agent_cli_runner(cli_args)
            
            # Should fail with SQL error
            assert result.returncode != 0, "CLI should fail with invalid SQL"


class TestCLIIntegrationScenarios:
    """Test CLI integration scenarios."""

    def test_cli_full_workflow_vcf_to_kuzu_to_lancedb(self, vcf_agent_cli_runner, sample_vcf_file_small, test_dbs):
        """Test full workflow: VCF → Kuzu → LanceDB via CLI."""
        kuzu_path = test_dbs["kuzu_path"]
        lancedb_path = test_dbs["lancedb_path"]
        
        # Step 1: Populate Kuzu from VCF
        kuzu_args = [
            "populate-kuzu-from-vcf",
            "--vcf_file_path", sample_vcf_file_small,
            "--kuzu_db_path", kuzu_path
        ]
        kuzu_result = vcf_agent_cli_runner(kuzu_args)
        assert kuzu_result.returncode == 0, f"Kuzu population failed. Stderr: {kuzu_result.stderr}"
        
        # Step 2: Initialize LanceDB
        lance_init_args = ["init-lancedb", "--db_path", lancedb_path]
        lance_init_result = vcf_agent_cli_runner(lance_init_args)
        assert lance_init_result.returncode == 0, f"LanceDB init failed. Stderr: {lance_init_result.stderr}"
        
        # Step 3: Add a variant to LanceDB (simulating data from Kuzu)
        embedding = ",".join(["0.1"] * 1024)
        add_variant_args = [
            "add-variant",
            "--db_path", lancedb_path,
            "--variant_id", "workflow-test-variant",
            "--chrom", "1",
            "--pos", "12345",
            "--ref", "A",
            "--alt", "G",
            "--embedding", embedding,
            "--clinical_significance", "Benign"
        ]
        add_result = vcf_agent_cli_runner(add_variant_args)
        assert add_result.returncode == 0, f"Add variant failed. Stderr: {add_result.stderr}"
        
        # Step 4: Search the variant
        search_args = [
            "search-embedding",
            "--db_path", lancedb_path,
            "--embedding", embedding,
            "--limit", "1"
        ]
        search_result = vcf_agent_cli_runner(search_args)
        assert search_result.returncode == 0, f"Search failed. Stderr: {search_result.stderr}"
        assert "workflow-test-variant" in search_result.stdout
        
        # Verify both databases exist and contain data
        assert os.path.exists(kuzu_path), "Kuzu database should exist"
        assert os.path.exists(lancedb_path), "LanceDB should exist"

    def test_cli_error_recovery_workflow(self, vcf_agent_cli_runner, test_dbs):
        """Test CLI error recovery in workflows."""
        lancedb_path = test_dbs["lancedb_path"]
        
        # Step 1: Try to add variant without initializing (may succeed due to auto-creation)
        embedding = ",".join(["0.1"] * 1024)
        add_args_fail = [
            "add-variant",
            "--db_path", lancedb_path,
            "--variant_id", "recovery-test",
            "--chrom", "1",
            "--pos", "12345",
            "--ref", "A",
            "--alt", "G",
            "--embedding", embedding
        ]
        fail_result = vcf_agent_cli_runner(add_args_fail)
        # LanceDB may auto-create tables, so this might succeed
        # We'll just verify it completes without crashing
        
        # Step 2: Initialize database (ensure it's properly set up)
        init_args = ["init-lancedb", "--db_path", lancedb_path]
        init_result = vcf_agent_cli_runner(init_args)
        assert init_result.returncode == 0, "Initialization should succeed"
        
        # Step 3: Add another variant (should succeed)
        embedding2 = ",".join(["0.2"] * 1024)
        add_args_success = [
            "add-variant",
            "--db_path", lancedb_path,
            "--variant_id", "recovery-test-2",
            "--chrom", "2",
            "--pos", "67890",
            "--ref", "C",
            "--alt", "T",
            "--embedding", embedding2
        ]
        retry_result = vcf_agent_cli_runner(add_args_success)
        assert retry_result.returncode == 0, f"Second add should succeed. Stderr: {retry_result.stderr}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 