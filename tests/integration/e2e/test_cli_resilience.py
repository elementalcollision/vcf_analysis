import pytest
from pathlib import Path
import os
import shutil # For reliably removing directories
from typing import Iterator, Union, Type # Added for type hint
import numpy as np # For mock embeddings
from cyvcf2 import VCF # For VCF parsing
import kuzu # For Kuzu specific exceptions if needed
from _pytest.python_api import RaisesContext # For explicit typing of excinfo

# Assuming vcf_agent.lancedb_integration and vcf_agent.graph_integration are importable
from vcf_agent.lancedb_integration import get_db as get_lancedb_connection, get_or_create_table, add_variants, Variant
from vcf_agent import graph_integration # For monkeypatching
from vcf_agent import vcf_utils # For populate_kuzu_from_vcf
from vcf_agent.graph_integration import get_kuzu_db_connection, create_schema as create_kuzu_schema

# Helper function to run CLI commands (can be adapted from other test files or made more generic)
# from vcf_agent.cli import main as cli_main # Or use subprocess if preferred for full isolation

# Fixture for a temporary, populated LanceDB for testing search resilience
@pytest.fixture
def temp_populated_lancedb(tmp_path: Path) -> Iterator[Path]:
    db_path = tmp_path / "temp_lancedb_for_search"
    db_path.mkdir(parents=True, exist_ok=True)
    lancedb_conn = None
    try:
        lancedb_conn = get_lancedb_connection(str(db_path))
        table = get_or_create_table(lancedb_conn, "variants")
        variants_data = [
            Variant(variant_id="rs123", chrom="1", pos=1000, ref="A", alt="T", embedding=[0.1]*1024).model_dump(),
            Variant(variant_id="rs456", chrom="2", pos=2000, ref="C", alt="G", embedding=[0.2]*1024).model_dump(),
        ]
        add_variants(table, variants_data)
        yield db_path
    finally:
        # LanceDB connection doesn't have a close method in the typical sense for file-based DBs.
        # Cleanup is handled by removing the directory.
        if db_path.exists():
            shutil.rmtree(db_path)

# --- Error Handling Tests (LanceDB Down - Refactored to use read_only_dir) ---

def test_cli_init_lancedb_down_read_only(vcf_agent_cli_runner, read_only_dir: Path):
    db_path_str = str(read_only_dir / "init_lancedb_ro_test")
    result = vcf_agent_cli_runner(["init-lancedb", "--db_path", db_path_str])
    assert result.returncode != 0
    assert any(msg.lower() in result.stderr.lower() for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError"])

def test_cli_search_embedding_lancedb_down_read_only(vcf_agent_cli_runner, read_only_dir: Path):
    db_path_str = str(read_only_dir / "search_lancedb_ro_test")
    # No need to create the dir, as it's read-only and connection should fail
    dummy_embedding = ",".join(["0.1"] * 1024) 
    result = vcf_agent_cli_runner([
        "search-embedding", 
        "--db_path", db_path_str, 
        "--embedding", dummy_embedding
    ])
    assert result.returncode != 0
    assert any(msg.lower() in result.stderr.lower() for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError", "No such file or directory"])

def test_cli_filter_lancedb_down_read_only(vcf_agent_cli_runner, read_only_dir: Path):
    db_path_str = str(read_only_dir / "filter_lancedb_ro_test")
    result = vcf_agent_cli_runner([
        "filter-lancedb", 
        "--db_path", db_path_str,
        "--filter_sql", "chrom = '1'"
    ])
    assert result.returncode != 0
    assert any(msg.lower() in result.stderr.lower() for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError", "No such file or directory"])

def test_cli_add_variant_lancedb_down_read_only(vcf_agent_cli_runner, read_only_dir: Path):
    db_path_str = str(read_only_dir / "add_variant_lancedb_ro_test")
    dummy_embedding = ",".join(["0.1"] * 1024)
    result = vcf_agent_cli_runner([
        "add-variant", 
        "--db_path", db_path_str, 
        "--variant_id", "rs123", 
        "--chrom", "1", 
        "--pos", "1000", 
        "--ref", "A", 
        "--alt", "T", 
        "--embedding", dummy_embedding
    ])
    assert result.returncode != 0
    assert any(msg.lower() in result.stderr.lower() for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError"])

def test_cli_update_variant_lancedb_down_read_only(vcf_agent_cli_runner, read_only_dir: Path):
    db_path_str = str(read_only_dir / "update_variant_lancedb_ro_test")
    result = vcf_agent_cli_runner([
        "update-variant", 
        "--db_path", db_path_str, 
        "--variant_id", "rs123", 
        "--updates", '{"pos": 1001}'
    ])
    assert result.returncode != 0
    assert any(msg.lower() in result.stderr.lower() for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError"])

def test_cli_delete_variants_lancedb_down_read_any(vcf_agent_cli_runner, read_only_dir: Path):
    db_path_str = str(read_only_dir / "delete_variants_lancedb_ro_test")
    result = vcf_agent_cli_runner([
        "delete-variants", 
        "--db_path", db_path_str, 
        "--filter_sql", "variant_id = 'rs123'"
    ])
    assert result.returncode != 0
    assert any(msg.lower() in result.stderr.lower() for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError"])

def test_cli_create_lancedb_index_down_read_only(vcf_agent_cli_runner, read_only_dir: Path):
    db_path_str = str(read_only_dir / "create_index_lancedb_ro_test")
    result = vcf_agent_cli_runner([
        "create-lancedb-index", 
        "--db_path", db_path_str, 
        "--column", "variant_id"
    ])
    assert result.returncode != 0
    assert any(msg.lower() in result.stderr.lower() for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError"])

# Keep the original test_cli_init_lancedb_down that uses unavailable_lancedb_dir
# as it tests a different scenario (recovery by creation).

def test_cli_init_lancedb_down(vcf_agent_cli_runner, unavailable_lancedb_dir: Path):
    """
    Test Scenario 1.1 (Error Handling - Recovery):
    CLI command 'init-lancedb' when the target LanceDB directory is made unavailable by renaming.
    `init-lancedb` should succeed by recreating the directory.
    """
    db_path_str = str(unavailable_lancedb_dir) 
    assert not unavailable_lancedb_dir.exists(), "conftest fixture should have made this path unavailable."
    result = vcf_agent_cli_runner(["init-lancedb", "--db_path", db_path_str])
    assert result.returncode == 0, (
        f"CLI init-lancedb should succeed by creating the directory. "
        f"Stderr: {result.stderr}, Stdout: {result.stdout}"
    )
    assert unavailable_lancedb_dir.exists(), "init-lancedb should have recreated the directory."

# --- Error Handling Tests (Kuzu Down) ---

def test_cli_search_embedding_kuzu_down(vcf_agent_cli_runner, temp_populated_lancedb: Path, unavailable_kuzu_dir: Path, monkeypatch):
    """
    Test Scenario (Error Handling K.1):
    CLI command \'search-embedding\' when Kuzu is unavailable for enrichment.
    LanceDB search should succeed, but Kuzu enrichment should fail gracefully with a warning.
    `unavailable_kuzu_dir` from conftest creates a Kuzu DB then renames its directory.
    """
    lancedb_path_str = str(temp_populated_lancedb)
    kuzu_db_path_str = str(unavailable_kuzu_dir) 
    assert not unavailable_kuzu_dir.is_dir(), "conftest fixture should have made Kuzu path not a directory."
    
    monkeypatch.setattr(graph_integration, "DEFAULT_KUZU_DB_PATH", kuzu_db_path_str)
    
    query_embedding = ",".join(["0.1"] * 1024)
    result = vcf_agent_cli_runner([
        "search-embedding", 
        "--db_path", lancedb_path_str, 
        "--embedding", query_embedding,
        "--table_name", "variants"
    ])

    assert result.returncode == 0, f"CLI should exit with 0 on LanceDB success even if Kuzu fails. Got: {result.returncode}\\nStderr: {result.stderr}\\nStdout: {result.stdout}"
    assert "rs123" in result.stdout, f"LanceDB result 'rs123' not found in stdout. Got: {result.stdout}"
    assert "Kuzu connection not available, skipping enrichment" in result.stderr or            "Error during Kuzu context enrichment" in result.stderr or            "Kuzu setup failed" in result.stderr or            "No such file or directory" in result.stderr or            "runtime error" in result.stderr.lower(),         f"Stderr should indicate a Kuzu enrichment failure. Got: {result.stderr}"

def test_api_populate_kuzu_from_vcf_kuzu_down(unavailable_kuzu_dir: Path, sample_vcf_file_small: Path, capsys, monkeypatch):
    """
    Test Scenario (Error Handling K.2):
    Python API vcf_utils.populate_kuzu_from_vcf when Kuzu DB dir is unavailable during writes.
    `unavailable_kuzu_dir` from conftest creates a Kuzu DB then renames its directory.
    Expects failure at kuzu.Database() or kuzu.Connection() time.
    """
    kuzu_db_path_str = str(unavailable_kuzu_dir) 
    assert not unavailable_kuzu_dir.is_dir(), "conftest fixture should have made Kuzu path not a directory."
    kuzu_conn = None

    monkeypatch.setattr(graph_integration, "DEFAULT_KUZU_DB_PATH", kuzu_db_path_str)
    
    try:
        with pytest.raises((RuntimeError, OSError)) as excinfo: # type: ignore
            kuzu_db_obj_down = kuzu.Database(kuzu_db_path_str) # This should fail
            kuzu_conn = kuzu.Connection(kuzu_db_obj_down)      # Or this if Database() call was lenient
            # The following lines should not be reached
            # graph_integration.create_schema(kuzu_conn)
            # vcf_utils.populate_kuzu_from_vcf(kuzu_conn, str(sample_vcf_file_small))

        assert excinfo is not None
        # Check for messages indicating the path is bad or connection failed
        error_msg = str(excinfo.value).lower()
        assert "no such file or directory" in error_msg or \
               "failed to connect" in error_msg or \
               "database not found" in error_msg or \
               "runtime error" in error_msg or \
               "io error" in error_msg or \
               "not a directory" in error_msg or \
               "cannot open file" in error_msg, \
               f"Expected DB connection/creation error due to unavailable path, got: {error_msg}"

    finally:
        if kuzu_conn: # Should be None if the above `raises` worked as expected
            try:
                kuzu_conn.close()
            except Exception:
                pass

# --- Security Check Tests ---

@pytest.fixture
def populated_lancedb_for_sqli_tests(tmp_path: Path) -> Iterator[Path]:
    """
    Creates and populates a temporary LanceDB for SQL injection tests.
    Yields the path to the database directory.
    """
    db_path = tmp_path / "sqli_test_db"
    db_path.mkdir(parents=True, exist_ok=True)
    lancedb_conn = None
    table_name = "variants"

    try:
        lancedb_conn = get_lancedb_connection(str(db_path))
        table = get_or_create_table(lancedb_conn, table_name)
        sample_variants = [
            Variant(variant_id='record1', chrom='1', pos=100, ref='A', alt='T', embedding=[0.1]*1024).model_dump(),
            Variant(variant_id='record2', chrom='2', pos=200, ref='C', alt='G', embedding=[0.2]*1024).model_dump(),
            Variant(variant_id='sensitive_record_test', chrom='3', pos=300, ref='G', alt='A', embedding=[0.3]*1024).model_dump(),
        ]
        add_variants(table, sample_variants)
        yield db_path
    finally:
        if db_path.exists():
            shutil.rmtree(db_path)

def test_cli_filter_lancedb_sqli_always_true(vcf_agent_cli_runner, populated_lancedb_for_sqli_tests: Path):
    """
    Test Scenario (Security S.1):
    Attempts an 'OR \'1\'=\'1\'' type SQL injection via filter-lancedb CLI.
    LanceDB/DuckDB should treat this as a literal or error out, not dump the table.
    """
    db_path_str = str(populated_lancedb_for_sqli_tests)
    table_name = "variants"
    filter_string = "variant_id = 'nonexistent' OR '1'='1'"
    
    result = vcf_agent_cli_runner([
        "filter-lancedb", 
        "--db_path", db_path_str, 
        "--table_name", table_name,
        "--filter_sql", filter_string
    ])
    
    # With the new heuristic, this type of SQLi should be rejected by the application layer.
    assert result.returncode != 0, "CLI should fail due to potentially unsafe SQL filter."
    assert "unsafe sql filter detected" in result.stderr.lower() or "query rejected" in result.stderr.lower(), \
        f"Expected an error message about unsafe SQL. Got: Stderr: {result.stderr}, Stdout: {result.stdout}"

def test_cli_filter_lancedb_sqli_comment_attack(vcf_agent_cli_runner, populated_lancedb_for_sqli_tests: Path):
    """
    Test Scenario (Security S.2):
    Attempts a SQL injection using comments (e.g., to append DROP TABLE) via filter-lancedb CLI.
    LanceDB/DuckDB should reject multi-statement queries or treat comments appropriately.
    The critical check is that the table is NOT dropped.
    """
    db_path_str = str(populated_lancedb_for_sqli_tests)
    table_name = "variants"
    malicious_filter = rf"variant_id = 'record1'; -- DROP TABLE {table_name}; --"

    result_injection_attempt = vcf_agent_cli_runner([
        "filter-lancedb",
        "--db_path", db_path_str,
        "--table_name", table_name,
        "--filter_sql", malicious_filter
    ])

    if result_injection_attempt.returncode == 0:
        assert "record1" in result_injection_attempt.stdout,             "If SQLi comment attack succeeded without error, it should have found 'record1' based on the benign part."
        assert "DROP TABLE" not in result_injection_attempt.stderr.upper(),             "No DROP TABLE error should appear if the command succeeded, implying it was ignored."
    else:
        assert "Error" in result_injection_attempt.stderr or                "syntax error" in result_injection_attempt.stderr.lower() or                "Cannot execute multiple statements" in result_injection_attempt.stderr or                "Invalid Input" in result_injection_attempt.stderr,             f"Expected an error for SQLi comment attack. Got: {result_injection_attempt.stderr}"

    result_verify = vcf_agent_cli_runner([
        "filter-lancedb",
        "--db_path", db_path_str,
        "--table_name", table_name,
        "--filter_sql", "variant_id = 'record2'"
    ])

    assert result_verify.returncode == 0, "Verification query failed (exit code)."
    assert "record2" in result_verify.stdout, "Verification query failed (record2 not found)."
    
    verify_sensitive_result = vcf_agent_cli_runner([
        "filter-lancedb", 
        "--db_path", db_path_str, 
        "--table_name", table_name, 
        "--filter_sql", "variant_id = 'sensitive_record_test'"
    ])
    assert verify_sensitive_result.returncode == 0, "Second verification query failed (exit code)."
    assert "sensitive_record_test" in verify_sensitive_result.stdout, "Second verification query failed (sensitive_record_test not found)."

# --- Filesystem/Permission Error Tests ---

@pytest.fixture
def read_only_dir(tmp_path: Path) -> Iterator[Path]:
    """
    Creates a temporary directory, makes it read-only, and yields its path.
    Restores write permissions for cleanup.
    """
    ro_dir = tmp_path / "read_only_test_dir"
    ro_dir.mkdir(parents=True, exist_ok=True)
    original_mode = ro_dir.stat().st_mode
    read_only_mode = 0o555  # Read and execute for all
    
    try:
        ro_dir.chmod(read_only_mode)
        yield ro_dir
    finally:
        ro_dir.chmod(original_mode)
        if ro_dir.exists():
            shutil.rmtree(ro_dir, ignore_errors=True)

def test_cli_add_variant_read_only_path(vcf_agent_cli_runner, read_only_dir: Path):
    """
    Test Scenario (Error Handling FS.2):
    CLI command \'add-variant\' when the target LanceDB directory is read-only.
    """
    db_path_str = str(read_only_dir / "test_db_for_add") 
    dummy_embedding = ",".join(["0.1"] * 1024)
    result = vcf_agent_cli_runner([
        "add-variant", 
        "--db_path", db_path_str, 
        "--variant_id", "rs123", 
        "--chrom", "1", 
        "--pos", "1000", 
        "--ref", "A", 
        "--alt", "T", 
        "--embedding", dummy_embedding
    ])
    assert result.returncode != 0
    assert any(msg in result.stderr for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError"]), (
        f"Stderr should indicate a permission/read-only error. Got: {result.stderr}"
    )

def test_cli_update_variant_read_only_path(vcf_agent_cli_runner, read_only_dir: Path):
    """
    Test Scenario (Error Handling FS.3):
    CLI command \'update-variant\' when the target LanceDB directory is read-only.
    """
    db_path_str = str(read_only_dir / "test_db_for_update")
    result = vcf_agent_cli_runner([
        "update-variant", 
        "--db_path", db_path_str, 
        "--variant_id", "rs123", 
        "--updates", '{"pos": 1001}'
    ])
    assert result.returncode != 0
    assert any(msg in result.stderr for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError"]), (
        f"Stderr should indicate a permission/read-only error. Got: {result.stderr}"
    )

def test_cli_delete_variants_read_only_path(vcf_agent_cli_runner, read_only_dir: Path):
    """
    Test Scenario (Error Handling FS.4):
    CLI command \'delete-variants\' when the target LanceDB directory is read-only.
    """
    db_path_str = str(read_only_dir / "test_db_for_delete")
    result = vcf_agent_cli_runner([
        "delete-variants", 
        "--db_path", db_path_str, 
        "--filter_sql", "variant_id = 'rs123'"
    ])
    assert result.returncode != 0
    assert any(msg in result.stderr for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError"]), (
        f"Stderr should indicate a permission/read-only error. Got: {result.stderr}"
    )

def test_cli_create_lancedb_index_read_only_path(vcf_agent_cli_runner, read_only_dir: Path):
    """
    Test Scenario (Error Handling FS.5):
    CLI command \'create-lancedb-index\' when the target LanceDB directory is read-only.
    """
    db_path_str = str(read_only_dir / "test_db_for_index")
    result = vcf_agent_cli_runner([
        "create-lancedb-index", 
        "--db_path", db_path_str, 
        "--column", "variant_id"
    ])
    assert result.returncode != 0
    assert any(msg in result.stderr for msg in ["Error", "Failed", "Permission denied", "Read-only", "OperationalError"]), (
        f"Stderr should indicate a permission/read-only error. Got: {result.stderr}"
    )

def test_api_populate_kuzu_read_only_path(read_only_dir: Path, sample_vcf_file_small: Path, capsys):
    """
    Test Scenario (Error Handling FS.6):
    Python API vcf_utils.populate_kuzu_from_vcf when Kuzu DB dir is read-only.
    Expects graceful failure (low/zero counts) and logged errors/exceptions.
    """
    kuzu_db_path_str = str(read_only_dir / "test_kuzu_data_ro")

    kuzu_conn = None
    try:
        with pytest.raises((RuntimeError, OSError)) as excinfo: # type: ignore
            kuzu_db_obj = kuzu.Database(kuzu_db_path_str) 
            kuzu_conn = kuzu.Connection(kuzu_db_obj)
            create_kuzu_schema(kuzu_conn)
            vcf_utils.populate_kuzu_from_vcf(kuzu_conn, str(sample_vcf_file_small))
        
        assert excinfo is not None 
        assert any(keyword in str(excinfo.value).lower() for keyword in 
                   ["read-only", "permission denied", "os error", "runtime error", "io error", "storage error", "cannot open file", "failed to create directory"]), (
            f"Exception message should indicate a read-only/permission issue. Got: {str(excinfo.value)}"
        )

    except Exception as e:
        assert any(keyword in str(e).lower() for keyword in 
                   ["read-only", "permission denied", "os error", "runtime error", "io error", "storage error", "cannot open file", "failed to create directory"]), (
            f"Exception message should indicate a read-only/permission issue. Got: {str(e)}"
        )
    finally:
        if kuzu_conn:
            try:
                kuzu_conn.close()
            except Exception:
                pass 

def test_api_ingest_vcf_lancedb_read_only(read_only_dir: Path, tmp_path: Path, sample_vcf_file_small: Path, capsys):
    """
    Test Scenario (Error Handling FS.7 - Refactored to API):
    API ingestion when the LanceDB path is read-only.
    Kuzu path is writable.
    """
    lancedb_ro_path_str = str(read_only_dir / "lancedb_data_ro_fs7")
    kuzu_writable_path_str = str(tmp_path / "kuzu_data_writable_fs7")
    os.makedirs(kuzu_writable_path_str, exist_ok=True)

    lancedb_conn = None
    kuzu_conn_obj = None
    try:
        kuzu_db = kuzu.Database(kuzu_writable_path_str)
        kuzu_conn_obj = kuzu.Connection(kuzu_db)
        create_kuzu_schema(kuzu_conn_obj)

        with pytest.raises((OSError, RuntimeError)) as excinfo_lance: # type: ignore
            lancedb_conn = get_lancedb_connection(lancedb_ro_path_str)
            table_obj = get_or_create_table(lancedb_conn, "variants") 
            dummy_variant_model = Variant(variant_id="rsTest", chrom="1", pos=100, ref="A", alt="T", embedding=np.array([0.1]*1024))
            add_variants(table_obj, [dummy_variant_model.model_dump()])

        assert excinfo_lance is not None
        assert any(keyword in str(excinfo_lance.value).lower() for keyword in 
                   ["read-only", "permission denied", "os error", "runtime error", "io error"]), (
            f"LanceDB error should indicate read-only/permission issue. Got: {str(excinfo_lance.value)}"
        )
        
    finally:
        if lancedb_conn:
            pass 
        if kuzu_conn_obj:
            try:
                kuzu_conn_obj.close()
            except Exception: 
                pass

def test_api_ingest_vcf_kuzu_read_only(tmp_path: Path, read_only_dir: Path, sample_vcf_file_small: Path, capsys):
    """
    Test Scenario (Error Handling FS.8 - Refactored to API):
    API ingestion when the Kuzu path is read-only.
    LanceDB path is writable.
    """
    lancedb_writable_path_str = str(tmp_path / "lancedb_data_writable_fs8")
    kuzu_ro_path_str = str(read_only_dir / "kuzu_data_ro_fs8")

    lancedb_conn = None
    kuzu_conn_obj = None
    try:
        lancedb_conn = get_lancedb_connection(lancedb_writable_path_str)
        table = get_or_create_table(lancedb_conn, "variants")

        with pytest.raises((RuntimeError, OSError)) as excinfo_kuzu: # type: ignore
            kuzu_db = kuzu.Database(kuzu_ro_path_str)
            kuzu_conn_obj = kuzu.Connection(kuzu_db)
            create_kuzu_schema(kuzu_conn_obj)
            vcf_utils.populate_kuzu_from_vcf(kuzu_conn_obj, str(sample_vcf_file_small))

        assert excinfo_kuzu is not None
        assert any(keyword in str(excinfo_kuzu.value).lower() for keyword in 
                   ["read-only", "permission denied", "os error", "runtime error", "io error", "storage error", "cannot open file", "failed to create directory"]), (
            f"Kuzu error should indicate read-only/permission issue. Got: {str(excinfo_kuzu.value)}"
        )
        
        dummy_variant_model = Variant(variant_id="rsTestFS8", chrom="1", pos=100, ref="A", alt="T", embedding=np.array([0.1]*1024))
        add_variants(table, [dummy_variant_model.model_dump()])
        assert table.search().limit(1).to_pandas()["variant_id"][0] == "rsTestFS8"

    finally:
        if lancedb_conn:
            pass 
        if kuzu_conn_obj:
            try:
                kuzu_conn_obj.close()
            except Exception: 
                pass

def test_api_ingest_vcf_both_read_only(read_only_dir: Path, sample_vcf_file_small: Path, capsys):
    """
    Test Scenario (Error Handling FS.9 - Refactored to API):
    API ingestion when both LanceDB and Kuzu paths are read-only.
    """
    lancedb_ro_path_str = str(read_only_dir / "lancedb_data_ro_fs9")
    kuzu_ro_path_str = str(read_only_dir / "kuzu_data_ro_fs9") 

    lancedb_conn = None
    kuzu_conn_obj = None
    expected_exceptions_both_ro = (OSError, RuntimeError) 

    with pytest.raises(expected_exceptions_both_ro) as excinfo: # type: ignore
        try:
            # Attempt Kuzu first as it might fail earlier with stricter path checks for DB creation
            try:
                kuzu_db = kuzu.Database(kuzu_ro_path_str)
                kuzu_conn_obj = kuzu.Connection(kuzu_db)
                create_kuzu_schema(kuzu_conn_obj)
                # If Kuzu somehow succeeds, try populating
                if kuzu_conn_obj:
                    vcf_utils.populate_kuzu_from_vcf(kuzu_conn_obj, str(sample_vcf_file_small))
            except (OSError, RuntimeError) as kuzu_e:
                # If Kuzu fails, this is an expected outcome for a read-only path.
                # We re-raise to be caught by the outer pytest.raises
                raise kuzu_e

            # If Kuzu part didn't raise, proceed to LanceDB (which should also fail)
            lancedb_conn = get_lancedb_connection(lancedb_ro_path_str)
            table_obj_lance = get_or_create_table(lancedb_conn, "variants") # This should fail
            
            if table_obj_lance: # Should not be reached if get_or_create_table fails
                dummy_variant_model = Variant(variant_id="rsTestFS9", chrom="1", pos=100, ref="A", alt="T", embedding=np.array([0.1]*1024))
                add_variants(table_obj_lance, [dummy_variant_model.model_dump()]) 
            
            # If neither Kuzu nor LanceDB operations raised an error (highly unlikely)
            if not excinfo: 
                 pytest.fail("Expected an exception due to read-only paths, but none was raised.")

        finally: # This finally is for the inner try related to DB operations
            if lancedb_conn:
                pass
            if kuzu_conn_obj:
                try:
                    kuzu_conn_obj.close()
                except Exception: 
                    pass
    
    assert excinfo is not None
    assert any(keyword in str(excinfo.value).lower() for keyword in 
               ["read-only", "permission denied", "os error", "runtime error", "io error", "storage error", "cannot open file", "failed to create directory"]), (
        f"Error should indicate read-only/permission issue for one or both DBs. Got: {str(excinfo.value)}"
    )

def test_api_ingest_vcf_unreadable_vcf(tmp_path: Path, unreadable_file: Path, capsys):
    """
    Test Scenario (Error Handling FS.10 - Refactored to API):
    API ingestion when the VCF file itself is unreadable.
    """
    lancedb_path_str = str(tmp_path / "lancedb_data_fs10")
    kuzu_path_str = str(tmp_path / "kuzu_data_fs10")
    os.makedirs(kuzu_path_str, exist_ok=True)
    os.makedirs(lancedb_path_str, exist_ok=True)

    lancedb_conn = None
    kuzu_conn_obj = None
    
    expected_exceptions_unreadable_vcf = (OSError, ValueError, RuntimeError)

    with pytest.raises(expected_exceptions_unreadable_vcf) as excinfo: # type: ignore
        try:
            kuzu_db = kuzu.Database(kuzu_path_str)
            kuzu_conn_obj = kuzu.Connection(kuzu_db)
            create_kuzu_schema(kuzu_conn_obj)

            lancedb_conn = get_lancedb_connection(lancedb_path_str)
            _ = get_or_create_table(lancedb_conn, "variants")

            vcf_utils.populate_kuzu_from_vcf(kuzu_conn_obj, str(unreadable_file))
            
        finally:
            if lancedb_conn:
                pass
            if kuzu_conn_obj:
                try:
                    kuzu_conn_obj.close()
                except Exception: 
                    pass

    assert excinfo is not None
    assert any(msg in str(excinfo.value).lower() for msg in ["error", "failed", "permission denied", "cannot open", "oserror", "ioerror"]),(
        f"Stderr should indicate a file read/permission error for the VCF. Got: {str(excinfo.value)}"
    )

# --- CLI Argument/Input Data Error Tests ---

def test_cli_add_variant_invalid_pos_type(vcf_agent_cli_runner, tmp_path: Path):
    """
    Test Scenario (Error Handling AR.1):
    CLI command \'add-variant\' with an invalid data type for \'--pos\'
    Argparse should catch this and return a user-friendly error.
    """
    db_path_str = str(tmp_path / "test_db_for_add_invalid_pos")
    dummy_embedding = ",".join(["0.1"] * 1024)
    result = vcf_agent_cli_runner([
        "add-variant", 
        "--db_path", db_path_str, 
        "--variant_id", "rs_invalid_pos", 
        "--chrom", "1", 
        "--pos", "not_an_integer", 
        "--ref", "A", 
        "--alt", "T", 
        "--embedding", dummy_embedding
    ])
    assert result.returncode != 0, "CLI should fail due to invalid argument type for --pos"
    assert "invalid" in result.stderr.lower() and ("int" in result.stderr.lower() or "argument type" in result.stderr.lower() or "--pos" in result.stderr.lower()), (
        f"Stderr should indicate an invalid argument type error for --pos. Got: {result.stderr}"
    )

def test_cli_update_variant_malformed_json_updates(vcf_agent_cli_runner, tmp_path: Path):
    """
    Test Scenario (Error Handling AR.2):
    CLI command \'update-variant\' with a malformed JSON string for \'--updates\'
    The CLI handler for this command should catch the JSONDecodeError.
    """
    db_path_str = str(tmp_path / "test_db_for_update_malformed_json")
    os.makedirs(db_path_str, exist_ok=True)

    malformed_json = "{'pos': 1001, 'broken_json" 
    result = vcf_agent_cli_runner([
        "update-variant", 
        "--db_path", db_path_str, 
        "--variant_id", "rs_malformed_json", 
        "--updates", malformed_json
    ])
    assert result.returncode != 0, "CLI should fail due to malformed JSON for --updates"
    assert "invalid json" in result.stderr.lower() or "jsondecodeerror" in result.stderr.lower() or "expecting property name" in result.stderr.lower(), (
        f"Stderr should indicate a JSON decoding error. Got: {result.stderr}" # stdout might also contain it
    )