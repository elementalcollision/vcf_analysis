import os
import shutil
import subprocess
import tempfile
import numpy as np
import pytest
import logging
import json # For CLI update test

from vcf_agent.lancedb_integration import (
    get_db, 
    get_or_create_table, 
    add_variants, 
    search_by_embedding, 
    update_variant, 
    delete_variants,
    Variant,
    create_scalar_index,
    filter_variants_by_metadata
)

@pytest.fixture(scope="module")
def lancedb_tmpdir():
    tmpdir = tempfile.mkdtemp(prefix="lancedb_test_")
    yield tmpdir
    shutil.rmtree(tmpdir)

@pytest.fixture
def test_variant_data():
    return {
        "variant_id": "rsTest1",
        "chrom": "1",
        "pos": 12345,
        "ref": "A",
        "alt": "G",
        "embedding": np.random.rand(1024).astype(np.float32).tolist(),
        "clinical_significance": "Pathogenic",
    }

@pytest.fixture
def caplog_fixture(caplog):
    caplog.set_level(logging.INFO)
    return caplog

def test_lancedb_python_api_happy_path(lancedb_tmpdir, test_variant_data, caplog_fixture):
    db = get_db(lancedb_tmpdir)
    assert f"Successfully connected to LanceDB at {lancedb_tmpdir}" in caplog_fixture.text
    table = get_or_create_table(db)
    assert "Created new table: 'variants' with schema Variant." in caplog_fixture.text
    
    add_variants(table, [test_variant_data])
    assert f"Successfully added 1 variants to table '{table.name}'." in caplog_fixture.text

    results = search_by_embedding(table, test_variant_data["embedding"], limit=1)
    assert f"Search in table '{table.name}' completed. Found 1 results for limit: 1, filter_sql: '[MASKED_SQL]'." in caplog_fixture.text
    assert not results.empty
    assert results.iloc[0]["variant_id"] == test_variant_data["variant_id"]

      # db.close()

def test_lancedb_python_api_update_variant(lancedb_tmpdir, test_variant_data, caplog_fixture):
    db = get_db(lancedb_tmpdir)
    table_name = "update_test_table"
    if table_name in db.table_names():
        db.drop_table(table_name)
    table = get_or_create_table(db, table_name=table_name)
    
    add_variants(table, [test_variant_data])
    
    variant_id_to_update = test_variant_data["variant_id"]
    updates = {"clinical_significance": "Benign", "pos": 54321}
    update_variant(table, variant_id_to_update, updates)
    assert f"Successfully updated variant '{variant_id_to_update}' in table '{table_name}' with keys: ['clinical_significance', 'pos']" in caplog_fixture.text
    
    # Verify update by searching for the record and checking updated fields
    results = table.search().where(f"variant_id = '{variant_id_to_update}'").to_pandas()
    assert not results.empty
    assert results.iloc[0]["clinical_significance"] == "Benign"
    assert results.iloc[0]["pos"] == 54321

def test_lancedb_python_api_delete_variant(lancedb_tmpdir, test_variant_data, caplog_fixture):
    db = get_db(lancedb_tmpdir)
    table_name = "delete_test_table"
    if table_name in db.table_names():
        db.drop_table(table_name)
    table = get_or_create_table(db, table_name=table_name)
    
    # Add multiple variants for deletion test
    variant1 = test_variant_data.copy()
    variant2 = test_variant_data.copy()
    variant2["variant_id"] = "rsTest2"
    variant2["clinical_significance"] = "Uncertain"
    add_variants(table, [variant1, variant2])
    assert table.count_rows() == 2
    
    # Delete one variant by ID
    delete_variants(table, f"variant_id = '{variant1['variant_id']}'")
    assert f"Successfully deleted variants from table '{table.name}' matching filter: 'variant_id = \'{variant1['variant_id']}\''." in caplog_fixture.text
    assert table.count_rows() == 1
    
    # Verify only the correct variant was deleted
    remaining_results = table.search().to_pandas()
    assert not remaining_results.empty
    assert remaining_results.iloc[0]["variant_id"] == variant2["variant_id"]

    # Test deleting with a non-existent ID (should not fail, just log)
    delete_variants(table, "variant_id = 'non_existent_id'")
    assert table.count_rows() == 1 # Count should remain the same

def test_add_variants_empty_list(lancedb_tmpdir, caplog_fixture):
    db = get_db(lancedb_tmpdir)
    # Use a unique table name for this test to ensure it starts empty
    empty_test_table_name = "empty_test_table"
    # Ensure the table is created fresh for this test
    if empty_test_table_name in db.table_names():
        db.drop_table(empty_test_table_name)
    table = get_or_create_table(db, table_name=empty_test_table_name)
    
    initial_count = table.count_rows()
    add_variants(table, [])
    assert f"add_variants called for table '{empty_test_table_name}' with an empty list. No data will be added." in caplog_fixture.text
    assert table.count_rows() == initial_count  # No rows should be added

def test_add_variants_schema_mismatch(lancedb_tmpdir, caplog_fixture):
    db = get_db(lancedb_tmpdir)
    table = get_or_create_table(db)
    bad_variant_data = {"wrong_field": "value"}
    with pytest.raises(Exception): # LanceDB might raise a specific error, adjust if known
        add_variants(table, [bad_variant_data])
    assert f"Failed to add batch of 1 variants to table '{table.name}': Field 'wrong_field' not found in target schema" in caplog_fixture.text

# Simplified CLI tests, focusing on core functionality with less stdout parsing
def test_lancedb_cli_init(tmp_path):
    db_path = tmp_path / "lancedb_cli_init"
    result = subprocess.run([
        "python", "-m", "vcf_agent.cli", "init-lancedb",
        "--db_path", str(db_path)
    ], capture_output=True, text=True, check=True)
    # Make assertion more flexible, just check for key phrases
    assert "Initialized LanceDB table" in result.stdout
    assert "variants" in result.stdout # Check that the default table name is mentioned

def test_lancedb_cli_add_and_search(tmp_path, test_variant_data):
    db_path = tmp_path / "lancedb_cli_as"
    subprocess.run([
        "python", "-m", "vcf_agent.cli", "init-lancedb",
        "--db_path", str(db_path)
    ], check=True)

    embedding_str = ",".join(map(str, test_variant_data["embedding"]))
    add_command = [
        "python", "-m", "vcf_agent.cli", "add-variant",
        "--db_path", str(db_path),
        "--variant_id", test_variant_data["variant_id"],
        "--chrom", test_variant_data["chrom"],
        "--pos", str(test_variant_data["pos"]),
        "--ref", test_variant_data["ref"],
        "--alt", test_variant_data["alt"],
        "--embedding", embedding_str,
        "--clinical_significance", test_variant_data["clinical_significance"]
    ]
    result = subprocess.run(add_command, capture_output=True, text=True, check=True)
    assert f"Added variant {test_variant_data['variant_id']}" in result.stdout

    search_command = [
        "python", "-m", "vcf_agent.cli", "search-embedding",
        "--db_path", str(db_path),
        "--embedding", embedding_str,
        "--limit", "1"
    ]
    result = subprocess.run(search_command, capture_output=True, text=True, check=True)
    assert test_variant_data["variant_id"] in result.stdout 

def test_lancedb_cli_update_and_delete(tmp_path, test_variant_data):
    db_path = tmp_path / "lancedb_cli_upd_del"
    subprocess.run([
        "python", "-m", "vcf_agent.cli", "init-lancedb",
        "--db_path", str(db_path)
    ], check=True)

    # Add initial variant via CLI
    embedding_str = ",".join(map(str, test_variant_data["embedding"]))
    add_command = [
        "python", "-m", "vcf_agent.cli", "add-variant",
        "--db_path", str(db_path),
        "--variant_id", test_variant_data["variant_id"],
        "--chrom", test_variant_data["chrom"],
        "--pos", str(test_variant_data["pos"]),
        "--ref", test_variant_data["ref"],
        "--alt", test_variant_data["alt"],
        "--embedding", embedding_str,
        "--clinical_significance", test_variant_data["clinical_significance"]
    ]
    subprocess.run(add_command, check=True)

    # Update variant via CLI
    variant_id_to_update = test_variant_data["variant_id"]
    updates_dict = {"clinical_significance": "Likely Benign", "pos": 11111}
    updates_json_str = json.dumps(updates_dict)
    update_command = [
        "python", "-m", "vcf_agent.cli", "update-variant",
        "--db_path", str(db_path),
        "--variant_id", variant_id_to_update,
        "--updates", updates_json_str
    ]
    result = subprocess.run(update_command, capture_output=True, text=True, check=True)
    assert f"Update command processed for variant {variant_id_to_update}" in result.stdout

    # Verify update by searching (Python API for simplicity here, or complex CLI search)
    db = get_db(db_path)
    table = get_or_create_table(db)
    updated_results = table.search().where(f"variant_id = '{variant_id_to_update}'").to_pandas()
    assert not updated_results.empty
    assert updated_results.iloc[0]["clinical_significance"] == "Likely Benign"
    assert updated_results.iloc[0]["pos"] == 11111

    # Delete variant via CLI
    delete_filter = f"variant_id = '{variant_id_to_update}'"
    delete_command = [
        "python", "-m", "vcf_agent.cli", "delete-variants",
        "--db_path", str(db_path),
        "--filter_sql", delete_filter
    ]
    result = subprocess.run(delete_command, capture_output=True, text=True, check=True)
    assert f"Delete command processed with filter: {delete_filter}" in result.stdout

    # Verify deletion
    # Re-open the table to get the latest state after CLI modification
    db_after_delete = get_db(db_path)
    table_after_delete = get_or_create_table(db_after_delete) # Assuming default table name 'variants'
    assert table_after_delete.count_rows() == 0

def test_lancedb_python_api_create_scalar_index(lancedb_tmpdir, test_variant_data, caplog_fixture):
    db = get_db(lancedb_tmpdir)
    table_name = "index_test_table_py"
    if table_name in db.table_names():
        db.drop_table(table_name)
    table = get_or_create_table(db, table_name=table_name)
    add_variants(table, [test_variant_data])

    column_to_index = "clinical_significance"
    create_scalar_index(table, column_to_index)
    assert f"Successfully created/updated scalar index on column '{column_to_index}'" in caplog_fixture.text
    
    # Verify by trying to create again, this time with replace=True
    caplog_fixture.clear() # Clear previous log messages
    create_scalar_index(table, column_to_index, replace=True)
    assert f"Successfully created/updated scalar index on column '{column_to_index}'" in caplog_fixture.text
    assert f"in table '{table_name}' using type 'BTREE'" in caplog_fixture.text # Verify default type was used

    # Test creating an index with a specific type and ensure it's logged
    specific_index_type = "BITMAP" # Example, assuming BITMAP is a valid scalar index type
    if specific_index_type in ['BTREE', 'BITMAP', 'LABEL_LIST']: # Check if it's a valid type for the test
        caplog_fixture.clear()
        create_scalar_index(table, column_to_index, index_type=specific_index_type, replace=True) # Use the Literal type
        assert f"Successfully created/updated scalar index on column '{column_to_index}'" in caplog_fixture.text
        assert f"in table '{table_name}' using type '{specific_index_type}'" in caplog_fixture.text
    else:
        pytest.skip(f"Skipping specific index type test for {specific_index_type} as it's not in Literal")

    # Test attempting to create an index on a non-existent column (should fail gracefully)
    caplog_fixture.clear()
    with pytest.raises(Exception):
        create_scalar_index(table, "non_existent_column")
    assert "Failed to create scalar index" in caplog_fixture.text

def test_lancedb_cli_create_index(tmp_path, test_variant_data):
    db_path = tmp_path / "lancedb_cli_idx"
    table_name = "cli_idx_table"

    # Init and add data first
    subprocess.run([
        "python", "-m", "vcf_agent.cli", "init-lancedb",
        "--db_path", str(db_path), "--table_name", table_name
    ], check=True)
    embedding_str = ",".join(map(str, test_variant_data["embedding"]))
    add_command = [
        "python", "-m", "vcf_agent.cli", "add-variant",
        "--db_path", str(db_path), "--table_name", table_name,
        "--variant_id", test_variant_data["variant_id"],
        "--chrom", test_variant_data["chrom"],
        "--pos", str(test_variant_data["pos"]),
        "--ref", test_variant_data["ref"],
        "--alt", test_variant_data["alt"],
        "--embedding", embedding_str,
        "--clinical_significance", test_variant_data["clinical_significance"]
    ]
    subprocess.run(add_command, check=True)

    column_to_index = "pos"
    index_command = [
        "python", "-m", "vcf_agent.cli", "create-lancedb-index",
        "--db_path", str(db_path), "--table_name", table_name,
        "--column", column_to_index,
        "--replace" # Test with replace
    ]
    result = subprocess.run(index_command, capture_output=True, text=True, check=True)
    assert f"Index creation command processed for column '{column_to_index}'" in result.stdout

    # Verify by trying to query (Python API for simplicity here)
    db = get_db(db_path)
    table = db.open_table(table_name) # Open existing table
    results = table.search().where(f"{column_to_index} = {test_variant_data['pos']}").limit(1).to_pandas()
    assert not results.empty
    assert results.iloc[0]["variant_id"] == test_variant_data["variant_id"] 