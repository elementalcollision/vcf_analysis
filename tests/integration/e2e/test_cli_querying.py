import pytest
from pathlib import Path
import numpy as np # For mock embeddings
import json # For parsing CLI output if it's JSON

# Utility functions from other modules if needed (e.g., for setting up DB)
from vcf_agent.lancedb_integration import get_db as get_lancedb_connection, get_or_create_table as get_or_create_lancedb_table, add_variants, Variant
from cyvcf2 import VCF # For parsing VCF to get data for setup

# Test Scenario 3.1 (Part 1): Querying LanceDB by embedding via CLI
def test_search_lancedb_embedding_cli(test_dbs, sample_vcf_file_small, vcf_agent_cli_runner):
    """
    E2E Test for 'search-embedding' CLI command.
    Verifies searching by embedding vector, with and without SQL filters.
    """
    lancedb_path = test_dbs["lancedb_path"]
    # Use a unique table name to avoid conflicts if tests run in parallel or don't clean up perfectly
    lancedb_table_name = "variants_cli_search_test"

    # 1. Setup: Populate LanceDB with known variants and embeddings
    variants_to_add = []
    variant_embeddings_map = {}
    try:
        vcf_parser = VCF(sample_vcf_file_small)
        for i, record in enumerate(vcf_parser):
            variant_id = f"{record.CHROM}-{record.POS}-{record.REF}-{record.ALT[0]}"
            mock_embedding = np.zeros(1024, dtype=np.float32) # Default embedding
            if i == 0: # Make the first variant's embedding distinct (all 1s)
                mock_embedding = np.ones(1024, dtype=np.float32)
            
            variant_embeddings_map[variant_id] = mock_embedding.tolist()
            
            variant_data = {
                "variant_id": variant_id,
                "embedding": mock_embedding.tolist(),
                "chrom": record.CHROM,
                "pos": record.POS,
                "ref": record.REF,
                "alt": record.ALT[0],
                "clinical_significance": f"SignificanceCLI_{i}"
            }
            variants_to_add.append(variant_data)
        vcf_parser.close()
    except Exception as e:
        pytest.fail(f"Failed to parse VCF or prepare variant data for CLI search test: {e}")

    assert len(variants_to_add) == 2
    target_variant_id_to_search = variants_to_add[0]["variant_id"]
    target_embedding_list = variant_embeddings_map[target_variant_id_to_search]
    target_embedding_cli_str = ",".join(map(str, target_embedding_list))

    try:
        db_conn = get_lancedb_connection(lancedb_path)
        if lancedb_table_name in db_conn.table_names():
            db_conn.drop_table(lancedb_table_name)
        # Ensure the schema used by get_or_create_lancedb_table matches what cli.py's add-variant would use.
        # The default is Variant, which is correct.
        lance_table = get_or_create_lancedb_table(db_conn, table_name=lancedb_table_name)
        add_variants(lance_table, variants_to_add) # Populate the table
        print(f"LanceDB table '{lancedb_table_name}' populated for CLI search test.")
    except Exception as e:
        pytest.fail(f"Failed to set up LanceDB for CLI search test: {e}")

    # 2. Test 'search-embedding' CLI command (no filter)
    cli_args_search = [
        "search-embedding",
        "--db_path", lancedb_path,
        "--table_name", lancedb_table_name,
        "--embedding", target_embedding_cli_str,
        "--limit", "1"
    ]
    result_search = vcf_agent_cli_runner(cli_args_search)
    print(f"CLI search stdout: {result_search.stdout}")
    print(f"CLI search stderr: {result_search.stderr}")
    assert result_search.returncode == 0, f"CLI 'search-embedding' failed. Stderr: {result_search.stderr}"
    # The CLI prints the DataFrame string representation. We need to check if the variant_id is in it.
    assert target_variant_id_to_search in result_search.stdout, \
        f"Expected variant {target_variant_id_to_search} not found in CLI search output."

    # 3. Test 'search-embedding' CLI command (with a matching filter)
    filter_sql_match = f"clinical_significance = 'SignificanceCLI_0'"
    cli_args_search_filter_match = [
        "search-embedding",
        "--db_path", lancedb_path,
        "--table_name", lancedb_table_name,
        "--embedding", target_embedding_cli_str,
        "--limit", "1",
        "--filter_sql", filter_sql_match
    ]
    result_search_filter_match = vcf_agent_cli_runner(cli_args_search_filter_match)
    print(f"CLI search (filter match) stdout: {result_search_filter_match.stdout}")
    assert result_search_filter_match.returncode == 0, f"CLI 'search-embedding' with filter failed. Stderr: {result_search_filter_match.stderr}"
    assert target_variant_id_to_search in result_search_filter_match.stdout, \
        f"Expected variant {target_variant_id_to_search} not found in CLI search output with matching filter."

    # 4. Test 'search-embedding' CLI command (with a non-matching filter)
    filter_sql_no_match = f"clinical_significance = 'NonExistentSignificance'"
    cli_args_search_filter_no_match = [
        "search-embedding",
        "--db_path", lancedb_path,
        "--table_name", lancedb_table_name,
        "--embedding", target_embedding_cli_str,
        "--limit", "1",
        "--filter_sql", filter_sql_no_match
    ]
    result_search_filter_no_match = vcf_agent_cli_runner(cli_args_search_filter_no_match)
    print(f"CLI search (filter no match) stdout: {result_search_filter_no_match.stdout}")
    assert result_search_filter_no_match.returncode == 0, f"CLI 'search-embedding' with non-matching filter failed. Stderr: {result_search_filter_no_match.stderr}"
    # Expecting empty DataFrame or a message indicating no results. The current CLI prints the DF, which might show headers only.
    # A robust check is that the target_variant_id is NOT in the output IF the output format for empty is just headers or minimal.
    # If it prints "Empty DataFrame", that's also a good sign.
    assert target_variant_id_to_search not in result_search_filter_no_match.stdout, \
        f"Variant {target_variant_id_to_search} unexpectedly found in CLI search output with non-matching filter."
    # A more specific check could be for pandas empty dataframe representation if known
    # For example: `assert "Empty DataFrame" in result_search_filter_no_match.stdout` or check row count if printed. 

# Test Scenario 3.1 (Part 2): Querying LanceDB by metadata filter via CLI
def test_filter_lancedb_metadata_cli(test_dbs, sample_vcf_file_small, vcf_agent_cli_runner):
    """
    E2E Test for 'filter-lancedb' CLI command.
    Verifies filtering by metadata, column selection, and limit.
    """
    lancedb_path = test_dbs["lancedb_path"]
    lancedb_table_name = "variants_cli_filter_test" # Distinct table

    # 1. Setup: Populate LanceDB (embeddings are not critical for this filter test)
    variants_to_add = []
    try:
        vcf_parser = VCF(sample_vcf_file_small)
        for i, record in enumerate(vcf_parser):
            variant_id = f"{record.CHROM}-{record.POS}-{record.REF}-{record.ALT[0]}"
            mock_embedding = np.random.rand(1024).astype(np.float32).tolist() # Still need embedding for schema
            variant_data = {
                "variant_id": variant_id,
                "embedding": mock_embedding,
                "chrom": record.CHROM,
                "pos": record.POS,
                "ref": record.REF,
                "alt": record.ALT[0],
                "clinical_significance": f"SignificanceFilterCLI_{i}"
            }
            variants_to_add.append(variant_data)
        vcf_parser.close()
    except Exception as e:
        pytest.fail(f"Failed to parse VCF or prepare variant data for CLI filter test: {e}")
    
    target_variant_1 = variants_to_add[0]
    target_variant_2 = variants_to_add[1]

    try:
        db_conn = get_lancedb_connection(lancedb_path)
        if lancedb_table_name in db_conn.table_names():
            db_conn.drop_table(lancedb_table_name)
        lance_table = get_or_create_lancedb_table(db_conn, table_name=lancedb_table_name)
        add_variants(lance_table, variants_to_add)
        print(f"LanceDB table '{lancedb_table_name}' populated for CLI filter test.")
    except Exception as e:
        pytest.fail(f"Failed to set up LanceDB for CLI filter test: {e}")

    # 2. Test filter matching one variant by ID
    filter_by_id = f"variant_id = '{target_variant_1['variant_id']}'"
    cli_args_filter_id = [
        "filter-lancedb",
        "--db_path", lancedb_path,
        "--table_name", lancedb_table_name,
        "--filter_sql", filter_by_id
    ]
    result_filter_id = vcf_agent_cli_runner(cli_args_filter_id)
    print(f"CLI filter by ID stdout: {result_filter_id.stdout}")
    assert result_filter_id.returncode == 0, f"CLI 'filter-lancedb' by ID failed. Stderr: {result_filter_id.stderr}"
    assert target_variant_1['variant_id'] in result_filter_id.stdout
    assert target_variant_2['variant_id'] not in result_filter_id.stdout

    # 3. Test filter matching by other metadata (e.g., chrom and pos for target_variant_2)
    # pos is an int, so no quotes in SQL filter
    filter_by_meta = f"chrom = '{target_variant_2['chrom']}' AND pos = {target_variant_2['pos']}"
    cli_args_filter_meta = [
        "filter-lancedb",
        "--db_path", lancedb_path,
        "--table_name", lancedb_table_name,
        "--filter_sql", filter_by_meta
    ]
    result_filter_meta = vcf_agent_cli_runner(cli_args_filter_meta)
    print(f"CLI filter by metadata stdout: {result_filter_meta.stdout}")
    assert result_filter_meta.returncode == 0, f"CLI 'filter-lancedb' by metadata failed. Stderr: {result_filter_meta.stderr}"
    assert target_variant_2['variant_id'] in result_filter_meta.stdout
    assert target_variant_1['variant_id'] not in result_filter_meta.stdout

    # 4. Test filter that should return no results
    filter_no_match = f"clinical_significance = 'ThisShouldNotExist'"
    cli_args_filter_no_match = [
        "filter-lancedb",
        "--db_path", lancedb_path,
        "--table_name", lancedb_table_name,
        "--filter_sql", filter_no_match
    ]
    result_filter_no_match = vcf_agent_cli_runner(cli_args_filter_no_match)
    print(f"CLI filter no match stdout: {result_filter_no_match.stdout}")
    assert result_filter_no_match.returncode == 0, f"CLI 'filter-lancedb' with no match failed. Stderr: {result_filter_no_match.stderr}"
    # The CLI output for no results usually includes "Found 0 variants matching filter:" and an empty DataFrame representation.
    assert "Found 0 variants" in result_filter_no_match.stdout or "Empty DataFrame" in result_filter_no_match.stdout
    assert target_variant_1['variant_id'] not in result_filter_no_match.stdout
    assert target_variant_2['variant_id'] not in result_filter_no_match.stdout

    # 5. Test --select_columns and --limit
    filter_select_limit = f"chrom = '{target_variant_1['chrom']}'" # Should match both variants if not for limit
    select_cols = "variant_id,pos"
    cli_args_select_limit = [
        "filter-lancedb",
        "--db_path", lancedb_path,
        "--table_name", lancedb_table_name,
        "--filter_sql", filter_select_limit,
        "--select_columns", select_cols,
        "--limit", "1"
    ]
    result_select_limit = vcf_agent_cli_runner(cli_args_select_limit)
    print(f"CLI filter with select/limit stdout: {result_select_limit.stdout}")
    assert result_select_limit.returncode == 0, f"CLI 'filter-lancedb' with select/limit failed. Stderr: {result_select_limit.stderr}"
    # Check that only selected columns are present (variant_id, pos) and not others (e.g., ref, alt, clinical_significance)
    assert "variant_id" in result_select_limit.stdout
    assert "pos" in result_select_limit.stdout
    assert "chrom" not in result_select_limit.stdout # 'chrom' was in filter but not select_columns for output DF
    assert "ref" not in result_select_limit.stdout
    assert "alt" not in result_select_limit.stdout
    assert "clinical_significance" not in result_select_limit.stdout
    assert "embedding" not in result_select_limit.stdout # Embedding should definitely not be there
    # Check that only 1 record is returned due to limit=1. 
    # This is harder to check directly from string output if multiple records could match the variant_id string. 
    # A more robust check would be to count lines or parse the output if it were structured (e.g., JSON).
    # For now, we assume if it contains one of the IDs, it's likely limited correctly.
    # A simple check: ensure it contains one ID but not the other if they would have different string representations
    # (which they do, as variant_id is unique).
    # This implicitly tests the limit if we know both variants match the filter_select_limit.
    found_target1 = target_variant_1['variant_id'] in result_select_limit.stdout
    found_target2 = target_variant_2['variant_id'] in result_select_limit.stdout
    assert (found_target1 and not found_target2) or (not found_target1 and found_target2), \
        "Limit=1 should result in exactly one of the two possible matching variants being returned."
    assert "Found 1 variants" in result_select_limit.stdout # CLI usually prints the count 