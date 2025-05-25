import pytest
from pathlib import Path
import kuzu
import numpy as np # For mock embeddings
from cyvcf2 import VCF # For VCF parsing
from typing import Union, Type, List # Import Union, Type, and List
from _pytest.python_api import RaisesContext # For explicit typing of excinfo
import gc # For garbage collection

# Imports from the VCF agent codebase
from vcf_agent.vcf_utils import populate_kuzu_from_vcf
from vcf_agent.graph_integration import get_kuzu_db_connection, create_schema # For setup if needed, or direct kuzu.Connection

# Helper function to count nodes/rels (can be shared or adapted from test_cli_ingestion.py)
def count_kuzu_entities(kuzu_db_path: str, cypher_query: str) -> int:
    """Counts entities in Kuzu based on a Cypher query."""
    try:
        with kuzu.Database(kuzu_db_path) as db:
            with kuzu.Connection(db) as conn:
                query_result_union = conn.execute(cypher_query)
                single_query_result = None
                if isinstance(query_result_union, list):
                    if not query_result_union:
                        return 0
                    if isinstance(query_result_union[0], kuzu.query_result.QueryResult):
                        single_query_result = query_result_union[0]
                    else:
                        pytest.fail(f"Kuzu conn.execute returned list with non-QueryResult: {query_result_union[0]}")
                        return -2
                elif isinstance(query_result_union, kuzu.query_result.QueryResult):
                    single_query_result = query_result_union
                else:
                    pytest.fail(f"Kuzu conn.execute returned unexpected type: {type(query_result_union)}")
                    return -3

                count = 0
                if single_query_result and single_query_result.has_next():
                    result_tuple = single_query_result.get_next()
                    if result_tuple and len(result_tuple) > 0:
                        count = result_tuple[0]
                # Explicitly delete all QueryResult objects
                if single_query_result:
                    del single_query_result
                if query_result_union:
                    if isinstance(query_result_union, list):
                        for item in query_result_union:
                            if item: del item
                    else:
                        del query_result_union
                gc.collect() # Force GC after deletion
                return count
    except ImportError:
        pytest.fail("Kuzu library not installed. It is required for these E2E tests.")
    except Exception as e:
        pytest.fail(f"Error querying Kuzu in count_kuzu_entities (query: '{cypher_query}'): {e}")
    return -1

@pytest.mark.integration
def test_successful_vcf_to_kuzu_api(test_dbs, sample_vcf_file_small):
    """
    E2E Test Scenario 1.4 (Part 1 - Kuzu): Successful VCF data population into Kuzu via Python API.
    Verifies node and relationship creation in Kuzu.
    """
    kuzu_path = test_dbs["kuzu_path"]
    db = None
    conn = None
    print(f"Starting test_successful_vcf_to_kuzu_api with kuzu_path: {kuzu_path}")
    try:
        Path(kuzu_path).mkdir(parents=True, exist_ok=True) 
        print(f"Attempting to create Kuzu Database at {kuzu_path}")
        db = kuzu.Database(kuzu_path)
        print("Kuzu Database object created.")
        conn = kuzu.Connection(db)
        print("Kuzu Connection object created.")
        
        print("Calling create_schema...")
        create_schema(conn) 
        print("create_schema completed.")

        print("Calling populate_kuzu_from_vcf...")
        counts = populate_kuzu_from_vcf(conn, str(sample_vcf_file_small))
        print(f"populate_kuzu_from_vcf returned: {counts}")
        
        assert counts.get("variants", 0) == 2, "Incorrect number of variants reported as added by API." 
        assert counts.get("samples", 0) == 1, "Incorrect number of samples reported as added by API."  
        assert counts.get("links", 0) == 2, "Incorrect number of links reported as added by API."      
        
        print("Assertions on counts passed.")
        # Force garbage collection before closing Kuzu objects, as a safeguard
        print("Forcing garbage collection before finally block...")
        gc.collect()
        print("Original Kuzu operations completed successfully in try block.")

    except Exception as e:
        pytest.fail(f"Kuzu operation failed: {e}")
    finally:
        print("In finally block of test_successful_vcf_to_kuzu_api")
        if conn:
            print("Closing Kuzu connection...")
            conn.close()
            print("Kuzu connection closed.")
        else:
            print("Kuzu connection was None in finally block.")
        if db:
            print("Closing Kuzu database...")
            db.close()
            print("Kuzu database closed.")
        else:
            print("Kuzu database was None in finally block.")
    print("Test function test_successful_vcf_to_kuzu_api finished.")

    print("Performing final Kuzu node count verification...")
    assert count_kuzu_entities(kuzu_path, "MATCH (n:Variant) RETURN count(n)") == 2, "Kuzu DB should have 2 Variant nodes after population."
    assert count_kuzu_entities(kuzu_path, "MATCH (n:Sample) RETURN count(n)") == 1, "Kuzu DB should have 1 Sample node after population."
    print("Final Kuzu node count verification completed.")
    gc.collect() # Force GC after all Kuzu operations

# Test for LanceDB population via Python API (Scenario 1.4 - Part 2)
from vcf_agent.lancedb_integration import get_db as get_lancedb_connection, get_or_create_table as get_or_create_lancedb_table, add_variants
from cyvcf2 import VCF # For parsing VCF in the test
import numpy as np # For generating mock embeddings

def test_successful_vcf_to_lancedb_api(test_dbs, sample_vcf_file_small):
    """
    E2E Test Scenario 1.4 (Part 2 - LanceDB): Successful VCF data population into LanceDB via Python API.
    This test manually parses the VCF, creates mock embeddings, and calls add_variants.
    Verifies data in LanceDB.
    """
    lancedb_path = test_dbs["lancedb_path"]
    lancedb_table_name = "variants_api_test" # Use a distinct table name for this test

    # 1. Prepare Variant Data (Parse VCF and create mock embeddings)
    variants_to_add = []
    try:
        vcf_parser = VCF(sample_vcf_file_small)
        for record in vcf_parser:
            variant_id = f"{record.CHROM}-{record.POS}-{record.REF}-{record.ALT[0]}"
            # Create a mock 1024-dimensional embedding
            mock_embedding = np.random.rand(1024).astype(np.float32).tolist()
            variant_data = {
                "variant_id": variant_id,
                "embedding": mock_embedding,
                "chrom": record.CHROM,
                "pos": record.POS,
                "ref": record.REF,
                "alt": record.ALT[0],
                # 'clinical_significance' is optional in the Pydantic model, so not strictly needed here
            }
            variants_to_add.append(variant_data)
        vcf_parser.close()
    except Exception as e:
        pytest.fail(f"Failed to parse VCF or prepare variant data for LanceDB: {e}")
    
    assert len(variants_to_add) == 2, "Should have prepared 2 variants from the sample VCF."

    # 2. Get LanceDB connection and table
    try:
        db_conn = get_lancedb_connection(lancedb_path)
        # Ensure the table is created fresh for this test to avoid interference
        if lancedb_table_name in db_conn.table_names():
            db_conn.drop_table(lancedb_table_name)
        lance_table = get_or_create_lancedb_table(db_conn, table_name=lancedb_table_name)
    except Exception as e:
        pytest.fail(f"Failed to set up LanceDB for API test: {e}")

    # 3. Call the Python API to add variants to LanceDB
    try:
        add_variants(lance_table, variants_to_add)
        print(f"add_variants API called successfully for {len(variants_to_add)} variants.")
    except Exception as e:
        pytest.fail(f"API call to add_variants failed: {e}")

    # 4. Verify LanceDB content
    try:
        variants_in_db_df = lance_table.to_pandas()
        assert len(variants_in_db_df) == 2, "Incorrect number of variants in LanceDB table."

        # Verify a couple of fields for the first variant from the VCF
        expected_variant_id_1 = f"chr1-100-A-G" # Based on sample_vcf_file_small
        db_variant1_df = variants_in_db_df[variants_in_db_df["variant_id"] == expected_variant_id_1]
        assert not db_variant1_df.empty, f"Variant {expected_variant_id_1} not found in LanceDB."
        assert db_variant1_df.iloc[0]["chrom"] == "chr1"
        assert db_variant1_df.iloc[0]["pos"] == 100
        assert "embedding" in db_variant1_df.columns
        assert len(db_variant1_df.iloc[0]["embedding"]) == 1024

        print(f"LanceDB verification successful. Found {len(variants_in_db_df)} variants.")

    except Exception as e:
        pytest.fail(f"Failed to verify LanceDB content after API call: {e}")

# Scenario 3.2: Querying and Data Retrieval (LanceDB - Python API)
from vcf_agent.lancedb_integration import search_by_embedding

def test_search_lancedb_by_embedding_api(test_dbs, sample_vcf_file_small):
    """
    E2E Test Scenario 3.2 (LanceDB): Successful search by embedding via Python API.
    Verifies that `search_by_embedding` returns the correct variant.
    """
    lancedb_path = test_dbs["lancedb_path"]
    lancedb_table_name = "variants_search_api_test" # Distinct table for this test

    # 1. Setup: Populate LanceDB with known variants and embeddings
    variants_to_add = []
    variant_embeddings = {}
    try:
        vcf_parser = VCF(sample_vcf_file_small)
        for i, record in enumerate(vcf_parser):
            variant_id = f"{record.CHROM}-{record.POS}-{record.REF}-{record.ALT[0]}"
            # Create a predictable mock embedding for searching
            # For simplicity, make one embedding very distinct for easy search
            mock_embedding = np.zeros(1024, dtype=np.float32)
            if i == 0: # Make the first variant's embedding unique (e.g., all 1s)
                mock_embedding = np.ones(1024, dtype=np.float32)
            
            variant_embeddings[variant_id] = mock_embedding.tolist() # Store for searching
            
            variant_data = {
                "variant_id": variant_id,
                "embedding": mock_embedding.tolist(),
                "chrom": record.CHROM,
                "pos": record.POS,
                "ref": record.REF,
                "alt": record.ALT[0],
                "clinical_significance": f"Significance_{i}" # Add some metadata for potential filter tests
            }
            variants_to_add.append(variant_data)
        vcf_parser.close()
    except Exception as e:
        pytest.fail(f"Failed to parse VCF or prepare variant data for LanceDB search test: {e}")

    assert len(variants_to_add) == 2
    target_variant_id_to_search = variants_to_add[0]["variant_id"]
    target_embedding_to_search = variant_embeddings[target_variant_id_to_search]

    try:
        db_conn = get_lancedb_connection(lancedb_path)
        if lancedb_table_name in db_conn.table_names():
            db_conn.drop_table(lancedb_table_name)
        lance_table = get_or_create_lancedb_table(db_conn, table_name=lancedb_table_name)
        add_variants(lance_table, variants_to_add)
    except Exception as e:
        pytest.fail(f"Failed to set up LanceDB for search API test: {e}")

    # 2. Call search_by_embedding API
    try:
        # Search for the first variant's embedding
        search_results_df = search_by_embedding(lance_table, target_embedding_to_search, limit=1)
        
        assert not search_results_df.empty, "Search returned no results."
        assert len(search_results_df) == 1, "Search should return exactly one result with limit=1 and unique embedding."
        assert search_results_df.iloc[0]["variant_id"] == target_variant_id_to_search, \
            f"Search returned wrong variant. Expected {target_variant_id_to_search}, got {search_results_df.iloc[0]['variant_id']}"
        print(f"LanceDB search_by_embedding for {target_variant_id_to_search} successful.")

        # Test with a filter that should match
        filter_sql_match = f"clinical_significance = 'Significance_0'"
        search_results_filtered_df = search_by_embedding(
            lance_table, 
            target_embedding_to_search, 
            limit=1, 
            filter_sql=filter_sql_match
        )
        assert not search_results_filtered_df.empty, "Filtered search returned no results."
        assert search_results_filtered_df.iloc[0]["variant_id"] == target_variant_id_to_search
        print(f"LanceDB search_by_embedding with filter '{filter_sql_match}' successful.")

        # Test with a filter that should NOT match the top result for this embedding
        filter_sql_no_match = f"clinical_significance = 'NonExistentSignificance'"
        search_results_no_match_df = search_by_embedding(
            lance_table, 
            target_embedding_to_search, 
            limit=1, 
            filter_sql=filter_sql_no_match
        )
        assert search_results_no_match_df.empty, \
            f"Filtered search should have returned no results for filter '{filter_sql_no_match}', but got {len(search_results_no_match_df)}."
        print(f"LanceDB search_by_embedding with filter '{filter_sql_no_match}' correctly returned no results.")

    except Exception as e:
        pytest.fail(f"API call to search_by_embedding failed: {e}")

# Scenario 3.2: Querying and Data Retrieval (Kuzu - Python API)
from vcf_agent.graph_integration import get_variant_context, create_schema # create_schema already imported but good to be explicit

def test_get_kuzu_variant_context_api(test_dbs, sample_vcf_file_small):
    """
    E2E Test Scenario 3.2 (Kuzu): Successful retrieval of variant context from Kuzu via Python API.
    Verifies that `get_variant_context` (with a passed connection) returns correct sample and zygosity information.
    """
    kuzu_path = test_dbs["kuzu_path"]

    # 1. Setup: Populate Kuzu with data from sample_vcf_file_small
    Path(kuzu_path).mkdir(parents=True, exist_ok=True)
    db_setup = kuzu.Database(kuzu_path)
    conn_setup = kuzu.Connection(db_setup)
    try:
        create_schema(conn_setup) # Ensure schema exists
        populate_kuzu_from_vcf(conn_setup, sample_vcf_file_small) 
        print(f"Kuzu database at {kuzu_path} populated for context test.")
    except Exception as e:
        pytest.fail(f"Failed to set up Kuzu for variant context test: {e}")
    finally:
        conn_setup.close()
        db_setup.close()

    # 2. Define variant IDs to query
    variant_ids_to_query = [
        "chr1-100-A-G", 
        "chr1-200-C-T"
    ]
    non_existent_variant_id = "chrX-999-N-N"
    all_query_ids = variant_ids_to_query + [non_existent_variant_id]

    # 3. Call get_variant_context API with a test-specific connection
    retrieved_contexts = None
    db_query = None
    conn_query = None
    try:
        db_query = kuzu.Database(kuzu_path) # Connect to the DB populated in step 1
        conn_query = kuzu.Connection(db_query)
        
        # Call the refactored get_variant_context, passing our test connection
        retrieved_contexts = get_variant_context(conn_query, all_query_ids)

    except Exception as e:
        pytest.fail(f"Call to get_variant_context failed: {e}")
    finally:
        if conn_query:
            conn_query.close()
        if db_query:
            db_query.close()

    # 4. Assert the retrieved context (same assertions as before)
    assert retrieved_contexts is not None, "Context retrieval failed."
    
    rs1_id = "chr1-100-A-G"
    assert rs1_id in retrieved_contexts, f"{rs1_id} not found in retrieved contexts."
    rs1_context = retrieved_contexts[rs1_id]
    assert len(rs1_context) == 1, f"Expected 1 sample context for {rs1_id}, got {len(rs1_context)}."
    assert rs1_context[0]["sample_id"] == "sample1"
    assert rs1_context[0]["zygosity"] == "HET", f"Incorrect zygosity for {rs1_id}. Expected HET, got {rs1_context[0]['zygosity']}."

    rs2_id = "chr1-200-C-T"
    assert rs2_id in retrieved_contexts, f"{rs2_id} not found in retrieved contexts."
    rs2_context = retrieved_contexts[rs2_id]
    assert len(rs2_context) == 1, f"Expected 1 sample context for {rs2_id}, got {len(rs2_context)}."
    assert rs2_context[0]["sample_id"] == "sample1"
    assert rs2_context[0]["zygosity"] == "HOM_ALT", f"Incorrect zygosity for {rs2_id}. Expected HOM_ALT, got {rs2_context[0]['zygosity']}."

    assert non_existent_variant_id in retrieved_contexts, f"Non-existent ID {non_existent_variant_id} should be in keys."
    assert len(retrieved_contexts[non_existent_variant_id]) == 0, \
        f"Expected no context for non-existent variant {non_existent_variant_id}, got {len(retrieved_contexts[non_existent_variant_id])}."

    print(f"Kuzu get_variant_context API test successful for IDs: {all_query_ids}")

# TODO: Add test for LanceDB population via Python API (Scenario 1.4 - Part 2)
# This will require:
# 1. A VCF parsing function (e.g., from vcf_utils or a new one).
# 2. An embedding generation step (mocked or real if available as an API).
# 3. Formatting data to match lancedb_integration.Variant schema.
# 4. Calling lancedb_integration.add_variants. 