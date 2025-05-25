import pytest
from pathlib import Path
import kuzu
import numpy as np # For mock embeddings
from cyvcf2 import VCF # For VCF parsing
from typing import Union, Type, List # Import Union, Type, and List
from _pytest.python_api import RaisesContext # For explicit typing of excinfo
import gc # For garbage collection

# Assuming your lancedb_integration and graph_integration (kuzu) modules
# have functions to connect and query, which will be used for assertions.
from vcf_agent.lancedb_integration import get_db as get_lancedb_connection, get_or_create_table, add_variants, Variant
from vcf_agent.vcf_utils import populate_kuzu_from_vcf
from vcf_agent.graph_integration import create_schema as create_kuzu_schema # Renamed for clarity

# Define a helper function to count nodes in Kuzu
def count_kuzu_nodes(kuzu_db_path: str, node_label: str) -> int:
    """
    Counts nodes with a specific label in Kuzu using context managers.
    Fails the test via pytest.fail() on error.
    """
    try:
        with kuzu.Database(kuzu_db_path) as db:
            with kuzu.Connection(db) as conn:
                query_result_union = conn.execute(f"MATCH (n:{node_label}) RETURN count(n)")
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
                    count_tuple = single_query_result.get_next()
                    if count_tuple and len(count_tuple) > 0:
                        count = int(count_tuple[0])
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
    except Exception as e:
        pytest.fail(f"Error counting Kuzu nodes for label {node_label} at {kuzu_db_path}: {e}")
        return -1 

def test_successful_vcf_ingestion_api(test_dbs, sample_vcf_file_small):
    """
    E2E Test Scenario 1.1 (API): Successful ingestion of a valid VCF file via Python APIs.
    Verifies LanceDB storage and Kuzu node/relationship creation.
    """
    lancedb_path = test_dbs["lancedb_path"]
    kuzu_path = test_dbs["kuzu_path"]

    # --- Kuzu Population ---
    kuzu_conn = None
    kuzu_db_obj = None
    try:
        Path(kuzu_path).mkdir(parents=True, exist_ok=True)
        kuzu_db_obj = kuzu.Database(kuzu_path)
        kuzu_conn = kuzu.Connection(kuzu_db_obj)
        create_kuzu_schema(kuzu_conn) # Ensure schema exists

        kuzu_counts = populate_kuzu_from_vcf(kuzu_conn, str(sample_vcf_file_small))
        assert kuzu_counts.get("variants", 0) == 2
        assert kuzu_counts.get("samples", 0) == 1
        assert kuzu_counts.get("links", 0) == 2
    except Exception as e:
        pytest.fail(f"Kuzu population failed during API ingestion test: {e}")
    finally:
        if kuzu_conn: kuzu_conn.close()
        if kuzu_db_obj: kuzu_db_obj.close() # Close DB object if it was created

    # --- LanceDB Population (Simplified example of adding one variant) ---
    lancedb_conn = None
    try:
        Path(lancedb_path).mkdir(parents=True, exist_ok=True)
        lancedb_conn = get_lancedb_connection(lancedb_path)
        table = get_or_create_table(lancedb_conn, "variants")
        
        # Create some mock Variant data (as Pydantic models, then dump to dict for add_variants)
        # This part would typically parse the VCF and create these.
        # For this combined test, we assume Kuzu population handled VCF details.
        # We just add a couple of mock variants to LanceDB to verify its part.
        mock_embeddings = [[0.1] * 1024, [0.2] * 1024]
        variants_to_add = [
            Variant(variant_id="chr1-100-A-G", chrom="1", pos=100, ref="A", alt="G", embedding=mock_embeddings[0]).model_dump(),
            Variant(variant_id="chr1-200-C-T", chrom="1", pos=200, ref="C", alt="T", embedding=mock_embeddings[1]).model_dump(),
        ]
        add_variants(table, variants_to_add)
        assert table.count_rows() == 2, "LanceDB table should have 2 variants after population."
    except Exception as e:
        pytest.fail(f"LanceDB population failed during API ingestion test: {e}")
    finally:
        # LanceDB connection object (db) does not have a close method for file-based DBs.
        # Cleanup is done by the test_dbs fixture by removing the directory.
        pass

    # Final verification: query Kuzu to ensure data is there as expected
    assert count_kuzu_nodes(kuzu_path, "Variant") == 2, "Kuzu DB should have 2 Variant nodes."
    assert count_kuzu_nodes(kuzu_path, "Sample") == 1, "Kuzu DB should have 1 Sample node."
    gc.collect() # Force GC after all Kuzu operations

def test_ingestion_non_existent_vcf_api(test_dbs):
    """
    E2E Test Scenario 1.2 (API): Attempt ingestion of a non-existent VCF file via API.
    Verifies appropriate error (e.g., FileNotFoundError).
    """
    kuzu_path = test_dbs["kuzu_path"]
    non_existent_vcf = "/path/to/surely/non_existent_file.vcf"

    kuzu_conn = None
    try:
        Path(kuzu_path).mkdir(parents=True, exist_ok=True)
        kuzu_db_obj = kuzu.Database(kuzu_path)
        kuzu_conn = kuzu.Connection(kuzu_db_obj)
        create_kuzu_schema(kuzu_conn)
        
        # Define the expected exception types
        expected_exceptions_non_existent: Type[Union[FileNotFoundError, RuntimeError]] = (FileNotFoundError, RuntimeError) # type: ignore
        with pytest.raises(expected_exceptions_non_existent) as excinfo: # type: ignore
            populate_kuzu_from_vcf(kuzu_conn, non_existent_vcf)
        assert "not found" in str(excinfo.value).lower() or "no such file" in str(excinfo.value).lower()
    except Exception as e:
        pytest.fail(f"Test setup for non-existent VCF API failed unexpectedly: {e}")
    finally:
        if kuzu_conn:
            kuzu_conn.close()

def test_ingestion_invalid_vcf_format_api(test_dbs, invalid_vcf_file):
    """
    E2E Test Scenario 1.3 (API): Attempt ingestion of an invalid format VCF file via API.
    Verifies appropriate parsing error.
    """
    kuzu_path = test_dbs["kuzu_path"]
    kuzu_conn = None
    try:
        Path(kuzu_path).mkdir(parents=True, exist_ok=True)
        kuzu_db_obj = kuzu.Database(kuzu_path)
        kuzu_conn = kuzu.Connection(kuzu_db_obj)
        create_kuzu_schema(kuzu_conn)

        # cyvcf2.VCF often raises OSError for badly malformed VCFs or issues opening/reading them.
        # It can also raise ValueError for some format issues. kuzu.RuntimeException is also possible if error propagates.
        expected_exceptions_invalid_vcf: Type[Union[OSError, ValueError, RuntimeError]] = (OSError, ValueError, RuntimeError) # type: ignore
        with pytest.raises(expected_exceptions_invalid_vcf) as excinfo: # type: ignore
            populate_kuzu_from_vcf(kuzu_conn, str(invalid_vcf_file))
        
        # Check for common error substrings related to VCF parsing issues
        error_str = str(excinfo.value).lower()
        assert "parse" in error_str or "invalid" in error_str or "format" in error_str or "malformed" in error_str or "error reading vcf" in error_str, \
            f"Expected a VCF parsing related error, but got: {error_str}"
    except Exception as e:
        pytest.fail(f"Test setup for invalid VCF API failed unexpectedly: {e}")
    finally:
        if kuzu_conn:
            kuzu_conn.close()

# Remove original CLI-based tests or comment them out if they are to be reactivated later
# when ingest-vcf CLI command is implemented. 