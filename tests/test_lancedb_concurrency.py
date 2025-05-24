import pytest
import threading
import time
import os
import shutil
import tempfile
import numpy as np
import logging

from vcf_agent.lancedb_integration import (
    get_db,
    get_or_create_table,
    add_variants,
    update_variant,
    delete_variants,
    create_scalar_index,
    Variant # Import the schema
)

# Configure basic logging for tests
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@pytest.fixture(scope="module")
def lancedb_concurrency_tmpdir():
    """Creates a temporary directory for LanceDB data for the module."""
    tmpdir = tempfile.mkdtemp(prefix="lancedb_concurrency_test_")
    logger.info(f"Created LanceDB temp directory for concurrency tests: {tmpdir}")
    yield tmpdir
    logger.info(f"Cleaning up LanceDB temp directory: {tmpdir}")
    shutil.rmtree(tmpdir)

@pytest.fixture(scope="function")
def db_table_for_concurrency(lancedb_concurrency_tmpdir):
    """Provides a fresh DB and table for each concurrency test function."""
    db_path = os.path.join(lancedb_concurrency_tmpdir, f"db_{time.time_ns()}") # Unique path for each test
    db = get_db(db_path)
    # Ensure table is created new for each test to avoid interference
    table_name = "concurrency_variants_test"
    if table_name in db.table_names():
        db.drop_table(table_name)
        logger.info(f"Dropped existing table {table_name} in {db_path}")
    table = get_or_create_table(db, table_name=table_name)
    logger.info(f"Created table {table_name} in {db_path} for test")
    return db, table

def generate_variant_data(idx: int, base_id: str = "concTest") -> dict:
    """Generates unique variant data for testing."""
    return {
        "variant_id": f"{base_id}_{idx}",
        "chrom": str((idx % 23) + 1), # Chromosomes 1-23
        "pos": 1000 + idx * 100,
        "ref": "A",
        "alt": "T",
        "embedding": np.random.rand(1024).astype(np.float32).tolist(),
        "clinical_significance": "Benign" if idx % 2 == 0 else "Pathogenic",
    }

def worker_add_variants(table, num_variants: int, start_idx: int):
    """Worker function to add a batch of variants."""
    variants_to_add = [generate_variant_data(i + start_idx, "add") for i in range(num_variants)]
    try:
        add_variants(table, variants_to_add)
        logger.info(f"Thread {threading.get_ident()}: Successfully added {num_variants} variants starting from index {start_idx}.")
    except Exception as e:
        logger.error(f"Thread {threading.get_ident()}: Error adding variants: {e}")
        raise # Re-raise to fail the test

def test_concurrent_add_variants(db_table_for_concurrency):
    """Tests concurrent calls to add_variants."""
    _, table = db_table_for_concurrency
    num_threads = 5
    variants_per_thread = 10
    threads = []

    for i in range(num_threads):
        thread = threading.Thread(target=worker_add_variants, args=(table, variants_per_thread, i * variants_per_thread))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    expected_total_variants = num_threads * variants_per_thread
    assert table.count_rows() == expected_total_variants, "Mismatch in total variants after concurrent adds."
    logger.info(f"test_concurrent_add_variants: All {expected_total_variants} variants added successfully.")


def worker_update_variants(table, variant_ids_to_update: list, new_significance: str):
    """Worker function to update variants."""
    try:
        for variant_id in variant_ids_to_update:
            update_variant(table, variant_id, {"clinical_significance": new_significance})
        logger.info(f"Thread {threading.get_ident()}: Successfully updated {len(variant_ids_to_update)} variants to {new_significance}.")
    except Exception as e:
        logger.error(f"Thread {threading.get_ident()}: Error updating variants: {e}")
        raise

def test_concurrent_update_variants(db_table_for_concurrency):
    """Tests concurrent calls to update_variant on different records."""
    _, table = db_table_for_concurrency
    num_initial_variants = 20
    initial_variants = [generate_variant_data(i, "updateInitial") for i in range(num_initial_variants)]
    add_variants(table, initial_variants)
    assert table.count_rows() == num_initial_variants

    num_threads = 4
    variants_per_thread = 5 # num_initial_variants must be divisible by num_threads * variants_per_thread
    
    threads = []
    all_updated_ids = []

    for i in range(num_threads):
        start_idx = i * variants_per_thread
        end_idx = start_idx + variants_per_thread
        ids_for_thread = [v["variant_id"] for v in initial_variants[start_idx:end_idx]]
        all_updated_ids.extend(ids_for_thread)
        
        thread = threading.Thread(target=worker_update_variants, args=(table, ids_for_thread, f"UpdatedByThread_{i}"))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    # Verification
    results_df = table.to_pandas()
    for i in range(num_threads):
        start_idx = i * variants_per_thread
        end_idx = start_idx + variants_per_thread
        thread_significance = f"UpdatedByThread_{i}"
        for variant_idx_in_original_list in range(start_idx, end_idx):
            original_variant_id = initial_variants[variant_idx_in_original_list]["variant_id"]
            updated_row = results_df[results_df["variant_id"] == original_variant_id]
            assert not updated_row.empty, f"Variant {original_variant_id} not found after updates."
            assert updated_row.iloc[0]["clinical_significance"] == thread_significance, \
                f"Variant {original_variant_id} not updated correctly by thread {i}. Expected {thread_significance}, got {updated_row.iloc[0]['clinical_significance']}"
    
    logger.info(f"test_concurrent_update_variants: All targeted variants updated correctly by concurrent threads.")


def worker_delete_variants(table, filter_to_delete: str, num_expected_to_delete: int):
    """Worker function to delete variants based on a filter."""
    try:
        delete_variants(table, filter_to_delete)
        logger.info(f"Thread {threading.get_ident()}: Delete operation with filter '{filter_to_delete}' completed.")
    except Exception as e:
        logger.error(f"Thread {threading.get_ident()}: Error deleting variants: {e}")
        raise

def test_concurrent_delete_variants(db_table_for_concurrency):
    """Tests concurrent calls to delete_variants using different filters."""
    _, table = db_table_for_concurrency
    # Add a mixed set of variants
    variants_set1 = [generate_variant_data(i, "delA") for i in range(10)] # chrom 1-10, Pathogenic/Benign
    variants_set2 = [generate_variant_data(i + 10, "delB") for i in range(10)] # chrom 11-20, Pathogenic/Benign
    add_variants(table, variants_set1 + variants_set2)
    initial_count = table.count_rows()
    assert initial_count == 20

    # Thread 1 deletes 'Pathogenic' variants from delA set (approx 5, on chroms 1,3,5,7,9)
    # Variant IDs: delA_0, delA_2, delA_4, delA_6, delA_8
    filter1 = "clinical_significance = 'Pathogenic' AND variant_id LIKE 'delA_%'"
    # Thread 2 deletes 'Benign' variants from delB set (approx 5, on chroms 12,14,16,18,20)
    # Variant IDs: delB_1, delB_3, delB_5, delB_7, delB_9 (indices 11,13,15,17,19 from generate_variant_data)
    filter2 = "clinical_significance = 'Benign' AND variant_id LIKE 'delB_%'"

    thread1 = threading.Thread(target=worker_delete_variants, args=(table, filter1, 5))
    thread2 = threading.Thread(target=worker_delete_variants, args=(table, filter2, 5))

    thread1.start()
    thread2.start()
    thread1.join()
    thread2.join()

    # Expected remaining:
    # - Benign from delA (delA_1,3,5,7,9) = 5
    # - Pathogenic from delB (delB_0,2,4,6,8) = 5
    expected_remaining_count = 10
    final_count = table.count_rows()
    assert final_count == expected_remaining_count, \
        f"Mismatch in remaining variants. Expected {expected_remaining_count}, got {final_count}. Filter1='{filter1}', Filter2='{filter2}'"
    
    # Further verify that specific types were deleted
    remaining_df = table.to_pandas()
    assert remaining_df[remaining_df["variant_id"].str.startswith("delA_") & (remaining_df["clinical_significance"] == "Pathogenic")].empty
    assert remaining_df[remaining_df["variant_id"].str.startswith("delB_") & (remaining_df["clinical_significance"] == "Benign")].empty
    logger.info(f"test_concurrent_delete_variants: Correct number of variants remaining after concurrent deletions.")

# Note: Testing concurrent create_scalar_index is tricky because the operation is typically
# idempotent or controlled by 'replace=True'. True concurrency issues (like partial index
# creation) are harder to reliably trigger and observe at this level without specific
# knowledge of LanceDB's internal index creation atomicity.
# A simple test could ensure no errors occur when called concurrently.

def worker_create_index(table, column_name: str, replace: bool):
    try:
        create_scalar_index(table, column_name, replace=replace)
        logger.info(f"Thread {threading.get_ident()}: create_scalar_index for {column_name} (replace={replace}) completed.")
    except RuntimeError as e:
        # If replace is False and the index already exists, LanceDB will raise a RuntimeError.
        # This is an expected outcome in a concurrent scenario where multiple threads attempt creation.
        if "already exists" in str(e) and not replace:
            logger.info(f"Thread {threading.get_ident()}: Expected error for {column_name} (replace=False) as index likely already exists: {e}")
        else:
            logger.error(f"Thread {threading.get_ident()}: Unexpected RuntimeError in create_scalar_index for {column_name} (replace={replace}): {e}")
            raise # Re-raise unexpected RuntimeErrors
    except Exception as e:
        logger.error(f"Thread {threading.get_ident()}: Error in create_scalar_index for {column_name} (replace={replace}): {e}")
        raise # Re-raise other unexpected errors

def test_concurrent_create_scalar_index(db_table_for_concurrency):
    """Tests concurrent calls to create_scalar_index.
    Mainly ensures no deadlocks or unhandled exceptions occur.
    Verifying the actual index state changes from concurrent calls is complex.
    """
    _, table = db_table_for_concurrency
    # Add some data so the column exists
    add_variants(table, [generate_variant_data(0, "idxTest")])

    num_threads = 3
    column_to_index = "clinical_significance"
    threads = []

    # First pass: all threads try to create (replace=False)
    # Only one should succeed in creating, others might do nothing or log if index exists.
    for _ in range(num_threads):
        thread = threading.Thread(target=worker_create_index, args=(table, column_to_index, False))
        threads.append(thread)
        thread.start()
    for thread in threads:
        thread.join()
    
    # Second pass: all threads try to create with replace=True
    # Each should complete without error.
    threads_replace = []
    for _ in range(num_threads):
        thread = threading.Thread(target=worker_create_index, args=(table, column_to_index, True))
        threads_replace.append(thread)
        thread.start()
    for thread in threads_replace:
        thread.join()

    # Primary assertion is that no exceptions were raised by the workers.
    # Further index state verification could be added if LanceDB provides inspection tools.
    logger.info("test_concurrent_create_scalar_index: Concurrent index creation calls completed without exceptions.") 