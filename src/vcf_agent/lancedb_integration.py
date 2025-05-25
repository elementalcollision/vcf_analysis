# NOTE: If you see linter errors for lancedb imports, ensure your editor is using the correct virtual environment. The import paths are correct for LanceDB >=0.22.1.
# This module's logging strategy aims to align with the Auditable Security Framework,
# focusing on capturing key operations, parameters (with masking for sensitive data),
# and success/failure status.
from lancedb.pydantic import LanceModel, Vector
import lancedb
from typing import List, Optional, Annotated, Type, cast, Literal, Dict, Tuple
import logging # Added for error logging
from lancedb.table import LanceTable # Ensure this import is present and correct
import threading # Added for locking
import re # Added for SQL literal masking
import pandas as pd # Ensure pandas is imported for DataFrame operations

# Import Kuzu integration functions
from . import graph_integration

logger = logging.getLogger(__name__) # Added for error logging

# Global lock for LanceDB write operations to ensure thread safety with local file-based DB
# This lock serializes write-related operations (add, update, delete, create_index)
# to prevent potential race conditions or corruption when accessing the underlying LanceDB files
# from multiple threads within a single VCF Agent process.
# For multi-process or distributed access, a more robust distributed locking mechanism or a
# database designed for high concurrency (like LanceDB Cloud) would be necessary.
lancedb_write_lock = threading.Lock()

# Helper function to mask string literals in SQL-like filter strings for safer logging
def mask_sql_literals(sql_string: Optional[str]) -> Optional[str]:
    """
    Masks string literals (single or double quoted) in a SQL-like string.
    Example: "name = \'John Doe\' AND age > 30" -> "name = \'{MASKED_STRING}\' AND age > 30"
    """
    if sql_string is None:
        return None
    # Mask single-quoted strings
    masked_sql = re.sub(r"'(?:\\\'|[^'])*'", "'{MASKED_STRING}'", sql_string)
    # Mask double-quoted strings (though less common in SQL, good to cover)
    masked_sql = re.sub(r'"(?:\\"|[^"])*"', '"{MASKED_STRING}"', masked_sql)
    return masked_sql

# LanceDB/Pydantic best practice: use Annotated for vector fields to specify dimension
# See: https://lancedb.github.io/lancedb/python/pydantic/
class Variant(LanceModel):
    """
    Pydantic model representing a variant record in LanceDB.

    Attributes:
        variant_id (str): Unique identifier for the variant (e.g., 'chr1-123-A-G').
        embedding (Vector): 1024-dimensional embedding vector for the variant.
        chrom (str): Chromosome (e.g., '1', 'X').
        pos (int): Genomic position.
        ref (str): Reference allele.
        alt (str): Alternate allele.
        clinical_significance (Optional[str]): Clinical significance annotation (e.g., 'Pathogenic').
    """
    variant_id: str
    embedding: Vector(1024)  # type: ignore[valid-type]
    chrom: str
    pos: int
    ref: str
    alt: str
    clinical_significance: Optional[str] = None

# LanceDB connection (local path for now)
def get_db(db_path: str = "./lancedb"):
    """
    Connects to a LanceDB database at the given path.

    Args:
        db_path (str): Path to the LanceDB directory. Defaults to './lancedb'.

    Returns:
        lancedb.DBConnection: LanceDB database connection object.

    Raises:
        Exception: If connection fails.

    Example:
        >>> db = get_db("./lancedb")
    """
    logger.info(f"Attempting to connect to LanceDB at {db_path}...")
    try:
        db = lancedb.connect(db_path)
        logger.info(f"Successfully connected to LanceDB at {db_path}.")
        return db
    except Exception as e:
        logger.error(f"Failed to connect to LanceDB at {db_path}: {e}")
        raise

def get_or_create_table(db, table_name: str = "variants"):
    """
    Opens an existing LanceDB table or creates a new one if it doesn't exist.

    Args:
        db (lancedb.DBConnection): LanceDB database connection object.
        table_name (str): Name of the table to open or create. Defaults to 'variants'.

    Returns:
        LanceTable: The opened or newly created LanceDB table.

    Raises:
        Exception: If table creation or opening fails.

    Example:
        >>> db = get_db()
        >>> table = get_or_create_table(db, "variants")
    """
    logger.info(f"Attempting to get or create table '{table_name}'...")
    try:
        if table_name in db.table_names():
            table = db.open_table(table_name)
            logger.info(f"Opened existing table: '{table_name}'.")
        else:
            table = db.create_table(table_name, schema=Variant, mode="overwrite")
            logger.info(f"Created new table: '{table_name}' with schema {Variant.__name__}.") # Log schema name
        return table
    except Exception as e:
        logger.error(f"Failed to get or create table '{table_name}': {e}")
        raise

def add_variants(table: LanceTable, variants: List[dict]):
    """
    Adds a batch of variant records to the specified LanceDB table.

    Args:
        table (LanceTable): The LanceDB table object.
        variants (List[dict]): List of variant dictionaries matching the Variant schema.

    Returns:
        None

    Raises:
        Exception: If the add operation fails.

    Example:
        >>> add_variants(table, [{"variant_id": "chr1-123-A-G", ...}])
    """
    num_variants = len(variants)
    logger.info(f"Attempting to add {num_variants} variants to table '{table.name}'...")
    if not variants:
        logger.warning(f"add_variants called for table '{table.name}' with an empty list. No data will be added.")
        return

    try:
        with lancedb_write_lock:
            table.add(variants)
        logger.info(f"Successfully added {num_variants} variants to table '{table.name}'.")
    except Exception as e:
        logger.error(f"Failed to add batch of {num_variants} variants to table '{table.name}': {e}")
        # Depending on requirements, you might want to:
        # 1. Raise the exception to let the caller handle it.
        # 2. Implement partial success/retry logic if applicable.
        # 3. Log problematic records if possible.
        raise

def search_by_embedding(table: LanceTable, query_embedding, limit: int = 10, filter_sql: Optional[str] = None) -> pd.DataFrame:
    """
    Searches the LanceDB table by a query embedding, with optional limit and SQL filter.
    Enriches results with context from the Kuzu graph database.

    Args:
        table (LanceTable): The LanceDB table object.
        query_embedding (List[float] or np.ndarray): Embedding vector to search by.
        limit (int): Maximum number of results to return. Defaults to 10.
        filter_sql (Optional[str]): Optional SQL-like filter string for metadata columns.

    Returns:
        pd.DataFrame: DataFrame of search results, enriched with Kuzu context in 'kuzu_observed_samples'.

    Raises:
        Exception: If the search or enrichment fails.

    Example:
        >>> df = search_by_embedding(table, [0.1]*1024, limit=5)
    """
    embedding_summary = f"shape: {len(query_embedding)}x{len(query_embedding[0])}" if isinstance(query_embedding, list) and query_embedding and isinstance(query_embedding[0], list) else f"length: {len(query_embedding)}" if isinstance(query_embedding, list) else "opaque_embedding"
    masked_filter_sql = mask_sql_literals(filter_sql) # Mask SQL literals for logging
    logger.info(f"Attempting to search table '{table.name}' with embedding ({embedding_summary}), limit: {limit}, filter_sql: '{masked_filter_sql}'.")
    try:
        q = table.search(query_embedding, vector_column_name='embedding').limit(limit)
        if filter_sql:
            q = q.where(filter_sql)
        results_df = q.to_df()
        logger.info(f"Search in table '{table.name}' completed. Found {len(results_df)} results for limit: {limit}, filter_sql: '{masked_filter_sql}'.")

        # Always add the 'kuzu_observed_samples' column, even if empty
        results_df['kuzu_observed_samples'] = [[] for _ in range(len(results_df))]

        # Enrich with Kuzu graph context if results are found
        if not results_df.empty and 'variant_id' in results_df.columns:
            variant_ids_to_enrich = results_df['variant_id'].tolist()
            if variant_ids_to_enrich:
                try:
                    logger.info(f"Fetching Kuzu context for {len(variant_ids_to_enrich)} variant IDs.")
                    kuzu_conn = graph_integration.get_managed_kuzu_connection()
                    if kuzu_conn:
                        variant_contexts = graph_integration.get_variant_context(kuzu_conn, variant_ids_to_enrich)
                        if isinstance(variant_contexts, list):
                            from typing import Dict, List, Tuple
                            context_map: Dict[str, List[Tuple[str, str]]] = {}
                            for row in variant_contexts:
                                v_id = row['variant_id']
                                tup: Tuple[str, str] = (row['sample_id'], row['zygosity'])
                                if v_id not in context_map:
                                    context_map[v_id] = []
                                context_map[v_id].append(tup)
                            variant_contexts = context_map
                        results_df['kuzu_observed_samples'] = results_df['variant_id'].apply(
                            lambda vid: variant_contexts.get(vid, [])
                        )
                        logger.info("Successfully enriched LanceDB results with Kuzu context.")
                    else:
                        logger.warning("Kuzu connection not available, skipping enrichment.")
                        # Already set to empty lists above
                except Exception as kuzu_e:
                    logger.error(f"Error during Kuzu context enrichment: {kuzu_e}")
                    import pandas as pd
                    results_df['kuzu_observed_samples'] = [pd.NA] * len(results_df)
            # else: already set to empty lists above
        elif 'variant_id' not in results_df.columns:
            logger.warning("'variant_id' column not found in LanceDB results, cannot enrich with Kuzu context.")
            # Already set to empty lists above
        # else: results_df is empty, already set to empty lists above

        return results_df
    except Exception as e:
        logger.error(f"Search by embedding failed in table '{table.name}' (limit: {limit}, filter_sql: '{masked_filter_sql}'): {e}")
        raise

def update_variant(table: LanceTable, variant_id: str, updates: dict):
    """
    Updates a specific variant record identified by variant_id.

    Args:
        table (LanceTable): The LanceDB table object.
        variant_id (str): The ID of the variant to update.
        updates (dict): Dictionary of fields and new values to update.

    Returns:
        None

    Raises:
        Exception: If the update operation fails.

    Example:
        >>> update_variant(table, "chr1-123-A-G", {"clinical_significance": "Benign"})
    """
    update_keys = list(updates.keys()) if updates else []
    logger.info(f"Attempting to update variant '{variant_id}' in table '{table.name}' with updates for keys: {update_keys}.")

    if not updates:
        logger.warning(f"update_variant called for variant '{variant_id}' in table '{table.name}' with no updates. Skipping.")
        return

    try:
        with lancedb_write_lock:
            # Ensure variant_id is properly quoted if it's a string type in the schema.
            # LanceDB's SQL-like .where() clauses expect string literals to be single-quoted.
            where_clause = f"variant_id = '{variant_id}'"
            table.update(where=where_clause, values=updates)
        logger.info(f"Successfully updated variant '{variant_id}' in table '{table.name}' with keys: {update_keys}.")
    except Exception as e:
        # Avoid logging the full 'updates' dict here directly in case of sensitive values.
        logger.error(f"Failed to update variant '{variant_id}' in table '{table.name}' with update keys {update_keys}: {e}")
        raise

def delete_variants(table: LanceTable, filter_sql: str):
    """
    Deletes variant records from the table based on a SQL filter expression.

    Args:
        table (LanceTable): The LanceDB table object.
        filter_sql (str): SQL-like filter string to identify records for deletion.

    Returns:
        None

    Raises:
        ValueError: If filter_sql is empty.
        Exception: If the delete operation fails.

    Example:
        >>> delete_variants(table, "variant_id = 'chr1-123-A-G'")
    """
    masked_filter_sql = mask_sql_literals(filter_sql) # Mask SQL literals for logging
    logger.info(f"Attempting to delete variants from table '{table.name}' matching filter: '{masked_filter_sql}'.")
    # TODO: Security Review: filter_sql is logged directly. If filter_sql could contain sensitive PII
    # (e.g., filtering on a direct identifier that is PII), it should be masked or summarized
    # in production logs according to the data handling policies. -> Addressed by masking.
    if not filter_sql: # Original filter_sql used for logic, masked_filter_sql for logging.
        logger.error(f"delete_variants called for table '{table.name}' with an empty filter_sql. No records will be deleted.")
        # Potentially raise an error or return early depending on desired behavior
        raise ValueError("filter_sql cannot be empty for delete_variants operation.")

    try:
        with lancedb_write_lock:
            table.delete(filter_sql)
        logger.info(f"Successfully deleted variants from table '{table.name}' matching filter: '{masked_filter_sql}'.")
    except Exception as e:
        logger.error(f"Failed to delete variants from table '{table.name}' using filter '{masked_filter_sql}': {e}")
        raise

LanceScalarIndexType = Literal['BTREE', 'BITMAP', 'LABEL_LIST'] # Define Literal for index types

def create_scalar_index(table: LanceTable, column_name: str, index_type: Optional[LanceScalarIndexType] = None, replace: bool = False, **kwargs):
    """
    Creates a scalar index on a specified column to improve filter performance.

    Args:
        table (LanceTable): The LanceDB table object.
        column_name (str): Name of the scalar column to index (e.g., 'chrom', 'pos').
        index_type (Optional[str]): Type of index ('BTREE', 'BITMAP', 'LABEL_LIST'). Defaults to 'BTREE'.
        replace (bool): If True, replaces an existing index. Defaults to False.
        **kwargs: Additional keyword arguments for index creation.

    Returns:
        None

    Raises:
        Exception: If index creation fails.

    Example:
        >>> create_scalar_index(table, "chrom", index_type="BTREE")
    """
    kwargs_summary = list(kwargs.keys()) # Log only keys of kwargs
    # Log the index_type that will actually be used.
    actual_index_type_to_log = index_type if index_type is not None else 'BTREE'
    logger.info(f"Attempting to create scalar index on column '{column_name}' in table '{table.name}' with type '{actual_index_type_to_log}', replace={replace}, kwargs_keys={kwargs_summary}.")
    try:
        with lancedb_write_lock:
            if index_type is not None:
                table.create_scalar_index(column_name, index_type=index_type, replace=replace, **kwargs)
            else:
                # Explicitly use 'BTREE' which is a valid Literal member
                table.create_scalar_index(column_name, index_type='BTREE', replace=replace, **kwargs)
        
        logger.info(f"Successfully created/updated scalar index on column '{column_name}' in table '{table.name}' using type '{actual_index_type_to_log}', replace={replace}.")
    except Exception as e:
        logger.error(f"Failed to create scalar index on column '{column_name}' in table '{table.name}' (type: {actual_index_type_to_log}, replace: {replace}): {e}")
        raise

def filter_variants_by_metadata(table: LanceTable, filter_sql: str, select_columns: Optional[List[str]] = None, limit: Optional[int] = None):
    """
    Filters variants in the table based on a SQL metadata filter.

    Args:
        table (LanceTable): The LanceDB table object.
        filter_sql (str): SQL-like filter string to apply to metadata columns.
        select_columns (Optional[List[str]]): List of columns to return. If None, all columns are returned.
        limit (Optional[int]): Maximum number of records to return.

    Returns:
        pd.DataFrame: DataFrame of filtered results.

    Raises:
        ValueError: If filter_sql is empty.
        Exception: If the filter operation fails.

    Example:
        >>> df = filter_variants_by_metadata(table, "chrom = '1' AND pos > 1000")
    """
    masked_filter_sql = mask_sql_literals(filter_sql) # Mask SQL literals for logging
    logger.info(f"Attempting to filter variants from table '{table.name}' using filter: '{masked_filter_sql}', select_columns: {select_columns}, limit: {limit}.")
    if not filter_sql: # Original filter_sql used for logic, masked_filter_sql for logging.
        logger.error(f"filter_variants_by_metadata called for table '{table.name}' with an empty filter_sql. This would return all rows (if not for safety checks).")
        raise ValueError("filter_sql cannot be empty for metadata filtering operations.")

    try:
        # logger.debug(f"Table object in filter_variants_by_metadata: type={{type(table)}}, dir={{dir(table)}}") # DEBUG LINE # Kept for potential future debugging
        # Hypothesis: Start with search(None) or search() for pure scalar queries
        # If search() requires a vector, this will fail. 
        # If search(None) works, or search() alone initiates a full scan, this is the path.
        # Try calling search() without any arguments, assuming it initiates a full table query context
        query_builder = table.search() 
        
        if select_columns:
            query_builder = query_builder.select(select_columns)
        
        query_builder = query_builder.where(filter_sql)

        if limit is not None:
            query_builder = query_builder.limit(limit)
        
        results_df = query_builder.to_pandas()
        logger.info(f"Successfully filtered variants from table '{table.name}' using filter '{masked_filter_sql}'. Found {len(results_df)} records (select_columns: {select_columns}, limit: {limit}).")
        return results_df
    except Exception as e:
        logger.error(f"Failed to filter variants from table '{table.name}' using filter '{masked_filter_sql}' (select_columns: {select_columns}, limit: {limit}): {e}")
        raise

# Example usage (conceptual, actual embedding generation and Pydantic model would be needed)
# from pydantic import BaseModel
# ... existing code ... 