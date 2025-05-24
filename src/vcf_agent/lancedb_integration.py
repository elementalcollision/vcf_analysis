# NOTE: If you see linter errors for lancedb imports, ensure your editor is using the correct virtual environment. The import paths are correct for LanceDB >=0.22.1.
# This module's logging strategy aims to align with the Auditable Security Framework,
# focusing on capturing key operations, parameters (with masking for sensitive data),
# and success/failure status.
from lancedb.pydantic import LanceModel, Vector
import lancedb
from typing import List, Optional, Annotated, Type, cast, Literal
import logging # Added for error logging
from lancedb.table import LanceTable # Ensure this import is present and correct
import threading # Added for locking
import re # Added for SQL literal masking

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
    variant_id: str
    embedding: Vector(1024)  # noqa: F821 - Suppress linter error for Vector(1024)
    chrom: str
    pos: int
    ref: str
    alt: str
    clinical_significance: Optional[str] = None

# LanceDB connection (local path for now)
def get_db(db_path: str = "./lancedb"):
    """Connects to a LanceDB database at the given path."""
    logger.info(f"Attempting to connect to LanceDB at {db_path}...")
    try:
        db = lancedb.connect(db_path)
        logger.info(f"Successfully connected to LanceDB at {db_path}.")
        return db
    except Exception as e:
        logger.error(f"Failed to connect to LanceDB at {db_path}: {e}")
        raise

def get_or_create_table(db, table_name: str = "variants"):
    """Opens an existing table or creates a new one if it doesn't exist."""
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

def add_variants(table: LanceTable, variants: List[dict]): # Added Type Hint for table
    """
    Adds a batch of variants to the specified LanceDB table.

    Args:
        table: The LanceDB table object.
        variants: A list of dictionaries, where each dictionary represents a variant conforming to the table schema.
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

def search_by_embedding(table: LanceTable, query_embedding, limit: int = 10, filter_sql: Optional[str] = None): # Added Type Hint
    """Searches the table by a query embedding, with optional limit and SQL filter."""
    # Avoid logging the query_embedding directly as it can be large, potentially sensitive (if inferable),
    # and adds verbosity. A summary of its structure is logged instead.
    embedding_summary = f"shape: {len(query_embedding)}x{len(query_embedding[0])}" if isinstance(query_embedding, list) and query_embedding and isinstance(query_embedding[0], list) else f"length: {len(query_embedding)}" if isinstance(query_embedding, list) else "opaque_embedding"
    masked_filter_sql = mask_sql_literals(filter_sql) # Mask SQL literals for logging
    logger.info(f"Attempting to search table '{table.name}' with embedding ({embedding_summary}), limit: {limit}, filter_sql: '{masked_filter_sql}'.")
    try:
        q = table.search(query_embedding, vector_column_name='embedding').limit(limit)
        if filter_sql:
            q = q.where(filter_sql)
        results = q.to_pandas()
        logger.info(f"Search in table '{table.name}' completed. Found {len(results)} results for limit: {limit}, filter_sql: '{masked_filter_sql}'.")
        return results
    except Exception as e:
        logger.error(f"Search by embedding failed in table '{table.name}' (limit: {limit}, filter_sql: '{masked_filter_sql}'): {e}")
        raise

def update_variant(table: LanceTable, variant_id: str, updates: dict): # Added Type Hint
    """Updates a specific variant record identified by variant_id.

    Args:
        table: The LanceDB table object.
        variant_id: The ID of the variant to update.
        updates: A dictionary containing the fields and new values to update.
                 Example: {"clinical_significance": "Benign", "pos": 12346}
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

def delete_variants(table: LanceTable, filter_sql: str): # Added Type Hint
    """Deletes variant records from the table based on a SQL filter expression.

    Args:
        table: The LanceDB table object.
        filter_sql: A SQL-like filter string to identify records for deletion.
                    Example: "variant_id = 'rs123'" or "chrom = \'1\' AND clinical_significance = \'Pathogenic\'"
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

def create_scalar_index(table: LanceTable, column_name: str, index_type: Optional[LanceScalarIndexType] = None, replace: bool = False, **kwargs): # Added Type Hint & use Literal
    """Creates a scalar index on a specified column to improve filter performance.

    Args:
        table: The LanceDB table object.
        column_name: The name of the scalar column to index (e.g., 'chrom', 'pos', 'clinical_significance').
        index_type: Optional. The type of index to create. Must be one of 'BTREE', 'BITMAP', 'LABEL_LIST'.
                    If None, LanceDB's default for scalar columns (BTREE) will be used.
        replace: Optional. If True, replaces an existing index on the column.
        **kwargs: Additional keyword arguments to pass to table.create_scalar_index().
                  (e.g., for specific index configurations if supported for scalar indexes)
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
    """Filters variants in the table based on a SQL metadata filter.

    Args:
        table: The LanceDB table object.
        filter_sql: A SQL-like filter string to apply to metadata columns.
                    Example: "chrom = '1' AND pos > 1000 AND clinical_significance = 'Pathogenic'"
        select_columns: Optional list of column names to return. If None, all columns are returned.
        limit: Optional maximum number of records to return.

    Returns:
        A pandas DataFrame containing the filtered results.
    """
    # TODO: Security Review: filter_sql is logged directly. If filter_sql could contain sensitive PII
    # (e.g., filtering on a direct identifier that is PII), it should be masked or summarized
    # in production logs according to the data handling policies. -> Addressed by masking.
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