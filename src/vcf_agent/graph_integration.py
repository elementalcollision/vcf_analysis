"""
Kuzu Graph Database Integration for VCF Agent

This module provides functions to interact with a Kuzu graph database,
including schema definition, data loading (variants, samples), and querying
using Cypher.
"""
import kuzu
import pandas as pd # Added for DataFrame conversion
import pyarrow as pa # Ensure pyarrow is imported
from typing import Dict, List, Any, Optional

# Default path for the Kuzu database
DEFAULT_KUZU_DB_PATH = "./kuzu_db"

# Global variable to hold the managed Kuzu connection
_kuzu_main_connection: Optional[kuzu.Connection] = None

def get_managed_kuzu_connection() -> kuzu.Connection:
    """
    Initializes and returns a managed Kuzu database connection.
    Ensures the connection is established once and the schema is created.

    Uses DEFAULT_KUZU_DB_PATH for the database location.

    Returns:
        A Kuzu database connection object.
    
    Raises:
        RuntimeError: If the Kuzu connection cannot be established or schema cannot be created.
    """
    global _kuzu_main_connection
    if _kuzu_main_connection is None:
        try:
            print(f"Initializing managed Kuzu connection to: {DEFAULT_KUZU_DB_PATH}")
            # Use the existing get_kuzu_db_connection which returns a new connection to a DB instance
            # This doesn't strictly create a singleton DB object across different calls to get_kuzu_db_connection
            # if db_path is different, but here we use a fixed path.
            # Kuzu's Database object itself handles the singleton nature for a given path.
            conn = get_kuzu_db_connection(db_path=DEFAULT_KUZU_DB_PATH)
            create_schema(conn) # Ensure schema exists on this connection
            _kuzu_main_connection = conn
            print("Managed Kuzu connection initialized and schema verified.")
        except Exception as e:
            # Log the error appropriately in a real application
            print(f"Failed to initialize managed Kuzu connection or create schema: {e}")
            raise RuntimeError(f"Kuzu setup failed: {e}") from e
    
    if _kuzu_main_connection is None: # Should not be reached if logic above is correct
        raise RuntimeError("Kuzu connection is unexpectedly None after initialization attempt.")
        
    return _kuzu_main_connection

def get_kuzu_db_connection(db_path: str = DEFAULT_KUZU_DB_PATH, read_only: bool = False) -> kuzu.Connection:
    """
    Initializes a Kuzu database at the given path and returns a connection.

    Args:
        db_path: Path to the Kuzu database directory.
        read_only: If True, opens the database in read-only mode.

    Returns:
        A Kuzu database connection object.
    """
    try:
        db = kuzu.Database(db_path)
        # TODO: Kuzu Python API for read_only mode needs to be confirmed.
        # As of recent versions, read_only might be set at Database level or not directly via Connection.
        # For now, we assume write access by default.
        # if read_only:
        #     conn = kuzu.Connection(db, access_mode=kuzu.AccessMode.READ_ONLY) # Example, check actual API
        # else:
        conn = kuzu.Connection(db)
        print(f"Successfully connected to Kuzu database at: {db_path}")
        return conn
    except Exception as e:
        print(f"Error connecting to Kuzu database at {db_path}: {e}")
        raise

def create_schema(conn: kuzu.Connection) -> None:
    """
    Creates the necessary node and relationship tables for VCF data if they don't exist.

    Schema:
    - Variant (Node): Represents a genetic variant.
        Properties: variant_id (STRING, PK), chrom (STRING), pos (INT64),
                    ref (STRING), alt (STRING), rs_id (STRING, optional)
    - Sample (Node): Represents a sample.
        Properties: sample_id (STRING, PK)
    - ObservedIn (Relationship): Links a Variant to a Sample.
        Properties: zygosity (STRING, e.g., HOM, HET)
    """
    schema_queries = [
        "CREATE NODE TABLE IF NOT EXISTS Variant(variant_id STRING, chrom STRING, pos INT64, ref STRING, alt STRING, rs_id STRING, PRIMARY KEY (variant_id))",
        "CREATE NODE TABLE IF NOT EXISTS Sample(sample_id STRING, PRIMARY KEY (sample_id))",
        "CREATE REL TABLE IF NOT EXISTS ObservedIn(FROM Variant TO Sample, zygosity STRING)"
    ]
    try:
        for query in schema_queries:
            print(f"Executing schema query: {query}")
            conn.execute(query)
        print("Kuzu schema created/verified successfully.")
    except Exception as e:
        print(f"Error creating Kuzu schema: {e}")
        raise

def add_variant(conn: kuzu.Connection, variant_data: Dict[str, Any]) -> None:
    """
    Adds a single variant node to the Kuzu database using parameterized queries.

    Args:
        conn: Active Kuzu connection.
        variant_data: Dictionary containing variant properties (variant_id, chrom, pos, ref, alt, rs_id).
                      Example: {'variant_id': 'chr1-123-A-G', 'chrom': '1', 'pos': 123,
                                'ref': 'A', 'alt': 'G', 'rs_id': 'rs12345'}
    """
    try:
        query = """
        CREATE (v:Variant {
            variant_id: $variant_id,
            chrom: $chrom,
            pos: $pos,
            ref: $ref,
            alt: $alt,
            rs_id: $rs_id
        })
        """
        params = {
            "variant_id": variant_data['variant_id'],
            "chrom": variant_data['chrom'],
            "pos": variant_data['pos'],
            "ref": variant_data['ref'],
            "alt": variant_data['alt'],
            "rs_id": variant_data.get('rs_id', '') # Ensure rs_id is provided, even if empty
        }
        conn.execute(query, parameters=params)
        print(f"Added variant: {variant_data['variant_id']}")
    except Exception as e:
        print(f"Error adding variant {variant_data.get('variant_id', 'UNKNOWN')}: {e}")
        raise

def add_sample(conn: kuzu.Connection, sample_data: Dict[str, Any]) -> None:
    """
    Adds a single sample node to the Kuzu database using parameterized queries.

    Args:
        conn: Active Kuzu connection.
        sample_data: Dictionary containing sample properties (sample_id).
                     Example: {'sample_id': 'SAMPLE_001'}
    """
    try:
        query = "CREATE (s:Sample {sample_id: $sample_id})"
        params = {"sample_id": sample_data['sample_id']}
        conn.execute(query, parameters=params)
        print(f"Added sample: {sample_data['sample_id']}")
    except Exception as e:
        print(f"Error adding sample {sample_data.get('sample_id', 'UNKNOWN')}: {e}")
        raise

def link_variant_to_sample(conn: kuzu.Connection, sample_id: str, variant_id: str, properties: Dict[str, Any]) -> None:
    """
    Creates an 'ObservedIn' relationship between a Sample and a Variant using parameterized queries.

    Args:
        conn: Active Kuzu connection.
        sample_id: ID of the Sample node.
        variant_id: ID of the Variant node.
        properties: Dictionary of properties for the relationship (e.g., {'zygosity': 'HET'}).
    """
    try:
        if 'zygosity' not in properties:
            raise ValueError("Missing 'zygosity' in properties for ObservedIn relationship.")

        query = """
        MATCH (v:Variant {variant_id: $v_id}), (s:Sample {sample_id: $s_id})
        CREATE (v)-[r:ObservedIn {zygosity: $zygosity}]->(s)
        """
        params = {
            "v_id": variant_id,
            "s_id": sample_id,
            "zygosity": properties['zygosity']
        }
        conn.execute(query, parameters=params)
        print(f"Linked sample {sample_id} to variant {variant_id} with properties: {properties}")
    except Exception as e:
        print(f"Error linking sample {sample_id} to variant {variant_id}: {e}")
        raise

def execute_query(conn: kuzu.Connection, cypher_query: str, params: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """
    Executes a given Cypher query using parameterized queries and returns the results
    as a list of dictionaries by converting to a Kuzu Arrow Table, then to a Pandas DataFrame.

    Args:
        conn: Active Kuzu connection.
        cypher_query: The Cypher query string to execute.
        params: Optional dictionary of parameters for the query. Keys should not include '$'.

    Returns:
        A list of dictionaries, where each dictionary represents a row in the result.
    """
    try:
        print(f"Executing Cypher query: {cypher_query} with params: {params}")
        query_result = conn.execute(cypher_query, parameters=params if params else {})

        # Convert Kuzu QueryResult directly to Pandas DataFrame
        df = query_result.get_as_df()
        
        if not isinstance(df, pd.DataFrame):
            # This case should ideally not be hit if get_as_df() works as expected.
            print(f"Warning: Expected Pandas DataFrame, but got {type(df)}.")
            raise TypeError(f"Conversion to DataFrame failed or returned unexpected type: {type(df)}.")

        results = df.to_dict(orient='records')
        print(f"Query returned {len(results)} results as dictionaries via DataFrame.")
        return results
    except Exception as e:
        print(f"Error executing Cypher query '{cypher_query}': {e}")
        raise

def get_variant_context(conn: kuzu.Connection, variant_ids: List[str]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Fetches graph context (e.g., samples, relationships) for a list of variant IDs from Kuzu.
    This function is intended to be called after a vector search in LanceDB provides relevant variant_ids.

    Args:
        conn: Active Kuzu connection.
        variant_ids: A list of variant IDs obtained from LanceDB vector search.

    Returns:
        A dictionary where keys are variant_ids and values are lists of dictionaries
        containing their graph context (e.g., associated samples and relationship properties).
        Example: {'rs123': [{'sample_id': 'S1', 'zygosity': 'HET'}]}}
    """
    if not variant_ids:
        return {}

    # Using UNWIND for efficient batch lookup if Kuzu supports it well for MATCH clauses.
    # Alternatively, loop and query one by one if UNWIND performance is not optimal for this pattern.
    query = """
    UNWIND $v_ids AS target_variant_id
    MATCH (v:Variant {variant_id: target_variant_id})-[o:ObservedIn]->(s:Sample)
    RETURN v.variant_id AS variant_id, s.sample_id AS sample_id, o.zygosity AS zygosity
    """
    params = {"v_ids": variant_ids}
    
    try:
        query_results = execute_query(conn, query, params=params)
        
        context_map: Dict[str, List[Dict[str, Any]]] = {v_id: [] for v_id in variant_ids}
        for row in query_results:
            v_id = row['variant_id']
            sample_info = {"sample_id": row['sample_id'], "zygosity": row['zygosity']}
            if v_id in context_map:
                context_map[v_id].append(sample_info)
        print(f"Fetched context for {len(context_map)} variants from Kuzu.")
        return context_map
    except Exception as e:
        print(f"Error fetching variant context from Kuzu: {e}")
        raise

if __name__ == '__main__':
    # Example Usage (illustrative, to be updated for get_variant_context)
    # db_connection = None
    # try:
    #     print("Attempting to connect to Kuzu and set up schema...")
    #     db_connection = get_kuzu_db_connection()
    #     create_schema(db_connection)
    #
    #     # Add some data (variant, sample, link)
    #     var_data1 = {'variant_id': 'rs123', 'chrom': '1', 'pos': 1000, 'ref': 'A', 'alt': 'T'}
    #     var_data2 = {'variant_id': 'rs456', 'chrom': '2', 'pos': 2000, 'ref': 'C', 'alt': 'G'}
    #     samp_data1 = {'sample_id': 'SAMPLE001'}
    #     samp_data2 = {'sample_id': 'SAMPLE002'}
    #     add_variant(db_connection, var_data1)
    #     add_variant(db_connection, var_data2)
    #     add_sample(db_connection, samp_data1)
    #     add_sample(db_connection, samp_data2)
    #     link_variant_to_sample(db_connection, 'SAMPLE001', 'rs123', {'zygosity': 'HET'})
    #     link_variant_to_sample(db_connection, 'SAMPLE001', 'rs456', {'zygosity': 'HOM'})
    #     link_variant_to_sample(db_connection, 'SAMPLE002', 'rs123', {'zygosity': 'HOM'})
    #
    #     print("\nQuerying all ObservedIn relationships...")
    #     all_links = execute_query(db_connection, "MATCH (v:Variant)-[o:ObservedIn]->(s:Sample) RETURN v.variant_id, s.sample_id, o.zygosity")
    #     for link in all_links:
    #         print(f"  Found Link: {link}")
    #
    #     # Simulate results from LanceDB vector search
    #     lancedb_results_variant_ids = ['rs123', 'rs456'] 
    #     print(f"\nFetching Kuzu context for variant IDs: {lancedb_results_variant_ids}")
    #     variant_contexts = get_variant_context(db_connection, lancedb_results_variant_ids)
    #     for v_id, context in variant_contexts.items():
    #         print(f"  Context for {v_id}: {context}")
    #
    # except NotImplementedError as nie:
    #     print(f"Feature not implemented: {nie}")
    # except Exception as e:
    #     print(f"An error occurred during Kuzu example usage: {e}")
    # finally:
    #     if db_connection:
    #         print("\nKuzu operations example finished.")
    pass 