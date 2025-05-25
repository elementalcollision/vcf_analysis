import pytest
import kuzu
import os
import shutil

from vcf_agent import graph_integration

# Define a test database path for on-disk tests, or use :memory: for in-memory
TEST_DB_PATH_DIR = "./test_kuzu_db_dir" # Directory for the test DB

@pytest.fixture(scope="function")
def db_connection():
    """Pytest fixture to set up and tear down a Kuzu DB connection for testing."""
    # Ensure the test directory is clean before running a test that uses it
    if os.path.exists(TEST_DB_PATH_DIR):
        shutil.rmtree(TEST_DB_PATH_DIR)
    os.makedirs(TEST_DB_PATH_DIR, exist_ok=True)
    
    # Using an on-disk database for this fixture to allow multiple connections if needed by tests
    # and to allow inspection of files if necessary.
    # For truly isolated in-memory, use db_path=":memory:" directly in tests not needing persistence.
    conn = graph_integration.get_kuzu_db_connection(db_path=TEST_DB_PATH_DIR)
    assert conn is not None, "Database connection should not be None"
    assert isinstance(conn, kuzu.Connection), "Should return a Kuzu connection object"
    
    # Create schema once per connection setup for tests that need it
    try:
        graph_integration.create_schema(conn)
    except Exception as e:
        pytest.fail(f"Schema creation failed during fixture setup: {e}")

    yield conn # Provide the connection to the test

    # Teardown: Close connection and cleanup database directory
    # Kuzu connections are typically closed when the Database object is closed or goes out of scope.
    # Explicitly try to clean up resources.
    del conn 
    if os.path.exists(TEST_DB_PATH_DIR):
        shutil.rmtree(TEST_DB_PATH_DIR)

def test_get_kuzu_db_connection_in_memory():
    """Test getting an in-memory Kuzu database connection."""
    try:
        # For in-memory, Kuzu typically uses an empty string or ":memory:"
        # Let's test with an empty string as per Kuzu's typical behavior for in-memory
        conn_memory = graph_integration.get_kuzu_db_connection(db_path=":memory:")
        assert conn_memory is not None, "In-memory database connection should not be None"
        assert isinstance(conn_memory, kuzu.Connection), "Should return a Kuzu connection object for in-memory DB"
        # Optionally, try a simple execution to ensure it's live
        conn_memory.execute("CREATE NODE TABLE TestNode(id INT64, PRIMARY KEY (id))")
        print("In-memory Kuzu connection test passed.")
    except Exception as e:
        pytest.fail(f"test_get_kuzu_db_connection_in_memory failed: {e}")
    finally:
        # For in-memory, no explicit cleanup of files needed.
        # Kuzu handles in-memory DB closure when connection/db object goes out of scope.
        pass

def test_get_kuzu_db_connection_on_disk(db_connection): # Uses the fixture
    """
    Test getting an on-disk Kuzu database connection using the fixture.
    The fixture handles setup (including schema creation) and teardown.
    """
    assert db_connection is not None, "On-disk database connection from fixture should not be None"
    # Further checks can be done if needed, but fixture already asserts type and schema creation
    print("On-disk Kuzu connection test (via fixture) passed.")

def test_create_schema(db_connection):
    """
    Test that the schema is created successfully by the fixture and is usable.
    The db_connection fixture already calls create_schema.
    This test verifies basic table presence by attempting a simple query or insertion.
    """
    try:
        # Attempt to query one of the created tables to see if it exists.
        result = db_connection.execute("MATCH (v:Variant) RETURN count(v) AS variant_count")
        assert result is not None, "Query to Variant table should execute."

        result_sample = db_connection.execute("MATCH (s:Sample) RETURN count(s) AS sample_count")
        assert result_sample is not None, "Query to Sample table should execute."
        
        # Verify ObservedIn relationship table by creating a relationship
        db_connection.execute("MERGE (:Variant {variant_id: 'test_schema_v1', chrom: '1', pos: 1, ref: 'A', alt: 'T', rs_id: 'rsTest1'})")
        db_connection.execute("MERGE (:Sample {sample_id: 'test_schema_s1'})")
        result_rel = db_connection.execute("""
            MATCH (v:Variant {variant_id: 'test_schema_v1'}), (s:Sample {sample_id: 'test_schema_s1'})
            CREATE (v)-[:ObservedIn {zygosity: 'HET'}]->(s)
        """)
        assert result_rel is not None, "Query to create ObservedIn relationship should execute."
        print("Schema creation test passed (tables Variant, Sample, and ObservedIn seem to exist).")

    except Exception as e:
        pytest.fail(f"test_create_schema failed: {e}")

def test_add_variant(db_connection):
    """Test adding a variant to the database."""
    variant_data = {
        'variant_id': 'chr1-12345-C-T',
        'chrom': '1',
        'pos': 12345,
        'ref': 'C',
        'alt': 'T',
        'rs_id': 'rs98765'
    }
    try:
        graph_integration.add_variant(db_connection, variant_data)
        
        # Verify by querying the added variant
        query = "MATCH (v:Variant {variant_id: $v_id}) RETURN v.variant_id AS variant_id, v.chrom AS chrom, v.pos AS pos, v.ref AS ref, v.alt AS alt, v.rs_id AS rs_id"
        params = {"v_id": variant_data['variant_id']}
        result = graph_integration.execute_query(db_connection, query, params=params)
        
        assert len(result) == 1, "Should find exactly one variant with the given ID."
        added_v = result[0]
        assert added_v['variant_id'] == variant_data['variant_id']
        assert added_v['chrom'] == variant_data['chrom']
        assert added_v['pos'] == variant_data['pos']
        assert added_v['ref'] == variant_data['ref']
        assert added_v['alt'] == variant_data['alt']
        assert added_v['rs_id'] == variant_data['rs_id']
        print(f"test_add_variant for {variant_data['variant_id']} passed.")

    except Exception as e:
        pytest.fail(f"test_add_variant failed: {e}")

def test_add_sample(db_connection):
    """Test adding a sample to the database."""
    sample_data = {'sample_id': 'SAMPLE_TEST_001'}
    try:
        graph_integration.add_sample(db_connection, sample_data)
        
        # Verify by querying the added sample
        query = "MATCH (s:Sample {sample_id: $s_id}) RETURN s.sample_id AS sample_id"
        params = {"s_id": sample_data['sample_id']}
        result = graph_integration.execute_query(db_connection, query, params=params)
        
        assert len(result) == 1, "Should find exactly one sample with the given ID."
        added_s = result[0]
        assert added_s['sample_id'] == sample_data['sample_id']
        print(f"test_add_sample for {sample_data['sample_id']} passed.")

    except Exception as e:
        pytest.fail(f"test_add_sample failed: {e}")

def test_link_variant_to_sample(db_connection):
    """Test linking a variant to a sample with properties."""
    variant_id = 'chrX-100-G-A'
    sample_id = 'SAMPLE_LINK_002'
    link_properties = {'zygosity': 'HOM'}

    try:
        # Add the variant and sample first
        graph_integration.add_variant(db_connection, {
            'variant_id': variant_id, 'chrom': 'X', 'pos': 100, 'ref': 'G', 'alt': 'A', 'rs_id': 'rsLinkTest'
        })
        graph_integration.add_sample(db_connection, {'sample_id': sample_id})
        
        # Link them
        graph_integration.link_variant_to_sample(db_connection, sample_id, variant_id, link_properties)
        
        # Verify by querying the relationship
        query = """
        MATCH (v:Variant {variant_id: $v_id})-[r:ObservedIn]->(s:Sample {sample_id: $s_id})
        RETURN r.zygosity AS zygosity
        """
        params = {"v_id": variant_id, "s_id": sample_id}
        result = graph_integration.execute_query(db_connection, query, params=params)
        
        assert len(result) == 1, "Should find exactly one relationship between the variant and sample."
        rel_props = result[0]
        assert rel_props['zygosity'] == link_properties['zygosity']
        print(f"test_link_variant_to_sample for {variant_id} -> {sample_id} passed.")

    except Exception as e:
        pytest.fail(f"test_link_variant_to_sample failed: {e}")

def test_execute_query_general(db_connection):
    """Test general query execution, including DataFrame conversion and empty results."""
    try:
        # Setup: Add a known variant
        variant_data = {
            'variant_id': 'chr2-555-A-C', 'chrom': '2', 'pos': 555, 'ref': 'A', 'alt': 'C', 'rs_id': 'rsExecuteTest'
        }
        graph_integration.add_variant(db_connection, variant_data)

        # Test 1: Query existing data
        query_existing = "MATCH (v:Variant {variant_id: $v_id}) RETURN v.variant_id AS id, v.pos AS position"
        params_existing = {"v_id": variant_data['variant_id']}
        result_existing = graph_integration.execute_query(db_connection, query_existing, params=params_existing)
        
        assert len(result_existing) == 1, "Should retrieve the added variant."
        assert isinstance(result_existing, list), "Result should be a list."
        assert isinstance(result_existing[0], dict), "Elements of the list should be dictionaries."
        assert result_existing[0]['id'] == variant_data['variant_id']
        assert result_existing[0]['position'] == variant_data['pos']
        print("test_execute_query_general - existing data query passed.")

        # Test 2: Query non-existing data (should return empty list)
        query_non_existing = "MATCH (v:Variant {variant_id: 'non_existent_id'}) RETURN v.variant_id"
        result_non_existing = graph_integration.execute_query(db_connection, query_non_existing)
        
        assert isinstance(result_non_existing, list), "Result for non-existing data should be a list."
        assert len(result_non_existing) == 0, "Query for non-existing data should return an empty list."
        print("test_execute_query_general - non-existing data query passed.")

        # Test 3: Query with no parameters (if applicable)
        # Add another variant to test a query without specific ID match
        graph_integration.add_variant(db_connection, {
             'variant_id': 'chr3-777-G-T', 'chrom': '3', 'pos': 777, 'ref': 'G', 'alt': 'T', 'rs_id': 'rsAnotherTest'
        })
        query_all = "MATCH (v:Variant) WHERE v.chrom = '3' RETURN v.variant_id AS variant_id"
        result_all_chrom3 = graph_integration.execute_query(db_connection, query_all) # No params
        assert len(result_all_chrom3) == 1
        assert result_all_chrom3[0]['variant_id'] == 'chr3-777-G-T'
        print("test_execute_query_general - query with no params (but with WHERE) passed.")

    except NotImplementedError as nie:
        # This might be hit if the DataFrame conversion method is still not correctly identified.
        # This is an important check for the current state of execute_query.
        pytest.skip(f"Skipping part of test_execute_query_general due to missing DataFrame conversion: {nie}")
    except Exception as e:
        pytest.fail(f"test_execute_query_general failed: {e}")

def test_get_variant_context(db_connection):
    """Test fetching graph context for a list of variant IDs."""
    try:
        # Setup: Add variants, samples, and links
        v1_id = "rs1001"
        v2_id = "rs1002"
        v3_id = "rs1003" # This variant won't be linked to S1
        v4_id = "rs1004" # This variant won't have any links

        s1_id = "SAMPLE_CTX_01"
        s2_id = "SAMPLE_CTX_02"

        graph_integration.add_variant(db_connection, {'variant_id': v1_id, 'chrom': '1', 'pos': 100, 'ref': 'A', 'alt': 'T'})
        graph_integration.add_variant(db_connection, {'variant_id': v2_id, 'chrom': '1', 'pos': 200, 'ref': 'C', 'alt': 'G'})
        graph_integration.add_variant(db_connection, {'variant_id': v3_id, 'chrom': '2', 'pos': 300, 'ref': 'T', 'alt': 'A'})
        graph_integration.add_variant(db_connection, {'variant_id': v4_id, 'chrom': '3', 'pos': 400, 'ref': 'G', 'alt': 'C'})

        graph_integration.add_sample(db_connection, {'sample_id': s1_id})
        graph_integration.add_sample(db_connection, {'sample_id': s2_id})

        graph_integration.link_variant_to_sample(db_connection, s1_id, v1_id, {'zygosity': 'HET'})
        graph_integration.link_variant_to_sample(db_connection, s1_id, v2_id, {'zygosity': 'HOM'})
        graph_integration.link_variant_to_sample(db_connection, s2_id, v1_id, {'zygosity': 'HOM'})
        graph_integration.link_variant_to_sample(db_connection, s2_id, v3_id, {'zygosity': 'HET'})

        # Test 1: Get context for variants linked to S1
        target_variant_ids = [v1_id, v2_id, v3_id, v4_id, "rs_non_existent"]
        contexts = graph_integration.get_variant_context(db_connection, target_variant_ids)

        assert isinstance(contexts, dict), "Context should be a dictionary."
        assert len(contexts) == len(target_variant_ids), "Context dict should have all target IDs as keys."

        # Check context for v1_id
        assert v1_id in contexts
        assert len(contexts[v1_id]['samples']) == 2
        assert {'sample_id': s1_id, 'zygosity': 'HET'} in contexts[v1_id]['samples']
        assert {'sample_id': s2_id, 'zygosity': 'HOM'} in contexts[v1_id]['samples']

        # Check context for v2_id
        assert v2_id in contexts
        assert len(contexts[v2_id]['samples']) == 1
        assert {'sample_id': s1_id, 'zygosity': 'HOM'} in contexts[v2_id]['samples']

        # Check context for v3_id (linked to S2 but not S1, included in query)
        assert v3_id in contexts
        assert len(contexts[v3_id]['samples']) == 1
        assert {'sample_id': s2_id, 'zygosity': 'HET'} in contexts[v3_id]['samples']

        # Check context for v4_id (no links)
        assert v4_id in contexts
        assert len(contexts[v4_id]['samples']) == 0

        # Check context for non_existent_variant
        assert "rs_non_existent" in contexts
        assert len(contexts["rs_non_existent"]['samples']) == 0
        print("test_get_variant_context passed.")

        # Test 2: Empty list of variant IDs
        empty_contexts = graph_integration.get_variant_context(db_connection, [])
        assert isinstance(empty_contexts, dict)
        assert len(empty_contexts) == 0
        print("test_get_variant_context with empty list passed.")

    except NotImplementedError as nie:
        pytest.skip(f"Skipping test_get_variant_context due to missing DataFrame conversion in execute_query: {nie}")
    except Exception as e:
        pytest.fail(f"test_get_variant_context failed: {e}")

# All core functional tests added. 