import pytest
import shutil
import tempfile
from pathlib import Path
import subprocess # For running CLI commands
import os
import kuzu
from vcf_agent.lancedb_integration import get_db, get_or_create_table
from vcf_agent.graph_integration import create_schema as create_kuzu_schema # Renamed for clarity
from typing import Iterator # Import Iterator for type hinting generators

# Define paths for test databases
# These could be further configured via environment variables or pytest options if needed
TEST_LANCEDB_PATH_BASE = Path(tempfile.gettempdir()) / "vcf_agent_tests" / "lancedb"
TEST_KUZU_PATH_BASE = Path(tempfile.gettempdir()) / "vcf_agent_tests" / "kuzu"

@pytest.fixture(scope="function") # Changed from session to function scope
def test_dbs(request):
    """
    Pytest fixture to set up and tear down unique test LanceDB and Kuzu databases
    for each test function.
    """
    # Create unique path for each test function to avoid collisions
    # Incorporate module and function name into path if possible, or use a simpler unique ID
    # For simplicity with tmp_path like behavior, pytest creates unique dirs per test func for tmp_path.
    # We will create subdirectories under our base paths, named after the test node ID.
    # Sanitize nodeid: replace characters not suitable for directory names.
    sanitized_nodeid = request.node.nodeid.replace("::", "_").replace("[", "_").replace("]", "").replace("/", "_")
    
    lancedb_path = TEST_LANCEDB_PATH_BASE / sanitized_nodeid
    kuzu_path = TEST_KUZU_PATH_BASE / sanitized_nodeid

    # Cleanup before test: Remove any existing test DB directories for this specific test
    if lancedb_path.exists():
        shutil.rmtree(lancedb_path)
    if kuzu_path.exists():
        shutil.rmtree(kuzu_path)

    # Create fresh directories
    lancedb_path.mkdir(parents=True, exist_ok=True)
    kuzu_path.mkdir(parents=True, exist_ok=True)

    yield {
        "lancedb_path": str(lancedb_path),
        "kuzu_path": str(kuzu_path)
    }

    # Teardown after test
    if lancedb_path.exists():
        shutil.rmtree(lancedb_path, ignore_errors=True) # ignore_errors for robustness
    if kuzu_path.exists():
        shutil.rmtree(kuzu_path, ignore_errors=True)   # ignore_errors for robustness

# Fixture for providing a sample VCF file path
@pytest.fixture(scope="session")
def sample_vcf_file_small():
    """
    Provides the path to a small, valid sample VCF file for testing.
    Assumes 'sample_test_data/small_valid.vcf' exists.
    Adjust path as per your project structure.
    """
    # TODO: Ensure this path is correct for your project or create a small VCF here
    sample_file = Path("sample_test_data/small_valid.vcf")
    if not sample_file.exists():
        # Create a minimal VCF for testing if it doesn't exist
        sample_file.parent.mkdir(parents=True, exist_ok=True)
        with open(sample_file, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
            f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
            f.write("chr1\t100\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\n")
            f.write("chr1\t200\trs2\tC\tT\t60\tPASS\tDP=20\tGT\t1/1\n")
    return str(sample_file)

@pytest.fixture(scope="function") # Function scope, as content might vary or be modified by tests if not careful
def invalid_vcf_file(tmp_path):
    """
    Creates a temporary, malformed VCF file for testing error handling.
    Example: Missing VCF header.
    """
    invalid_file = tmp_path / "invalid.vcf"
    with open(invalid_file, "w") as f:
        # Malformed content: e.g., data line without proper header
        f.write("#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tsample1\\n") # Header is okay
        f.write("chr1\\t100\\trs1\\tA\\tG\\t50\\tPASS\\tDP=10\\tGT\\t0/1\\tEXTRA_COLUMN\\n") # Extra column
    return str(invalid_file)

@pytest.fixture(scope="session")
def vcf_agent_cli_runner():
    """
    Provides a helper function to run VCF agent CLI commands.
    Sets PYTHONPATH and CWD to ensure the module is found.
    """
    # Determine project root: assuming conftest.py is in tests/integration/e2e/
    # So, project_root is three levels up from this file's directory.
    project_root = Path(__file__).resolve().parent.parent.parent.parent
    src_dir = project_root / "src"

    def runner(command_args: list):
        base_command = ["python", "-m", "vcf_agent.cli"]
        full_command = base_command + command_args
        
        # Get current environment and update PYTHONPATH
        env = os.environ.copy()
        current_python_path = env.get("PYTHONPATH", "")
        # Prepend src_dir to PYTHONPATH
        env["PYTHONPATH"] = f"{str(src_dir)}{os.pathsep}{current_python_path}" if current_python_path else str(src_dir)

        try:
            result = subprocess.run(
                full_command,
                capture_output=True,
                text=True,
                check=False,
                timeout=60,
                cwd=str(project_root), # Run from project root
                env=env # Pass the modified environment
            )
            return result
        except subprocess.TimeoutExpired:
            pytest.fail(f"CLI command timed out: {' '.join(full_command)}")
        except FileNotFoundError: # Should be less likely now with CWD and PYTHONPATH
            pytest.fail(f"CLI command not found. Ensure vcf_agent.cli is executable and python is in PATH. Command: {' '.join(full_command)}")
        except Exception as e: # Catch other potential subprocess errors
            pytest.fail(f"CLI command failed with an unexpected error: {e}. Command: {' '.join(full_command)}")
            
    return runner 

@pytest.fixture(scope="function")
def unavailable_lancedb_dir(tmp_path: Path) -> Iterator[Path]:
    original_path = tmp_path / "lancedb_test_data"
    unavailable_path = tmp_path / "lancedb_test_data_unavailable"

    # Ensure a clean state: remove both if they exist from previous runs
    if original_path.exists():
        shutil.rmtree(original_path)
    if unavailable_path.exists():
        shutil.rmtree(unavailable_path)

    # Create the original path, then rename it to make it "unavailable"
    original_path.mkdir(parents=True, exist_ok=True) 
    # You might want to initialize a dummy table here if tests expect an existing but inaccessible DB
    # For example:
    db = get_db(str(original_path))
    get_or_create_table(db, "variants") 
    
    os.rename(original_path, unavailable_path)
    
    yield original_path # This is the path the test should try to use

    # Cleanup
    # If the test operation recreated the original_path, remove it first
    if original_path.exists() and original_path.is_dir():
        shutil.rmtree(original_path)
    elif original_path.exists(): # If it's a file
        original_path.unlink()

    # Rename the unavailable_path back (if it still exists)
    if unavailable_path.exists():
        os.rename(unavailable_path, original_path)
    
    # Final cleanup of original_path in case renaming failed or it was recreated again
    if original_path.exists(): 
        shutil.rmtree(original_path, ignore_errors=True)


@pytest.fixture(scope="function")
def unavailable_kuzu_dir(tmp_path: Path) -> Iterator[Path]:
    original_path = tmp_path / "kuzu_test_data_for_unavailable_fixture"
    unavailable_path = tmp_path / "kuzu_test_data_unavailable_actual"

    # Ensure a clean state for these specific paths
    if original_path.exists():
        if original_path.is_dir():
            shutil.rmtree(original_path)
        else:
            original_path.unlink()
    if unavailable_path.exists():
        shutil.rmtree(unavailable_path)

    # Create and initialize Kuzu DB at original_path
    original_path.mkdir(parents=True, exist_ok=True)
    try:
        with kuzu.Database(str(original_path)) as db_init:
            with kuzu.Connection(db_init) as conn_init:
                create_kuzu_schema(conn_init)
        print(f"Successfully initialized Kuzu DB for fixture at {original_path}")
    except Exception as e:
        print(f"Warning: Failed to initialize dummy Kuzu DB in unavailable_kuzu_dir fixture at {original_path}: {e}")
        # If init fails, the rest of the fixture logic might not make sense, but proceed to see test behavior

    # Rename the directory where Kuzu DB was initialized
    if original_path.exists():
        os.rename(original_path, unavailable_path)
        print(f"Renamed {original_path} to {unavailable_path}")
    else:
        print(f"Warning: original_path {original_path} did not exist before rename in unavailable_kuzu_dir")

    # Now, original_path (the path to be yielded) should NOT be a Kuzu DB directory.
    # To make it even more certain Kuzu can't use it, create a file at original_path.
    if not original_path.exists(): # original_path should not exist after rename
        try:
            with open(original_path, 'w') as f:
                f.write("This is a file, not a Kuzu DB directory.")
            print(f"Created a dummy file at {original_path} to ensure Kuzu cannot use it as DB dir.")
        except OSError as e:
            print(f"Could not create dummy file at {original_path}: {e}") # e.g. if tmp_path itself became problematic
    else:
        # This case is problematic: original_path still exists after it was supposedly renamed.
        print(f"ERROR in unavailable_kuzu_dir: {original_path} still exists after attempting to rename it.")
        # Try to remove it if it's a dir, or unlink if file, to ensure test targets a non-DB path.
        if original_path.is_dir(): shutil.rmtree(original_path, ignore_errors=True)
        elif original_path.is_file(): original_path.unlink(missing_ok=True)


    yield original_path # This path should now be a file or non-existent, not a Kuzu DB dir

    # Cleanup
    if original_path.exists(): # The dummy file we might have created
        if original_path.is_file():
            original_path.unlink(missing_ok=True)
        elif original_path.is_dir(): # Should not be a dir if logic above worked
             shutil.rmtree(original_path, ignore_errors=True)
            
    if unavailable_path.exists() and unavailable_path.is_dir(): # The renamed actual Kuzu DB dir
        shutil.rmtree(unavailable_path, ignore_errors=True)

@pytest.fixture(scope="function")
def unreadable_file(tmp_path: Path) -> Iterator[Path]:
    """Creates a temporary file and makes it unreadable."""
    file_path = tmp_path / "unreadable.txt"
    file_path.write_text("This content should not be readable by the test.")
    original_mode = file_path.stat().st_mode
    # Remove all read permissions (owner, group, other)
    # Owner write is 0o200, owner read is 0o400. We want to remove 0o400.
    unreadable_mode = original_mode & ~0o444 # Remove r--r--r--
    try:
        file_path.chmod(unreadable_mode)
        yield file_path
    finally:
        # Restore original permissions to allow cleanup
        file_path.chmod(original_mode)
        if file_path.exists():
            file_path.unlink(missing_ok=True) 