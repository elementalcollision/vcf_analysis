import pytest
import requests
import time
import threading
from http.server import HTTPServer # To check if port is free, though not foolproof
import socket
from prometheus_client.parser import text_string_to_metric_families
import os # For creating a dummy file if needed
from typing import Optional

from vcf_agent import metrics
from vcf_agent.agent import validate_vcf # Example tool
from vcf_agent.bcftools_integration import run_bcftools_command
from vcf_agent.agent import run_llm_analysis_task # For AI metrics test
from vcf_agent.config import SessionConfig # For AI metrics test
# from vcf_agent.agent import run_llm_analysis_task # For AI metrics test
import prometheus_client
import pytest

@pytest.fixture(autouse=True)
def reset_prometheus_registry():
    # Unregister all custom collectors (except the default process/python_gc ones)
    collectors = list(prometheus_client.REGISTRY._collector_to_names.keys())
    for collector in collectors:
        try:
            prometheus_client.REGISTRY.unregister(collector)
        except KeyError:
            pass
# It's hard to reliably stop the prometheus_client.start_http_server as it runs in a separate thread
# and doesn't return a server object. We'll rely on daemon=True and pick a unique port.
TEST_METRICS_PORT = 8999 # Choose a port unlikely to be in use

@pytest.fixture(scope="module")
def metrics_server():
    """Starts the VCF Agent's metrics HTTP server on a test port."""
    # Clear the global registry before starting the server for this module
    # This helps avoid interference from other tests if they use the same global registry.
    # Get a list of collectors to unregister, as unregistering modifies the set
    # collectors_to_remove = list(metrics.REGISTRY._collector_to_names.keys())
    # for collector in collectors_to_remove:
    #     try:
    #         metrics.REGISTRY.unregister(collector)
    #     except KeyError: # pragma: no cover
    #         pass # May have already been unregistered by a dependent collector

    # Attempt to clear default collectors that might be auto-registered by prometheus_client
    # This is a bit of a blunt instrument.
    from prometheus_client import REGISTRY, PROCESS_COLLECTOR, PLATFORM_COLLECTOR
    # Unregister default collectors if they are registered
    # Need to do this carefully as unregistering non-existent collectors raises KeyError
    default_collectors = []
    if PROCESS_COLLECTOR in REGISTRY._collector_to_names: # type: ignore
        default_collectors.append(PROCESS_COLLECTOR)
    if PLATFORM_COLLECTOR in REGISTRY._collector_to_names: # type: ignore
        default_collectors.append(PLATFORM_COLLECTOR)
    
    # Unregister any application-specific collectors that might have been added to the global REGISTRY
    # This is hard to do generically without knowing their names.
    # For now, we'll rely on our metrics being prefixed with 'vcf_agent_'
    # and hope that `metrics.py` re-registers them.
    # A better approach is a custom registry for tests.

    # The above attempts to clear REGISTRY are complex due to its global nature.
    # A simpler, though less complete, approach for now is to ensure our server starts fresh.

    # Check if port is already in use (basic check)
    try:
        with socket.create_connection(("localhost", TEST_METRICS_PORT), timeout=0.1): # type: ignore
            pytest.skip(f"Port {TEST_METRICS_PORT} is already in use. Skipping metrics server tests.")
    except (socket.timeout, ConnectionRefusedError):
        pass # Port is likely free

    metrics.start_metrics_http_server(port=TEST_METRICS_PORT)
    
    # Give the server a moment to start, with retries
    server_url = f"http://localhost:{TEST_METRICS_PORT}"
    for _ in range(10): # Try for up to ~2.5 seconds
        try:
            response = requests.get(f"{server_url}/metrics", timeout=0.2)
            if response.status_code == 200:
                break
        except requests.exceptions.ConnectionError:
            time.sleep(0.25)
    else: # pragma: no cover
        pytest.fail(f"Metrics server did not become available at {server_url} after multiple retries.")
    
    yield server_url
    
    # No explicit teardown for start_http_server as it's a daemon thread.
    # If we had a server object, we'd call server.shutdown() here.

def test_metrics_endpoint_exposes_data(metrics_server):
    """Test that the /metrics endpoint is available and serves data."""
    metrics_server_url = metrics_server # Use the yielded value
    try:
        response = requests.get(f"{metrics_server_url}/metrics", timeout=1)
        response.raise_for_status() # Check for HTTP errors
    except requests.exceptions.ConnectionError:
        pytest.fail(f"Could not connect to metrics server at {metrics_server_url}. Ensure it started correctly.")

    assert response.status_code == 200
    assert "text/plain" in response.headers["Content-Type"]
    content = response.text
    assert "# HELP" in content # Basic check for Prometheus format
    assert "# TYPE" in content
    
    # Check for some known metric names (actual values not checked here)
    assert "vcf_agent_tool_requests_total" in content
    assert "vcf_agent_bcftools_commands_total" in content
    assert "vcf_agent_ai_requests_total" in content

from prometheus_client.parser import text_string_to_metric_families
import os # For creating a dummy file if needed

def get_metric_value(metrics_text: str, metric_name: str, labels: dict) -> Optional[float]:
    """Helper to parse metrics text and find a specific metric value. Adds debug output."""
    print(f"[DEBUG] Looking for metric: {repr(metric_name)} with labels: {repr(labels)}")
    
    for family in text_string_to_metric_families(metrics_text):
        print(f"[DEBUG] Encountered Family: name={repr(family.name)}, type={family.type}")
        
        # The metric_name from the test might be the full name (e.g., my_metric_total)
        # The family.name from the parser might be the base name (e.g., my_metric)
        # So, check if the queried metric_name starts with the family.name
        # This is relevant for counters (_total) and summaries/histograms (_count, _sum, _bucket)
        if not metric_name.startswith(family.name):
            # Also, family.name might be longer if metric_name is just a base for a suffixed family.
            # Example: metric_name = "request_duration", family.name = "request_duration_seconds"
            # This case is less common if tests use full names, but good to be aware of.
            # For now, primary focus is metric_name (e.g., "foo_total") starts with family.name (e.g., "foo")
            if not family.name.startswith(metric_name): # Handle cases where family.name might be more specific (e.g. _seconds)
                 print(f"[DEBUG] Skipping family {repr(family.name)} as it does not form base for {repr(metric_name)}")
                 continue
            # If family.name starts with metric_name, it implies metric_name is a base and family.name is suffixed
            # This means we should check samples in this family if their full name matches the *family.name*
            # This is complex. For now, the main hypothesis is metric_name is suffixed, family.name is base.

        print(f"[DEBUG] Potential Family Match: Querying {repr(metric_name)}, Family base is {repr(family.name)}. Checking samples.")

        for sample in family.samples:
            # Sample.name is the fully qualified name including suffixes like _total, _count, _sum, _bucket
            print(f"[DEBUG]  Sample: name={repr(sample.name)}, labels={repr(sample.labels)}, value={sample.value}")
            
            # We must match the FULL metric_name (e.g., 'vcf_agent_tool_requests_total')
            # with the sample's full name.
            if sample.name == metric_name:
                print(f"[DEBUG]    Exact Sample Name Match: {repr(sample.name)} == {repr(metric_name)}")
                sample_labels_as_dict = dict(sample.labels)
                if sample_labels_as_dict == labels:
                    print(f"[DEBUG]      Found LABEL match for {metric_name} with labels {labels}: {sample.value}")
                    return sample.value
                else:
                    print(f"[DEBUG]      Label MISMATCH for {metric_name}. Expected: {labels}, Got: {sample_labels_as_dict}")
            else:
                print(f"[DEBUG]    Sample Name Mismatch: {repr(sample.name)} != {repr(metric_name)}")
        
        # If we checked all samples in a potentially relevant family and didn't return, 
        # it means no sample within that family was a full match for both name and labels.
        # Do not return None here, continue to the next family.

    print(f"[DEBUG] Exhausted all families. No match found overall for {metric_name} with labels {labels}")
    return None

def test_tool_invocation_updates_metrics(metrics_server):
    """Test that invoking an agent tool updates relevant metrics."""
    metrics_server_url = metrics_server # Use the yielded value
    tool_name_to_test = "validate_vcf"
    # Use a non-existent file to ensure a predictable error status for the metric
    dummy_filepath = "/tmp/non_existent_dummy_file_for_test.vcf"
    
    # 1. Get metrics BEFORE tool invocation
    response_before = requests.get(f"{metrics_server_url}/metrics", timeout=1)
    response_before.raise_for_status()
    metrics_before_text = response_before.text
    # print(f"\n--- Metrics BEFORE tool ({tool_name_to_test}) ---\n{metrics_before_text}\n---------------------") # DEBUG

    # Expected labels for an error case in validate_vcf (returns error message, not exception)
    expected_labels = {"tool_name": tool_name_to_test, "status": "error"}

    # Get initial counter values (or 0 if not present yet)
    requests_total_before = get_metric_value(metrics_before_text, "vcf_agent_tool_requests_total", expected_labels) or 0.0
    duration_count_before = get_metric_value(metrics_before_text, "vcf_agent_tool_duration_seconds_count", expected_labels) or 0.0

    # 2. Invoke the tool
    try:
        # validate_vcf is already decorated and will use global metrics.http_registry
        validate_vcf(dummy_filepath) 
    except Exception:
        # The tool itself might raise an exception if the underlying validate_vcf_file does, 
        # or it might return an error string. The metric status should be 'error'.
        pass # We expect errors, metrics should capture this.

    # Give a very brief moment for metrics to be processed if there was any async aspect (unlikely for direct calls)
    time.sleep(0.1)

    # 3. Get metrics AFTER tool invocation
    response_after = requests.get(f"{metrics_server_url}/metrics", timeout=1)
    response_after.raise_for_status()
    metrics_after_text = response_after.text
    print(f"\n--- Metrics AFTER tool ({tool_name_to_test}) ---\n{metrics_after_text}\n---------------------") # DEBUG

    # 4. Assert metric increments
    requests_total_after = get_metric_value(metrics_after_text, "vcf_agent_tool_requests_total", expected_labels) or 0.0
    duration_count_after = get_metric_value(metrics_after_text, "vcf_agent_tool_duration_seconds_count", expected_labels) or 0.0
    
    assert requests_total_after == requests_total_before + 1, \
        f"Tool requests total did not increment correctly for {tool_name_to_test}. Before: {requests_total_before}, After: {requests_total_after}"
    assert duration_count_after == duration_count_before + 1, \
        f"Tool duration count did not increment correctly for {tool_name_to_test}. Before: {duration_count_before}, After: {duration_count_after}"

def test_bcftools_invocation_updates_metrics(metrics_server):
    """Test that invoking a bcftools command updates relevant metrics."""
    metrics_server_url = metrics_server # Use the yielded value
    bcftools_subcommand_to_test = "--version" # run_bcftools_command uses "--version"
    # 1. Get metrics BEFORE bcftools invocation
    response_before = requests.get(f"{metrics_server_url}/metrics", timeout=1)
    response_before.raise_for_status()
    metrics_before_text = response_before.text

    expected_labels = {"bcftools_subcommand": bcftools_subcommand_to_test, "status": "success"}

    commands_total_before = get_metric_value(metrics_before_text, "vcf_agent_bcftools_commands_total", expected_labels) or 0.0
    duration_count_before = get_metric_value(metrics_before_text, "vcf_agent_bcftools_duration_seconds_count", expected_labels) or 0.0
    errors_total_value_before = get_metric_value(metrics_before_text, "vcf_agent_bcftools_errors_total", {"bcftools_subcommand": bcftools_subcommand_to_test}) or 0.0

    # 2. Invoke the bcftools command wrapper
    try:
        return_code, stdout, stderr = run_bcftools_command(["--version"])
        assert return_code == 0, f"bcftools --version failed. stdout: {stdout}, stderr: {stderr}"
        assert "bcftools" in stdout.lower(), f"Expected 'bcftools' in version output. Got: {stdout}"
    except FileNotFoundError: # pragma: no cover
        pytest.skip("bcftools not found in PATH, skipping bcftools metrics test.")
    except Exception as e: # pragma: no cover
        pytest.fail(f"run_bcftools_command with --version raised an unexpected exception: {e}")

    time.sleep(0.1) # Brief pause for metrics processing if any async behavior

    # 3. Get metrics AFTER bcftools invocation
    response_after = requests.get(f"{metrics_server_url}/metrics", timeout=1)
    response_after.raise_for_status()
    metrics_after_text = response_after.text

    # 4. Assert metric increments
    commands_total_after = get_metric_value(metrics_after_text, "vcf_agent_bcftools_commands_total", expected_labels) or 0.0
    duration_count_after = get_metric_value(metrics_after_text, "vcf_agent_bcftools_duration_seconds_count", expected_labels) or 0.0
    errors_total_value_after = get_metric_value(metrics_after_text, "vcf_agent_bcftools_errors_total", {"bcftools_subcommand": bcftools_subcommand_to_test}) or 0.0

    assert commands_total_after == commands_total_before + 1, \
        f"BCFTools commands total did not increment. Before: {commands_total_before}, After: {commands_total_after}"
    assert duration_count_after == duration_count_before + 1, \
        f"BCFTools duration count did not increment. Before: {duration_count_before}, After: {duration_count_after}"
    assert errors_total_value_after == errors_total_value_before, \
        f"BCFTools errors total changed unexpectedly for a success case. Before: {errors_total_value_before}, After: {errors_total_value_after}"

def test_ai_interaction_updates_metrics(metrics_server):
    """Test that an AI interaction (even if it fails) updates relevant AI metrics."""
    metrics_server_url = metrics_server # Use the yielded value
    # This test assumes that an actual LLM call might fail if not configured,
    # which is fine as we're testing the metrics for error paths too.
    task_name = "test_ai_metric_task"
    model_provider_to_test = "ollama" # Could be any, ollama is default if no infra

    # Expected labels for an error case (most likely if Ollama isn't running)
    # The actual error type might vary.
    # We are primarily interested in the request count and duration observation.
    # Status will likely be "error".
    expected_labels = {"model_provider": model_provider_to_test, "endpoint_task": task_name, "status": "success"}
    
    # 1. Get metrics BEFORE AI call
    response_before = requests.get(f"{metrics_server_url}/metrics", timeout=1)
    response_before.raise_for_status()
    metrics_before_text = response_before.text
    print(f"\n--- Metrics BEFORE AI ({task_name}) ---\n{metrics_before_text}\n---------------------") # DEBUG

    requests_total_before = get_metric_value(metrics_before_text, "vcf_agent_ai_requests_total", expected_labels) or 0.0
    duration_count_before = get_metric_value(metrics_before_text, "vcf_agent_ai_response_seconds_count", expected_labels) or 0.0
    concurrent_gauge_value_before = get_metric_value(metrics_before_text, "vcf_agent_ai_concurrent_requests", {"model_provider": model_provider_to_test, "endpoint_task": task_name}) or 0.0

    # 2. Invoke the AI task runner
    session_cfg = SessionConfig(raw_mode=True, model_provider=model_provider_to_test)
    try:
        run_llm_analysis_task(
            task=task_name, 
            file_paths=[], 
            session_config=session_cfg,
            model_provider=model_provider_to_test
        )
    except Exception as e:
        print(f"AI task call failed as expected for metrics test: {e}")
        pass 

    time.sleep(0.2)

    # 3. Get metrics AFTER AI call
    response_after = requests.get(f"{metrics_server_url}/metrics", timeout=1)
    response_after.raise_for_status()
    metrics_after_text = response_after.text
    print(f"\n--- Metrics AFTER AI ({task_name}) ---\n{metrics_after_text}\n---------------------") # DEBUG

    # 4. Assert metric increments
    requests_total_after = get_metric_value(metrics_after_text, "vcf_agent_ai_requests_total", expected_labels) or 0.0
    duration_count_after = get_metric_value(metrics_after_text, "vcf_agent_ai_response_seconds_count", expected_labels) or 0.0
    concurrent_gauge_value_after = get_metric_value(metrics_after_text, "vcf_agent_ai_concurrent_requests", {"model_provider": model_provider_to_test, "endpoint_task": task_name}) or 0.0

    assert requests_total_after == requests_total_before + 1, \
        f"AI requests total did not increment. Before: {requests_total_before}, After: {requests_total_after}"
    assert duration_count_after == duration_count_before + 1, \
        f"AI duration count did not increment. Before: {duration_count_before}, After: {duration_count_after}"
    assert concurrent_gauge_value_before == 0.0, "Concurrent AI requests should be 0 before test call"
    assert concurrent_gauge_value_after == 0.0, "Concurrent AI requests should be 0 after test call completion/error"

# Test for AI interaction metrics would be next.
# More tests will be added here for specific metric updates. 