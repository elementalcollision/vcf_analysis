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
from vcf_agent.bcftools_integration import bcftools_view # Example bcftools call wrapper
from vcf_agent.agent import run_llm_analysis_task # For AI metrics test
from vcf_agent.config import SessionConfig # For AI metrics test
# from vcf_agent.agent import run_llm_analysis_task # For AI metrics test

# It's hard to reliably stop the prometheus_client.start_http_server as it runs in a separate thread
# and doesn't return a server object. We'll rely on daemon=True and pick a unique port.
TEST_METRICS_PORT = 8999 # Choose a port unlikely to be in use

@pytest.fixture(scope="module")
def metrics_server():
    """Starts the VCF Agent's metrics HTTP server on a test port."""
    # Check if port is already in use (basic check)
    try:
        with socket.create_connection(("localhost", TEST_METRICS_PORT), timeout=0.1):
            pytest.skip(f"Port {TEST_METRICS_PORT} is already in use. Skipping metrics server tests.")
    except (socket.timeout, ConnectionRefusedError):
        pass # Port is likely free

    # Ensure the http_registry is clean for module-scoped tests if needed, 
    # though typically metrics accumulate globally. For true isolation, a custom registry 
    # would be passed to start_http_server and used by instrumented code, 
    # but our current setup uses the global metrics.http_registry.
    
    # The start_metrics_http_server in metrics.py uses its own thread management.
    # We call it here to ensure it's running for the tests.
    metrics.start_metrics_http_server(port=TEST_METRICS_PORT)
    
    # Give the server a moment to start
    time.sleep(0.5) 
    
    yield f"http://localhost:{TEST_METRICS_PORT}"
    
    # No explicit teardown for start_http_server as it's a daemon thread.
    # If we had a server object, we'd call server.shutdown() here.

def test_metrics_endpoint_exposes_data(metrics_server_url):
    """Test that the /metrics endpoint is available and serves data."""
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
    """Helper to parse metrics text and find a specific metric value."""
    for family in text_string_to_metric_families(metrics_text):
        if family.name == metric_name:
            for sample in family.samples:
                if sample.labels == labels:
                    return sample.value
    return None

def test_tool_invocation_updates_metrics(metrics_server_url):
    """Test that invoking an agent tool updates relevant metrics."""
    tool_name_to_test = "validate_vcf"
    # Use a non-existent file to ensure a predictable error status for the metric
    dummy_filepath = "/tmp/non_existent_dummy_file_for_test.vcf"
    
    # 1. Get metrics BEFORE tool invocation
    response_before = requests.get(f"{metrics_server_url}/metrics", timeout=1)
    response_before.raise_for_status()
    metrics_before_text = response_before.text

    # Expected labels for an error case in validate_vcf
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

    # 4. Assert metric increments
    requests_total_after = get_metric_value(metrics_after_text, "vcf_agent_tool_requests_total", expected_labels) or 0.0
    duration_count_after = get_metric_value(metrics_after_text, "vcf_agent_tool_duration_seconds_count", expected_labels) or 0.0
    
    assert requests_total_after == requests_total_before + 1, \
        f"Tool requests total did not increment correctly for {tool_name_to_test}. Before: {requests_total_before}, After: {requests_total_after}"
    assert duration_count_after == duration_count_before + 1, \
        f"Tool duration count did not increment correctly for {tool_name_to_test}. Before: {duration_count_before}, After: {duration_count_after}"

def test_bcftools_invocation_updates_metrics(metrics_server_url):
    """Test that invoking a bcftools command updates relevant metrics."""
    bcftools_subcommand_to_test = "view" # bcftools_view uses "view"
    # Using --version should be a quick, successful command
    bcftools_args = ["--version"]
    
    # 1. Get metrics BEFORE bcftools invocation
    response_before = requests.get(f"{metrics_server_url}/metrics", timeout=1)
    response_before.raise_for_status()
    metrics_before_text = response_before.text

    expected_labels = {"bcftools_subcommand": bcftools_subcommand_to_test, "status": "success"}

    commands_total_before = get_metric_value(metrics_before_text, "vcf_agent_bcftools_commands_total", expected_labels) or 0.0
    duration_count_before = get_metric_value(metrics_before_text, "vcf_agent_bcftools_duration_seconds_count", expected_labels) or 0.0
    errors_total_value_before = get_metric_value(metrics_before_text, "vcf_agent_bcftools_errors_total", {"bcftools_subcommand": bcftools_subcommand_to_test}) or 0.0


    # 2. Invoke the bcftools command wrapper
    # bcftools_view calls run_bcftools_command which is instrumented
    try:
        return_code, stdout, stderr = bcftools_view(bcftools_args)
        assert return_code == 0, f"bcftools view --version failed. stderr: {stderr}"
    except FileNotFoundError:
        pytest.skip("bcftools not found in PATH, skipping bcftools metrics test.")
    except Exception as e:
        pytest.fail(f"bcftools_view raised an unexpected exception: {e}")

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

def test_ai_interaction_updates_metrics(metrics_server_url):
    """Test that an AI interaction (even if it fails) updates relevant AI metrics."""
    # This test assumes that an actual LLM call might fail if not configured,
    # which is fine as we're testing the metrics for error paths too.
    task_name = "test_ai_metric_task"
    model_provider_to_test = "ollama" # Could be any, ollama is default if no infra

    # Expected labels for an error case (most likely if Ollama isn't running)
    # The actual error type might vary.
    # We are primarily interested in the request count and duration observation.
    # Status will likely be "error".
    expected_labels = {"model_provider": model_provider_to_test, "endpoint_task": task_name, "status": "error"}
    
    # 1. Get metrics BEFORE AI call
    response_before = requests.get(f"{metrics_server_url}/metrics", timeout=1)
    response_before.raise_for_status()
    metrics_before_text = response_before.text

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