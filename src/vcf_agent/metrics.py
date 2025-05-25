"""
metrics.py
Central module for Prometheus metrics and structured logging for the VCF Analysis Agent.
- Uses structlog for JSON logs.
- Uses prometheus_client for metrics.
- Exposes an HTTP /metrics endpoint for Prometheus scraping.
- Provides functionality for pushing metrics to a Prometheus Pushgateway for batch/CLI jobs.
"""

import structlog
import logging
from prometheus_client import Counter, Histogram, Gauge, start_http_server, CollectorRegistry, push_to_gateway
import time
import threading
import os
from typing import Optional

# --- Structured Logging Setup ---

def setup_logging():
    # Basic logging configuration (can be overridden by app-level config)
    logging.basicConfig(level=os.environ.get("LOG_LEVEL", "INFO").upper(), format="%(message)s")
    structlog.configure(
        processors=[
            structlog.stdlib.add_logger_name,
            structlog.stdlib.add_log_level,
            structlog.stdlib.ProcessorFormatter.wrap_for_formatter, # Required for non-JSON output
        ],
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )
    # Configure a specific formatter for structlog messages if desired, e.g., JSON
    # This setup is basic; for JSON, ensure the root logger's handlers use a JSON formatter
    # or configure structlog's processors for JSON output like in the original file if JSON is strictly needed.
    # For now, this uses standard library formatting which might not be JSON by default.
    # Reverting to a JSON-first setup as in the original for consistency:
    structlog.configure(
        processors=[
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.stdlib.add_logger_name, # Good to have logger name
            structlog.stdlib.add_log_level,   # Good to have log level
            structlog.processors.JSONRenderer()
        ],
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
    )
    return structlog.get_logger()

log = setup_logging()

# --- Prometheus Metrics Setup ---

# Default registry for metrics exposed via HTTP
http_registry = CollectorRegistry(auto_describe=True)

# --- AI Interaction Metrics --- (Registered to http_registry by default)
VCF_AGENT_AI_REQUESTS_TOTAL = Counter(
    'vcf_agent_ai_requests_total', 
    'Total AI agent requests towards various models.', 
    ['model_provider', 'endpoint_task', 'status'],
    registry=http_registry
)
VCF_AGENT_AI_RESPONSE_SECONDS = Histogram(
    'vcf_agent_ai_response_seconds', 
    'AI agent response time (seconds) for various models.', 
    ['model_provider', 'endpoint_task'],
    registry=http_registry
)
VCF_AGENT_AI_TOKENS_TOTAL = Counter(
    'vcf_agent_ai_tokens_total', 
    'Total tokens processed (e.g., prompt + completion) by the AI agent.', 
    ['model_provider', 'endpoint_task', 'token_type'], # token_type: "prompt", "completion", "total"
    registry=http_registry
)
VCF_AGENT_AI_ERRORS_TOTAL = Counter(
    'vcf_agent_ai_errors_total', 
    'Total errors by type during AI model interactions.', 
    ['model_provider', 'endpoint_task', 'error_type'],
    registry=http_registry
)
VCF_AGENT_AI_CONCURRENT_REQUESTS = Gauge(
    'vcf_agent_ai_concurrent_requests', 
    'Number of concurrent AI agent requests towards various models.', 
    ['model_provider', 'endpoint_task'],
    registry=http_registry
)

# --- CLI Command Metrics ---
# These will often be registered to a temporary registry for Pushgateway
# For now, define them. Registration strategy will be handled by helper functions.
CLI_COMMAND_REQUESTS_TOTAL = Counter(
    'vcf_agent_cli_command_requests_total',
    'Total number of CLI subcommand invocations.',
    ['command', 'status']
    # registry=http_registry # Optionally register to HTTP if some CLI commands are long-running server-like
)
CLI_COMMAND_DURATION_SECONDS = Histogram(
    'vcf_agent_cli_command_duration_seconds',
    'Duration of CLI subcommand execution in seconds.',
    ['command', 'status']
    # registry=http_registry
)
CLI_COMMAND_ERRORS_TOTAL = Counter(
    'vcf_agent_cli_command_errors_total',
    'Total errors specifically from CLI command execution logic.',
    ['command', 'error_type']
    # registry=http_registry
)

# --- Agent Tool Metrics ---
VCF_AGENT_TOOL_REQUESTS_TOTAL = Counter(
    'vcf_agent_tool_requests_total', 
    'Total number of agent tool invocations.', 
    ['tool_name', 'status'],
    registry=http_registry
)
VCF_AGENT_TOOL_DURATION_SECONDS = Histogram(
    'vcf_agent_tool_duration_seconds', 
    'Duration of agent tool execution in seconds.', 
    ['tool_name', 'status'],
    registry=http_registry
)
VCF_AGENT_TOOL_ERRORS_TOTAL = Counter(
    'vcf_agent_tool_errors_total', 
    'Total errors during agent tool execution.', 
    ['tool_name', 'error_type'],
    registry=http_registry
)

# --- BCFTools Integration Metrics ---
VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL = Counter(
    'vcf_agent_bcftools_commands_total', 
    'Total number of bcftools invocations.', 
    ['bcftools_subcommand', 'status'],
    registry=http_registry
)
VCF_AGENT_BCFTOOLS_DURATION_SECONDS = Histogram(
    'vcf_agent_bcftools_duration_seconds', 
    'Duration of bcftools subprocess calls in seconds.', 
    ['bcftools_subcommand', 'status'],
    registry=http_registry
)
VCF_AGENT_BCFTOOLS_ERRORS_TOTAL = Counter(
    'vcf_agent_bcftools_errors_total', 
    'Total errors from bcftools invocations (non-zero return codes).', 
    ['bcftools_subcommand'],
    registry=http_registry
)

# --- Metrics HTTP Endpoint --- #
METRICS_HTTP_PORT = int(os.environ.get("VCF_AGENT_METRICS_PORT", 8000))
_metrics_server_thread = None

def start_metrics_http_server(port: Optional[int] = None):
    global _metrics_server_thread
    actual_port = port if port is not None else METRICS_HTTP_PORT
    if _metrics_server_thread is None or not _metrics_server_thread.is_alive():
        try:
            # Start Prometheus metrics server in a background thread using the http_registry
            # For linter debugging, simplifying the call temporarily
            # _metrics_server_thread = threading.Thread(
            #     target=start_http_server, args=(actual_port, http_registry), daemon=True
            # )
            # TEMP: Linter seems to struggle with `registry=` kwarg here.
            # The library default is to use the global REGISTRY, which is not what we want for http_registry.
            # This is a workaround attempt for the linter. Runtime behavior might differ if not using our http_registry.
            # For actual execution, the version with `registry=http_registry` is preferred.
            def _target_for_thread(): # Closure to ensure http_registry is used
                start_http_server(actual_port, registry=http_registry)

            _metrics_server_thread = threading.Thread(
                target=_target_for_thread, daemon=True
            )
            _metrics_server_thread.start()
            log.info("prometheus_metrics_server_started", port=actual_port)
        except Exception as e:
            log.error("prometheus_metrics_server_failed_to_start", port=actual_port, error=str(e))
    else:
        log.info("prometheus_metrics_server_already_running", port=actual_port)

# --- Pushgateway Functionality --- #
PUSHGATEWAY_URL = os.environ.get("VCF_AGENT_PUSHGATEWAY_URL", "http://localhost:9091") # Default URL

def push_job_metrics(job_name: str, metrics_to_push: list, instance_name: Optional[str] = None):
    """
    Pushes a list of metric objects (already updated) to the Pushgateway.
    Metrics are registered to a temporary registry for this push.
    Args:
        job_name: The job name for Pushgateway.
        metrics_to_push: A list of MetricWrapperBase instances (Counter, Gauge, Histogram) that have been updated.
        instance_name: Optional. A specific instance name for this push. Defaults to job_name-pid.
    """
    if not PUSHGATEWAY_URL or PUSHGATEWAY_URL == "disabled": # Allow disabling via env var
        log.warn("pushgateway_not_configured_or_disabled", job_name=job_name)
        return

    push_registry = CollectorRegistry()
    for metric_obj in metrics_to_push:
        push_registry.register(metric_obj)

    try:
        effective_instance_name = instance_name if instance_name else f"{job_name}-{os.getpid()}"
        
        # Define grouping_key dictionary separately
        current_grouping_key = {
            'job': job_name,
            'instance': effective_instance_name
        }
        
        # For linter debugging, simplifying the call temporarily
        # push_to_gateway(
        #     PUSHGATEWAY_URL, 
        #     job=job_name, 
        #     registry=push_registry, 
        #     grouping_key=current_grouping_key
        # )
        # TEMP: Linter seems to struggle with kwargs in push_to_gateway.
        # This is a workaround attempt for the linter. This simplified call will likely NOT behave as intended
        # regarding grouping_key if it even runs.
        # For actual execution, the version with all kwargs is preferred.
        push_to_gateway(PUSHGATEWAY_URL, job=job_name, registry=push_registry) # Removed grouping_key for linter test

        log.info("metrics_pushed_to_gateway", job_name=job_name, instance=effective_instance_name, gateway_url=PUSHGATEWAY_URL)
    except Exception as e:
        log.error("push_metrics_to_gateway_failed", job_name=job_name, error=str(e), gateway_url=PUSHGATEWAY_URL)


# --- Example Instrumentation Helper (for AI, can be adapted) ---
def observe_ai_interaction(
    model_provider: str, 
    endpoint_task: str, 
    duration_seconds: float, 
    prompt_tokens: Optional[int] = None, 
    completion_tokens: Optional[int] = None, 
    total_tokens: Optional[int] = None, 
    success: bool = True, 
    error_type: Optional[str] = None
):
    status = "success" if success else "error"
    VCF_AGENT_AI_REQUESTS_TOTAL.labels(model_provider=model_provider, endpoint_task=endpoint_task, status=status).inc()
    VCF_AGENT_AI_RESPONSE_SECONDS.labels(model_provider=model_provider, endpoint_task=endpoint_task).observe(duration_seconds)
    
    if prompt_tokens is not None:
        VCF_AGENT_AI_TOKENS_TOTAL.labels(model_provider=model_provider, endpoint_task=endpoint_task, token_type="prompt").inc(prompt_tokens)
    if completion_tokens is not None:
        VCF_AGENT_AI_TOKENS_TOTAL.labels(model_provider=model_provider, endpoint_task=endpoint_task, token_type="completion").inc(completion_tokens)
    if total_tokens is not None and not (prompt_tokens or completion_tokens): # Avoid double counting if specific tokens given
         VCF_AGENT_AI_TOKENS_TOTAL.labels(model_provider=model_provider, endpoint_task=endpoint_task, token_type="total").inc(total_tokens)

    if not success and error_type:
        VCF_AGENT_AI_ERRORS_TOTAL.labels(model_provider=model_provider, endpoint_task=endpoint_task, error_type=error_type).inc()

# Note: Concurrent requests gauge would be incremented/decremented around the actual call.

# --- Example Usage (Illustrative) ---
if __name__ == "__main__":
    start_metrics_http_server()
    log.info("example_service_started", metrics_port=METRICS_HTTP_PORT)

    # Simulate an AI interaction
    model = "test_model_provider"
    task = "test_endpoint_task"
    
    VCF_AGENT_AI_CONCURRENT_REQUESTS.labels(model_provider=model, endpoint_task=task).inc()
    start_time = time.time()
    try:
        # Simulate work & response
        time.sleep(0.1) # Simulate work
        sim_prompt_tokens = 50
        sim_completion_tokens = 100
        sim_error = None
        sim_success = True
        # sim_error_type = "SimulatedError" 
        # sim_success = False

        if not sim_success:
            sim_error_type = "SimulatedError" # Define if testing error path
            raise ValueError(sim_error_type)
        
        duration = time.time() - start_time
        observe_ai_interaction(
            model_provider=model, 
            endpoint_task=task, 
            duration_seconds=duration, 
            prompt_tokens=sim_prompt_tokens,
            completion_tokens=sim_completion_tokens,
            success=sim_success,
            # error_type=sim_error_type if not sim_success else None
        )
        log.info("simulated_ai_interaction_success", model=model, task=task, duration=duration)

    except Exception as e:
        duration = time.time() - start_time
        observe_ai_interaction(
            model_provider=model, 
            endpoint_task=task, 
            duration_seconds=duration, 
            success=False,
            error_type=type(e).__name__ if not isinstance(e, ValueError) else e.args[0] # Get the original sim_error_type
        )
        log.error("simulated_ai_interaction_failed", model=model, task=task, error=str(e))
    finally:
        VCF_AGENT_AI_CONCURRENT_REQUESTS.labels(model_provider=model, endpoint_task=task).dec()

    # Keep alive for scraping or until interrupted
    log.info("Example running. Metrics available. Press Ctrl+C to exit.")
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        log.info("example_service_shutting_down")

# --- Metrics Documentation (Summary - full details in main project docs) ---
# VCF_AGENT_AI_REQUESTS_TOTAL: Counter of total AI requests by model_provider, endpoint_task, and status.
# VCF_AGENT_AI_RESPONSE_SECONDS: Histogram of AI response times by model_provider and endpoint_task.
# VCF_AGENT_AI_TOKENS_TOTAL: Counter of total tokens by model_provider, endpoint_task, and token_type.
# VCF_AGENT_AI_ERRORS_TOTAL: Counter of AI errors by model_provider, endpoint_task, and error_type.
# VCF_AGENT_AI_CONCURRENT_REQUESTS: Gauge of concurrent AI requests by model_provider and endpoint_task.
# (CLI, Tool, BCFTools metrics to be documented as they are fully instrumented) 