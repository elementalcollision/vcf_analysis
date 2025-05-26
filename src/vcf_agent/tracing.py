"""
OpenTelemetry Distributed Tracing for VCF Analysis Agent

This module provides comprehensive distributed tracing capabilities using OpenTelemetry,
enabling end-to-end observability of requests through the VCF Analysis Agent system.

## Key Features

- **Automatic Instrumentation**: HTTP requests, logging, and asyncio operations
- **Custom Spans**: Tool executions, AI model interactions, VCF processing operations
- **Trace Context Propagation**: Maintains trace context across service boundaries
- **Multiple Export Protocols**: Support for gRPC and HTTP/protobuf to Jaeger
- **Resource Attribution**: Service metadata and version information

## Architecture

The tracing system consists of:
- **TracerProvider**: Central configuration and span management
- **OTLP Exporter**: Sends traces to Jaeger collector
- **BatchSpanProcessor**: Efficient batching and export of spans
- **Auto-Instrumentation**: Automatic tracing of common operations
- **Custom Tracers**: Application-specific span creation

## Usage Examples

### Basic Tracer Setup
```python
from vcf_agent.tracing import init_tracer, setup_auto_instrumentation

# Initialize tracer for your service
tracer = init_tracer("vcf-processor")

# Enable automatic instrumentation
setup_auto_instrumentation()
```

### Creating Custom Spans
```python
from vcf_agent.tracing import init_tracer

tracer = init_tracer("vcf-analysis")

# Basic span creation
with tracer.start_as_current_span("validate_vcf_file") as span:
    span.set_attribute("file.path", "/path/to/file.vcf")
    span.set_attribute("file.size_bytes", 1024000)
    
    # Your validation logic here
    result = validate_file()
    
    span.set_attribute("validation.result", "valid")
    span.set_attribute("variant.count", result.variant_count)

# Nested spans for detailed tracing
with tracer.start_as_current_span("process_vcf_pipeline") as parent_span:
    parent_span.set_attribute("pipeline.version", "1.0")
    
    with tracer.start_as_current_span("parse_header") as header_span:
        header_span.set_attribute("header.sample_count", 10)
        # Parse header logic
    
    with tracer.start_as_current_span("process_variants") as variant_span:
        variant_span.set_attribute("processing.batch_size", 1000)
        # Process variants logic
```

### Error Handling and Status
```python
from opentelemetry.trace import Status, StatusCode

with tracer.start_as_current_span("risky_operation") as span:
    try:
        # Operation that might fail
        result = perform_operation()
        span.set_attribute("operation.result", "success")
    except Exception as e:
        span.set_status(Status(StatusCode.ERROR, str(e)))
        span.set_attribute("error.type", type(e).__name__)
        span.set_attribute("error.message", str(e))
        raise
```

### Tool Instrumentation Pattern
```python
def my_analysis_tool(input_data: str) -> str:
    \"\"\"Example of instrumenting an agent tool.\"\"\"
    tracer = init_tracer("vcf-tools")
    
    with tracer.start_as_current_span("tool.my_analysis") as span:
        span.set_attribute("tool.name", "my_analysis")
        span.set_attribute("tool.version", "1.0")
        span.set_attribute("input.size", len(input_data))
        
        try:
            result = process_data(input_data)
            span.set_attribute("output.size", len(result))
            span.set_attribute("tool.status", "success")
            return result
        except Exception as e:
            span.set_status(Status(StatusCode.ERROR, str(e)))
            span.set_attribute("tool.status", "error")
            span.set_attribute("error.type", type(e).__name__)
            raise
```

## Environment Variables

### Required Configuration
- `OTEL_EXPORTER_OTLP_TRACES_ENDPOINT`: Jaeger collector endpoint
  - gRPC: `http://localhost:4317` (default)
  - HTTP: `http://localhost:4318/v1/traces`

### Optional Configuration
- `OTEL_EXPORTER_OTLP_TRACES_PROTOCOL`: Export protocol (`grpc` or `http/protobuf`)
- `OTEL_EXPORTER_OTLP_TRACES_HEADERS`: Additional headers for authentication
- `OTEL_PYTHON_LOG_LEVEL`: Set to `debug` for detailed OTel logging
- `OTEL_SERVICE_NAME`: Override service name (defaults to function parameter)

## Auto-Instrumentation

The module automatically instruments:
- **HTTP Requests**: All outbound HTTP calls via `requests` library
- **Logging**: Correlation between logs and traces
- **Asyncio**: Async/await operations and task scheduling

Auto-instrumentation is enabled once per application lifecycle and provides
zero-code tracing for common operations.

## Best Practices

### Span Naming
- Use descriptive, hierarchical names: `service.operation.sub_operation`
- Include the operation type: `tool.validate_vcf`, `api.get_variant`
- Be consistent across similar operations

### Attributes
- Use semantic conventions where possible
- Include relevant context: file paths, IDs, counts, durations
- Avoid high-cardinality values (user IDs, timestamps)
- Use structured attribute names: `file.path`, `variant.count`

### Error Handling
- Always set span status for errors
- Include error type and message as attributes
- Don't let tracing exceptions break application flow

### Performance
- Spans have minimal overhead (~microseconds)
- Batch processing optimizes export performance
- Use sampling in high-throughput scenarios

## Integration with Metrics

Tracing integrates seamlessly with the metrics system:
- Trace context is automatically added to structured logs
- Metrics can be correlated with specific traces
- Both systems share service metadata and configuration

## Troubleshooting

### No Traces in Jaeger
1. Check Jaeger collector is running on the configured endpoint
2. Verify network connectivity from agent to Jaeger
3. Enable debug logging: `OTEL_PYTHON_LOG_LEVEL=debug`
4. Check for export errors in application logs

### Missing Spans
1. Ensure `init_tracer()` is called before creating spans
2. Verify span creation is within an active trace context
3. Check that spans are properly closed (use context managers)

### Performance Issues
1. Adjust batch processor settings for your workload
2. Consider sampling for high-volume applications
3. Monitor export queue depth and processing times

## Example: Complete Service Setup
```python
from vcf_agent.tracing import init_tracer, setup_auto_instrumentation

def main():
    # Initialize tracing for the service
    tracer = init_tracer("vcf-analysis-service")
    setup_auto_instrumentation()
    
    # Your application code with automatic tracing
    with tracer.start_as_current_span("main_operation"):
        # All HTTP requests, logs, and async operations
        # will be automatically traced
        process_vcf_files()

if __name__ == "__main__":
    main()
```

For more information, see the observability documentation at:
docs/source/monitoring_with_prometheus.md
"""

import os
import importlib.metadata # Added for versioning
from opentelemetry import trace
from opentelemetry.sdk.trace import TracerProvider as SdkTracerProvider # Explicit SDK type
from opentelemetry.sdk.trace.export import BatchSpanProcessor

# Pick the OTLP exporter based on protocol environment variable if set
# Defaulting to gRPC if not specified or if http is not explicitly requested.
_otel_protocol = os.getenv("OTEL_EXPORTER_OTLP_TRACES_PROTOCOL", "grpc")

if _otel_protocol == "http/protobuf":
    from opentelemetry.exporter.otlp.proto.http.trace_exporter import OTLPSpanExporter
else:  # grpc is the default
    from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import OTLPSpanExporter

from opentelemetry.sdk.resources import Resource, SERVICE_NAME, TELEMETRY_SDK_LANGUAGE, TELEMETRY_SDK_NAME, TELEMETRY_SDK_VERSION
from opentelemetry.instrumentation.requests import RequestsInstrumentor
from opentelemetry.instrumentation.logging import LoggingInstrumentor
from opentelemetry.instrumentation.asyncio import AsyncioInstrumentor
# Ensure this import works relative to the `src` directory if `vcf_agent` is a top-level package in `src`
# If tracing.py is inside vcf_agent, then it should be `from . import __version__`
try:
    from vcf_agent import __version__ # Assuming __version__ is in vcf_agent/__init__.py
except ImportError:
    from . import __version__ # Fallback for relative import if tracing.py is inside vcf_agent package

# Configure root logger if OTEL_PYTHON_LOG_LEVEL is debug
_otel_python_log_level_env = os.getenv("OTEL_PYTHON_LOG_LEVEL")
if _otel_python_log_level_env and _otel_python_log_level_env.lower() == "debug":
    import logging
    # Configure the root logger
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        force=True # force=True will remove and re-add handlers to the root logger
    )
    print(f"[TRACING] Root logger configured to DEBUG because OTEL_PYTHON_LOG_LEVEL is set to '{_otel_python_log_level_env}'")

# Flag to ensure auto-instrumentation setup runs only once
_AUTO_INSTRUMENTATION_SETUP_DONE = False

def init_tracer(service_name: str) -> trace.Tracer:
    current_provider = trace.get_tracer_provider()
    is_sdk_provider_already_set = isinstance(current_provider, SdkTracerProvider)

    if not is_sdk_provider_already_set:
        resource = Resource(attributes={
            SERVICE_NAME: service_name,
            TELEMETRY_SDK_LANGUAGE: "python",
            TELEMETRY_SDK_NAME: "opentelemetry",
            TELEMETRY_SDK_VERSION: importlib.metadata.version('opentelemetry-api'), # Corrected SDK version
            "service.version": __version__,
        })
        
        # Explicitly create SdkTracerProvider
        provider = SdkTracerProvider(resource=resource) 
        
        exporter_endpoint = os.getenv("OTEL_EXPORTER_OTLP_TRACES_ENDPOINT")
        if not exporter_endpoint and _otel_protocol == "grpc":
            exporter_endpoint = "http://localhost:4317" # Default for gRPC
            print(f"[TRACING] WARNING: OTEL_EXPORTER_OTLP_TRACES_ENDPOINT not set. Defaulting to {exporter_endpoint} for gRPC protocol.")
        elif not exporter_endpoint and _otel_protocol == "http/protobuf":
            exporter_endpoint = "http://localhost:4318" # Default for HTTP/protobuf
            print(f"[TRACING] WARNING: OTEL_EXPORTER_OTLP_TRACES_ENDPOINT not set. Defaulting to {exporter_endpoint} for http/protobuf protocol.")
        elif not exporter_endpoint:
            print(f"[TRACING] WARNING: OTEL_EXPORTER_OTLP_TRACES_ENDPOINT not set and protocol is '{_otel_protocol}'. Exporter might not work.")

        print(f"[TRACING] OTLP Exporter Endpoint: {exporter_endpoint}")
        print(f"[TRACING] OTLP Exporter Protocol: {_otel_protocol}")

        otlp_exporter = OTLPSpanExporter(
            endpoint=exporter_endpoint,
            # headers can be set via OTEL_EXPORTER_OTLP_TRACES_HEADERS env var
            # protocol is handled by choice of exporter class above
        )
        provider.add_span_processor(BatchSpanProcessor(
            otlp_exporter,
            schedule_delay_millis=200, # Default 5000
            export_timeout_millis=10000 # Default 30000, but let's try a bit longer for flush
        ))
        
        # This call is guarded by Once() inside OTel 
        # and will issue the warning if overridden without a good reason.
        # Our check `if not is_sdk_provider_already_set` prevents unnecessary reconfiguration.
        trace.set_tracer_provider(provider) 
        print(f"OpenTelemetry TracerProvider SET by init_tracer for service: {service_name}")
    else:
        print(f"OpenTelemetry TracerProvider was already SET. init_tracer for '{service_name}' will use existing provider. Type: {type(current_provider)}")

    return trace.get_tracer(
        instrumenting_module_name=f"vcf_agent.{service_name.replace('-', '_')}", 
        instrumenting_library_version=__version__ # Corrected keyword
    )


def setup_auto_instrumentation():
    global _AUTO_INSTRUMENTATION_SETUP_DONE

    if _AUTO_INSTRUMENTATION_SETUP_DONE:
        print("Auto-instrumentation already configured.")
        return

    current_provider = trace.get_tracer_provider()
    if not isinstance(current_provider, SdkTracerProvider):
        # This ensures that if setup_auto_instrumentation is called before any init_tracer,
        # a provider is initialized.
        print("Warning: SDK TracerProvider not set before setup_auto_instrumentation. Initializing with default service name.")
        init_tracer("vcf-agent-default-setup")

    RequestsInstrumentor().instrument()
    print("Requests auto-instrumentation enabled.")
    
    LoggingInstrumentor().instrument(set_logging_format=False)
    print("Logging auto-instrumentation enabled.")

    AsyncioInstrumentor().instrument()
    print("Asyncio auto-instrumentation enabled.")
    
    _AUTO_INSTRUMENTATION_SETUP_DONE = True

# Example usage (typically called once at application startup):
# if __name__ == "__main__":
#     tracer = init_tracer(service_name="my-python-app")
#     setup_auto_instrumentation()
#
#     with tracer.start_as_current_span("example-main-operation") as main_span:
#         main_span.set_attribute("example.attribute", "example_value")
#         print("Doing some work...")
#         # Your application code here
#     print("Tracing initialized and example span created.") 