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