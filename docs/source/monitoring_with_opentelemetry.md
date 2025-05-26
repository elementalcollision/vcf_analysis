# Distributed Tracing with OpenTelemetry

The VCF Agent uses OpenTelemetry for distributed tracing to provide observability across its components. This document explains how to use and extend the tracing capabilities.

## Architecture

The VCF Agent's tracing system consists of the following components:

1. **OpenTelemetry SDK** - Core tracing functionality
2. **OTLP Exporters** - Send traces to Jaeger
3. **Auto-Instrumentation** - Automatic tracing for HTTP requests, logging, and asyncio
4. **Custom Spans** - Application-specific span creation
5. **Jaeger** - Trace visualization and analysis
6. **Prometheus** - Metrics collection and monitoring
7. **Grafana** - Dashboards for visualizing metrics

## Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `OTEL_EXPORTER_OTLP_TRACES_ENDPOINT` | Jaeger collector endpoint | `http://localhost:4317` |
| `OTEL_EXPORTER_OTLP_TRACES_PROTOCOL` | Protocol (`grpc` or `http/protobuf`) | `grpc` |
| `OTEL_PYTHON_LOG_LEVEL` | OpenTelemetry logging level | Not set |
| `OTEL_SERVICE_NAME` | Override service name | Uses function parameter |

### Docker Compose Setup

The `docker-compose.yml` file includes Jaeger, Prometheus, and Grafana services:

```yaml
services:
  jaeger:
    image: jaegertracing/all-in-one:latest
    ports:
      - "4317:4317"  # OTLP gRPC receiver
      - "4318:4318"  # OTLP HTTP receiver
      - "16686:16686"  # Jaeger UI
      - "16687:16687"  # Prometheus metrics
    environment:
      - COLLECTOR_OTLP_ENABLED=true
      
  prometheus:
    image: prom/prometheus:latest
    ports:
      - "9090:9090"
    volumes:
      - ./prometheus.yml:/etc/prometheus/prometheus.yml
      
  grafana:
    image: grafana/grafana-oss:latest
    ports:
      - "3000:3000"
    volumes:
      - grafana_data:/var/lib/grafana
    depends_on:
      - prometheus
```

## Usage

### Initializing a Tracer

```python
from vcf_agent.tracing import init_tracer

# Initialize a tracer for your component
tracer = init_tracer(service_name="my-component")
```

### Creating Spans

```python
# Basic span
with tracer.start_as_current_span("operation_name") as span:
    span.set_attribute("attribute_name", "value")
    # Your code here
    
# Error handling
from opentelemetry.trace import Status, StatusCode

with tracer.start_as_current_span("risky_operation") as span:
    try:
        # Code that might fail
        result = perform_operation()
        span.set_attribute("operation.result", "success")
    except Exception as e:
        span.set_status(Status(StatusCode.ERROR, str(e)))
        span.set_attribute("error.type", type(e).__name__)
        span.set_attribute("error.message", str(e))
        raise
```

### Instrumenting Tools

The VCF Agent uses a consistent pattern for instrumenting tools:

```python
@tools.tool
def my_tool(arg1, arg2):
    with agent_tracer.start_as_current_span(f"tool.{tool_name}") as tool_span:
        tool_span.set_attribute("tool.name", tool_name)
        tool_span.set_attribute("tool.args.arg1", arg1)
        tool_span.set_attribute("tool.args.arg2", arg2)
        tool_span.add_event("tool.execution.start")
        
        try:
            # Tool implementation
            result = process_data(arg1, arg2)
            
            tool_span.set_attribute("result.success", True)
            tool_span.add_event("tool.execution.end")
            return result
        except Exception as e:
            tool_span.record_exception(e)
            tool_span.set_status(trace.StatusCode.ERROR, str(e))
            tool_span.add_event("tool.execution.failed")
            raise
```

## Components with Tracing

The VCF Agent implements tracing in several key components:

1. **CLI** - Traces command execution
2. **Agent Core** - Traces tool usage and LLM interactions
3. **BCFtools Integration** - Traces command execution and error handling
4. **LanceDB Operations** - Traces database operations
5. **Kuzu Graph DB** - Traces graph database operations

## Auto-Instrumentation

The agent automatically instruments:

- **HTTP Requests** - All outbound HTTP calls via the `requests` library
- **Logging** - Correlation between logs and traces
- **Asyncio** - Async/await operations

Auto-instrumentation is enabled once during application startup:

```python
from vcf_agent.tracing import setup_auto_instrumentation

setup_auto_instrumentation()
```

## Viewing Traces

1. Start the Docker Compose services:
   ```
   docker-compose up -d
   ```

2. Access the Jaeger UI at http://localhost:16686

3. Select a service (e.g., `vcf-agent-cli`, `vcf-agent-core`) and click "Find Traces"

## Troubleshooting

### No Traces in Jaeger

1. Check that Jaeger is running: `docker-compose ps`
2. Verify the OTLP endpoint is correct: `OTEL_EXPORTER_OTLP_TRACES_ENDPOINT=http://localhost:4317`
3. Enable debug logging: `OTEL_PYTHON_LOG_LEVEL=debug`
4. Check for errors in the application logs

### Missing Spans

1. Ensure `init_tracer()` is called before creating spans
2. Verify span creation is within an active trace context
3. Check that spans are properly closed (use context managers)

## Extending the Tracing System

### Adding a New Instrumented Component

1. Import the tracer and initialize it:
   ```python
   from vcf_agent.tracing import init_tracer
   
   component_tracer = init_tracer(service_name="vcf-agent-new-component")
   ```

2. Add spans to your operations:
   ```python
   with component_tracer.start_as_current_span("operation_name") as span:
       # Operation code
   ```

3. Include relevant attributes:
   ```python
   span.set_attribute("operation.type", "validation")
   span.set_attribute("file.path", filepath)
   ```

### Creating Custom Span Attributes

Follow these naming conventions for span attributes:

- **File Operations**: `file.path`, `file.size_bytes`, `file.type`
- **Database Operations**: `db.operation`, `db.table`, `db.query`
- **Errors**: `error.type`, `error.message`
- **Tools**: `tool.name`, `tool.args.<arg_name>`, `tool.status`
- **AI**: `llm.model_provider`, `llm.task`, `llm.status`

## Best Practices

1. **Use Descriptive Span Names**: Include the component and operation 
2. **Add Relevant Attributes**: Include enough context to understand the operation
3. **Handle Errors Properly**: Set status and record exceptions
4. **Close Spans Appropriately**: Use context managers
5. **Propagate Context**: Ensure context flows through async operations
6. **Keep Spans Focused**: Create spans for logical operations, not every function

## Performance Considerations

1. **Span Creation Overhead**: Creating spans has minimal overhead (microseconds)
2. **Attribute Cardinality**: Avoid high-cardinality values (e.g., user IDs, timestamps)
3. **Batch Processing**: The system uses BatchSpanProcessor for efficient exporting
4. **Sampling**: Consider sampling for high-throughput applications 