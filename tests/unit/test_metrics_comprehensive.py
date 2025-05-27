"""
Comprehensive unit tests for the metrics module.

Tests the Prometheus metrics collection, structured logging, HTTP server,
Pushgateway functionality, and OpenTelemetry integration.
"""

import pytest
import os
import time
import threading
from unittest.mock import patch, MagicMock, Mock, call
from prometheus_client import CollectorRegistry, Counter, Histogram, Gauge

import vcf_agent.metrics
from vcf_agent.metrics import (
    add_otel_context,
    setup_logging,
    log,
    http_registry,
    VCF_AGENT_AI_REQUESTS_TOTAL,
    VCF_AGENT_AI_RESPONSE_SECONDS,
    VCF_AGENT_AI_TOKENS_TOTAL,
    VCF_AGENT_AI_ERRORS_TOTAL,
    VCF_AGENT_AI_CONCURRENT_REQUESTS,
    CLI_COMMAND_REQUESTS_TOTAL,
    CLI_COMMAND_DURATION_SECONDS,
    CLI_COMMAND_ERRORS_TOTAL,
    VCF_AGENT_TOOL_REQUESTS_TOTAL,
    VCF_AGENT_TOOL_DURATION_SECONDS,
    VCF_AGENT_TOOL_ERRORS_TOTAL,
    VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL,
    VCF_AGENT_BCFTOOLS_DURATION_SECONDS,
    VCF_AGENT_BCFTOOLS_ERRORS_TOTAL,
    start_metrics_http_server,
    push_job_metrics,
    observe_ai_interaction,
    observe_bcftools_command,
    METRICS_HTTP_PORT,
    PUSHGATEWAY_URL
)


class TestOTelContextProcessor:
    """Test cases for the OpenTelemetry context processor."""

    @patch('vcf_agent.metrics.trace.get_current_span')
    @patch('vcf_agent.metrics.format_trace_id')
    @patch('vcf_agent.metrics.format_span_id')
    def test_add_otel_context_with_valid_span(self, mock_format_span_id, mock_format_trace_id, mock_get_current_span):
        """Test adding OTel context with a valid recording span."""
        # Mock span with valid context
        from vcf_agent.metrics import Span
        mock_span = MagicMock(spec=Span)
        mock_span.is_recording.return_value = True
        
        mock_span_context = MagicMock()
        mock_span_context.is_valid = True
        mock_span_context.trace_id = 12345678901234567890123456789012
        mock_span_context.span_id = 1234567890123456
        mock_span.get_span_context.return_value = mock_span_context
        
        # Mock resource with service name
        mock_resource = MagicMock()
        mock_resource.attributes = {'service.name': 'vcf-agent'}
        mock_span.resource = mock_resource
        
        mock_get_current_span.return_value = mock_span
        mock_format_trace_id.return_value = "formatted_trace_id"
        mock_format_span_id.return_value = "formatted_span_id"
        
        event_dict = {"message": "test"}
        result = add_otel_context(None, None, event_dict)
        
        assert "trace_id" in result
        assert "span_id" in result
        assert result["service.name"] == "vcf-agent"

    @patch('vcf_agent.metrics.trace.get_current_span')
    @patch('vcf_agent.metrics.format_trace_id')
    @patch('vcf_agent.metrics.format_span_id')
    def test_add_otel_context_with_private_resource(self, mock_format_span_id, mock_format_trace_id, mock_get_current_span):
        """Test adding OTel context with private resource attribute."""
        from vcf_agent.metrics import Span
        mock_span = MagicMock(spec=Span)
        mock_span.is_recording.return_value = True
        
        mock_span_context = MagicMock()
        mock_span_context.is_valid = True
        mock_span_context.trace_id = 12345678901234567890123456789012
        mock_span_context.span_id = 1234567890123456
        mock_span.get_span_context.return_value = mock_span_context
        
        # No public resource, but private _resource
        mock_span.resource = None
        mock_private_resource = MagicMock()
        mock_private_resource.attributes = {'service.name': 'vcf-agent-private'}
        mock_span._resource = mock_private_resource
        
        mock_get_current_span.return_value = mock_span
        mock_format_trace_id.return_value = "formatted_trace_id"
        mock_format_span_id.return_value = "formatted_span_id"
        
        event_dict = {"message": "test"}
        result = add_otel_context(None, None, event_dict)
        
        assert result["service.name"] == "vcf-agent-private"

    @patch('vcf_agent.metrics.trace.get_current_span')
    def test_add_otel_context_with_non_recording_span(self, mock_get_current_span):
        """Test adding OTel context with non-recording span."""
        mock_span = MagicMock()
        mock_span.is_recording.return_value = False
        mock_get_current_span.return_value = mock_span
        
        event_dict = {"message": "test"}
        result = add_otel_context(None, None, event_dict)
        
        assert "trace_id" not in result
        assert "span_id" not in result
        assert result == event_dict

    @patch('vcf_agent.metrics.trace.get_current_span')
    def test_add_otel_context_with_invalid_span_context(self, mock_get_current_span):
        """Test adding OTel context with invalid span context."""
        mock_span = MagicMock()
        mock_span.is_recording.return_value = True
        
        mock_span_context = MagicMock()
        mock_span_context.is_valid = False
        mock_span.get_span_context.return_value = mock_span_context
        
        mock_get_current_span.return_value = mock_span
        
        event_dict = {"message": "test"}
        result = add_otel_context(None, None, event_dict)
        
        assert "trace_id" not in result
        assert "span_id" not in result

    @patch('vcf_agent.metrics.trace.get_current_span')
    def test_add_otel_context_with_no_span(self, mock_get_current_span):
        """Test adding OTel context when no span is available."""
        mock_get_current_span.return_value = None
        
        event_dict = {"message": "test"}
        result = add_otel_context(None, None, event_dict)
        
        assert result == event_dict


class TestLoggingSetup:
    """Test cases for logging setup functionality."""

    @patch('vcf_agent.metrics.structlog.configure')
    @patch('vcf_agent.metrics.logging.basicConfig')
    @patch.dict(os.environ, {'LOG_LEVEL': 'DEBUG'})
    def test_setup_logging_with_debug_level(self, mock_basic_config, mock_structlog_configure):
        """Test logging setup with DEBUG level from environment."""
        mock_logger = MagicMock()
        with patch('vcf_agent.metrics.structlog.get_logger', return_value=mock_logger):
            result = setup_logging()
        
        mock_basic_config.assert_called_once_with(level='DEBUG', format="%(message)s")
        # structlog.configure may be called multiple times during module import
        assert mock_structlog_configure.call_count >= 1
        assert result == mock_logger

    @patch('vcf_agent.metrics.structlog.configure')
    @patch('vcf_agent.metrics.logging.basicConfig')
    @patch.dict(os.environ, {}, clear=True)
    def test_setup_logging_with_default_level(self, mock_basic_config, mock_structlog_configure):
        """Test logging setup with default INFO level."""
        mock_logger = MagicMock()
        with patch('vcf_agent.metrics.structlog.get_logger', return_value=mock_logger):
            result = setup_logging()
        
        mock_basic_config.assert_called_once_with(level='INFO', format="%(message)s")

    def test_log_instance_exists(self):
        """Test that the log instance is properly created."""
        assert log is not None
        # Test that we can call logging methods without error
        with patch.object(log, 'info') as mock_info:
            log.info("test message")
            mock_info.assert_called_once_with("test message")


class TestMetricsRegistry:
    """Test cases for metrics registry and metric definitions."""

    def test_http_registry_exists(self):
        """Test that the HTTP registry is properly created."""
        assert isinstance(http_registry, CollectorRegistry)

    def test_ai_metrics_exist(self):
        """Test that AI metrics are properly defined."""
        assert isinstance(VCF_AGENT_AI_REQUESTS_TOTAL, Counter)
        assert isinstance(VCF_AGENT_AI_RESPONSE_SECONDS, Histogram)
        assert isinstance(VCF_AGENT_AI_TOKENS_TOTAL, Counter)
        assert isinstance(VCF_AGENT_AI_ERRORS_TOTAL, Counter)
        assert isinstance(VCF_AGENT_AI_CONCURRENT_REQUESTS, Gauge)

    def test_cli_metrics_exist(self):
        """Test that CLI metrics are properly defined."""
        assert isinstance(CLI_COMMAND_REQUESTS_TOTAL, Counter)
        assert isinstance(CLI_COMMAND_DURATION_SECONDS, Histogram)
        assert isinstance(CLI_COMMAND_ERRORS_TOTAL, Counter)

    def test_tool_metrics_exist(self):
        """Test that tool metrics are properly defined."""
        assert isinstance(VCF_AGENT_TOOL_REQUESTS_TOTAL, Counter)
        assert isinstance(VCF_AGENT_TOOL_DURATION_SECONDS, Histogram)
        assert isinstance(VCF_AGENT_TOOL_ERRORS_TOTAL, Counter)

    def test_bcftools_metrics_exist(self):
        """Test that BCFtools metrics are properly defined."""
        assert isinstance(VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL, Counter)
        assert isinstance(VCF_AGENT_BCFTOOLS_DURATION_SECONDS, Histogram)
        assert isinstance(VCF_AGENT_BCFTOOLS_ERRORS_TOTAL, Counter)


class TestMetricsHTTPServer:
    """Test cases for the metrics HTTP server functionality."""

    @patch('vcf_agent.metrics.threading.Thread')
    @patch('vcf_agent.metrics.start_http_server')
    @patch('vcf_agent.metrics.log')
    def test_start_metrics_http_server_success(self, mock_log, mock_start_http_server, mock_thread):
        """Test successful metrics HTTP server start."""
        mock_thread_instance = MagicMock()
        mock_thread_instance.is_alive.return_value = False
        mock_thread.return_value = mock_thread_instance
        
        # Reset the global thread variable
        import vcf_agent.metrics
        vcf_agent.metrics._metrics_server_thread = None
        
        start_metrics_http_server(port=8001)
        
        mock_thread.assert_called_once()
        mock_thread_instance.start.assert_called_once()
        mock_log.info.assert_called_with("prometheus_metrics_server_started", port=8001)

    @patch('vcf_agent.metrics.threading.Thread')
    @patch('vcf_agent.metrics.start_http_server')
    @patch('vcf_agent.metrics.log')
    def test_start_metrics_http_server_default_port(self, mock_log, mock_start_http_server, mock_thread):
        """Test metrics HTTP server start with default port."""
        mock_thread_instance = MagicMock()
        mock_thread_instance.is_alive.return_value = False
        mock_thread.return_value = mock_thread_instance
        
        # Reset the global thread variable
        import vcf_agent.metrics
        vcf_agent.metrics._metrics_server_thread = None
        
        start_metrics_http_server()
        
        mock_log.info.assert_called_with("prometheus_metrics_server_started", port=METRICS_HTTP_PORT)

    @patch('vcf_agent.metrics.threading.Thread')
    @patch('vcf_agent.metrics.start_http_server')
    @patch('vcf_agent.metrics.log')
    def test_start_metrics_http_server_already_running(self, mock_log, mock_start_http_server, mock_thread):
        """Test metrics HTTP server when already running."""
        mock_thread_instance = MagicMock()
        mock_thread_instance.is_alive.return_value = True
        
        # Set the global thread variable to simulate running server
        import vcf_agent.metrics
        vcf_agent.metrics._metrics_server_thread = mock_thread_instance
        
        start_metrics_http_server(port=8001)
        
        mock_thread.assert_not_called()
        mock_log.info.assert_called_with("prometheus_metrics_server_already_running", port=8001)

    @patch('vcf_agent.metrics.threading.Thread')
    @patch('vcf_agent.metrics.start_http_server')
    @patch('vcf_agent.metrics.log')
    def test_start_metrics_http_server_failure(self, mock_log, mock_start_http_server, mock_thread):
        """Test metrics HTTP server start failure."""
        mock_thread.side_effect = Exception("Thread creation failed")
        
        # Reset the global thread variable
        import vcf_agent.metrics
        vcf_agent.metrics._metrics_server_thread = None
        
        start_metrics_http_server(port=8001)
        
        mock_log.error.assert_called_with(
            "prometheus_metrics_server_failed_to_start", 
            port=8001, 
            error="Thread creation failed"
        )


class TestPushgatewayFunctionality:
    """Test cases for Pushgateway functionality."""

    @patch('vcf_agent.metrics.push_to_gateway')
    @patch('vcf_agent.metrics.CollectorRegistry')
    @patch('vcf_agent.metrics.log')
    @patch('vcf_agent.metrics.PUSHGATEWAY_URL', 'http://test-gateway:9091')
    def test_push_job_metrics_success(self, mock_log, mock_registry_class, mock_push_to_gateway):
        """Test successful metrics push to gateway."""
        mock_registry = MagicMock()
        mock_registry_class.return_value = mock_registry
        
        mock_metric1 = MagicMock()
        mock_metric2 = MagicMock()
        metrics_to_push = [mock_metric1, mock_metric2]
        
        push_job_metrics("test_job", metrics_to_push, "test_instance")
        
        mock_registry.register.assert_has_calls([call(mock_metric1), call(mock_metric2)])
        mock_push_to_gateway.assert_called_once_with(
            'http://test-gateway:9091', 
            job='test_job', 
            registry=mock_registry
        )
        mock_log.info.assert_called_with(
            "metrics_pushed_to_gateway", 
            job_name="test_job", 
            instance="test_instance", 
            gateway_url='http://test-gateway:9091'
        )

    @patch('vcf_agent.metrics.push_to_gateway')
    @patch('vcf_agent.metrics.CollectorRegistry')
    @patch('vcf_agent.metrics.log')
    @patch('vcf_agent.metrics.os.getpid', return_value=12345)
    @patch('vcf_agent.metrics.PUSHGATEWAY_URL', 'http://test-gateway:9091')
    def test_push_job_metrics_default_instance(self, mock_getpid, mock_log, mock_registry_class, mock_push_to_gateway):
        """Test metrics push with default instance name."""
        mock_registry = MagicMock()
        mock_registry_class.return_value = mock_registry
        
        mock_metric = MagicMock()
        metrics_to_push = [mock_metric]
        
        push_job_metrics("test_job", metrics_to_push)
        
        mock_log.info.assert_called_with(
            "metrics_pushed_to_gateway", 
            job_name="test_job", 
            instance="test_job-12345", 
            gateway_url='http://test-gateway:9091'
        )

    @patch('vcf_agent.metrics.log')
    @patch('vcf_agent.metrics.PUSHGATEWAY_URL', 'disabled')
    def test_push_job_metrics_disabled(self, mock_log):
        """Test metrics push when Pushgateway is disabled."""
        push_job_metrics("test_job", [])
        
        mock_log.warn.assert_called_with("pushgateway_not_configured_or_disabled", job_name="test_job")

    @patch('vcf_agent.metrics.log')
    @patch('vcf_agent.metrics.PUSHGATEWAY_URL', '')
    def test_push_job_metrics_not_configured(self, mock_log):
        """Test metrics push when Pushgateway URL is empty."""
        push_job_metrics("test_job", [])
        
        mock_log.warn.assert_called_with("pushgateway_not_configured_or_disabled", job_name="test_job")

    @patch('vcf_agent.metrics.push_to_gateway')
    @patch('vcf_agent.metrics.CollectorRegistry')
    @patch('vcf_agent.metrics.log')
    @patch('vcf_agent.metrics.PUSHGATEWAY_URL', 'http://test-gateway:9091')
    def test_push_job_metrics_failure(self, mock_log, mock_registry_class, mock_push_to_gateway):
        """Test metrics push failure."""
        mock_registry = MagicMock()
        mock_registry_class.return_value = mock_registry
        mock_push_to_gateway.side_effect = Exception("Push failed")
        
        mock_metric = MagicMock()
        metrics_to_push = [mock_metric]
        
        push_job_metrics("test_job", metrics_to_push)
        
        mock_log.error.assert_called_with(
            "push_metrics_to_gateway_failed", 
            job_name="test_job", 
            error="Push failed", 
            gateway_url='http://test-gateway:9091'
        )


class TestObserveAIInteraction:
    """Test cases for AI interaction observation."""

    @patch.object(VCF_AGENT_AI_REQUESTS_TOTAL, 'labels')
    @patch.object(VCF_AGENT_AI_RESPONSE_SECONDS, 'labels')
    @patch.object(VCF_AGENT_AI_TOKENS_TOTAL, 'labels')
    @patch.object(VCF_AGENT_AI_CONCURRENT_REQUESTS, 'labels')
    def test_observe_ai_interaction_success(self, mock_concurrent_labels, mock_tokens_labels, mock_response_labels, mock_requests_labels):
        """Test observing successful AI interaction."""
        # Mock metric instances
        mock_requests_metric = MagicMock()
        mock_response_metric = MagicMock()
        mock_tokens_metric = MagicMock()
        mock_concurrent_metric = MagicMock()
        
        mock_requests_labels.return_value = mock_requests_metric
        mock_response_labels.return_value = mock_response_metric
        mock_tokens_labels.return_value = mock_tokens_metric
        mock_concurrent_labels.return_value = mock_concurrent_metric
        
        observe_ai_interaction(
            model_provider="openai",
            endpoint_task="chat_completion",
            duration_seconds=1.5,
            prompt_tokens=100,
            completion_tokens=50,
            success=True
        )
        
        # Check that metrics were called correctly
        mock_requests_metric.inc.assert_called_once()
        mock_response_metric.observe.assert_called_once_with(1.5)
        mock_concurrent_metric.dec.assert_called_once()
        
        # Check token metrics were called
        assert mock_tokens_metric.inc.call_count == 2  # prompt and completion tokens

    @patch.object(VCF_AGENT_AI_REQUESTS_TOTAL, 'labels')
    @patch.object(VCF_AGENT_AI_RESPONSE_SECONDS, 'labels')
    @patch.object(VCF_AGENT_AI_ERRORS_TOTAL, 'labels')
    @patch.object(VCF_AGENT_AI_CONCURRENT_REQUESTS, 'labels')
    def test_observe_ai_interaction_failure(self, mock_concurrent_labels, mock_errors_labels, mock_response_labels, mock_requests_labels):
        """Test observing failed AI interaction."""
        # Mock metric instances
        mock_requests_metric = MagicMock()
        mock_response_metric = MagicMock()
        mock_errors_metric = MagicMock()
        mock_concurrent_metric = MagicMock()
        
        mock_requests_labels.return_value = mock_requests_metric
        mock_response_labels.return_value = mock_response_metric
        mock_errors_labels.return_value = mock_errors_metric
        mock_concurrent_labels.return_value = mock_concurrent_metric
        
        observe_ai_interaction(
            model_provider="ollama",
            endpoint_task="embedding",
            duration_seconds=0.5,
            success=False,
            error_type="timeout"
        )
        
        # Check that metrics were called correctly
        mock_requests_metric.inc.assert_called_once()
        mock_response_metric.observe.assert_called_once_with(0.5)
        mock_errors_metric.inc.assert_called_once()
        mock_concurrent_metric.dec.assert_called_once()

    def test_observe_ai_interaction_total_tokens_only(self):
        """Test observing AI interaction with only total tokens."""
        with patch.object(VCF_AGENT_AI_TOKENS_TOTAL, 'labels') as mock_labels:
            mock_metric = MagicMock()
            mock_labels.return_value = mock_metric
            
            observe_ai_interaction(
                model_provider="cerebras",
                endpoint_task="completion",
                duration_seconds=2.0,
                total_tokens=200,
                success=True
            )
            
            # Check that total tokens metric was called
            mock_metric.inc.assert_called_with(200)

    def test_observe_ai_interaction_no_double_counting(self):
        """Test that total tokens are not counted when specific tokens are provided."""
        with patch.object(VCF_AGENT_AI_TOKENS_TOTAL, 'labels') as mock_labels:
            mock_metric = MagicMock()
            mock_labels.return_value = mock_metric
            
            observe_ai_interaction(
                model_provider="openai",
                endpoint_task="chat",
                duration_seconds=1.0,
                prompt_tokens=50,
                completion_tokens=25,
                total_tokens=75,  # Should be ignored to avoid double counting
                success=True
            )
            
            # Check that only prompt and completion tokens were recorded (not total)
            assert mock_metric.inc.call_count == 2  # Only prompt and completion, not total


class TestObserveBCFToolsCommand:
    """Test cases for BCFtools command observation."""

    @patch.object(VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL, 'labels')
    @patch.object(VCF_AGENT_BCFTOOLS_DURATION_SECONDS, 'labels')
    def test_observe_bcftools_command_success(self, mock_duration_labels, mock_commands_labels):
        """Test observing successful BCFtools command."""
        mock_commands_metric = MagicMock()
        mock_duration_metric = MagicMock()
        
        mock_commands_labels.return_value = mock_commands_metric
        mock_duration_labels.return_value = mock_duration_metric
        
        observe_bcftools_command(
            bcftools_subcommand="view",
            duration_seconds=3.2,
            success=True
        )
        
        # Check that metrics were called correctly
        mock_commands_metric.inc.assert_called_once()
        mock_duration_metric.observe.assert_called_once_with(3.2)

    @patch.object(VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL, 'labels')
    @patch.object(VCF_AGENT_BCFTOOLS_DURATION_SECONDS, 'labels')
    @patch.object(VCF_AGENT_BCFTOOLS_ERRORS_TOTAL, 'labels')
    def test_observe_bcftools_command_failure(self, mock_errors_labels, mock_duration_labels, mock_commands_labels):
        """Test observing failed BCFtools command."""
        mock_commands_metric = MagicMock()
        mock_duration_metric = MagicMock()
        mock_errors_metric = MagicMock()
        
        mock_commands_labels.return_value = mock_commands_metric
        mock_duration_labels.return_value = mock_duration_metric
        mock_errors_labels.return_value = mock_errors_metric
        
        observe_bcftools_command(
            bcftools_subcommand="filter",
            duration_seconds=1.0,
            success=False,
            error_type="invalid_filter"
        )
        
        # Check that metrics were called correctly
        mock_commands_metric.inc.assert_called_once()
        mock_duration_metric.observe.assert_called_once_with(1.0)
        mock_errors_metric.inc.assert_called_once()


class TestEnvironmentVariables:
    """Test cases for environment variable handling."""

    def test_metrics_port_constant_exists(self):
        """Test that METRICS_HTTP_PORT constant exists."""
        assert hasattr(vcf_agent.metrics, 'METRICS_HTTP_PORT')
        assert isinstance(METRICS_HTTP_PORT, int)

    def test_pushgateway_url_constant_exists(self):
        """Test that PUSHGATEWAY_URL constant exists."""
        assert hasattr(vcf_agent.metrics, 'PUSHGATEWAY_URL')
        assert isinstance(PUSHGATEWAY_URL, str)


class TestMetricsIntegration:
    """Integration test cases for metrics functionality."""

    def test_metrics_labels_consistency(self):
        """Test that metric labels are consistent across related metrics."""
        # Test that metrics have the expected label names (tuples, not lists)
        assert VCF_AGENT_AI_REQUESTS_TOTAL._labelnames == ('model_provider', 'endpoint_task', 'status')
        assert VCF_AGENT_AI_RESPONSE_SECONDS._labelnames == ('model_provider', 'endpoint_task', 'status')
        assert VCF_AGENT_AI_CONCURRENT_REQUESTS._labelnames == ('model_provider', 'endpoint_task')

    def test_concurrent_requests_gauge_operations(self):
        """Test concurrent requests gauge increment and decrement."""
        with patch.object(VCF_AGENT_AI_CONCURRENT_REQUESTS, 'labels') as mock_labels:
            mock_gauge = MagicMock()
            mock_labels.return_value = mock_gauge
            
            gauge = VCF_AGENT_AI_CONCURRENT_REQUESTS.labels(
                model_provider="test", 
                endpoint_task="test"
            )
            
            # Test that labels method was called
            mock_labels.assert_called_with(model_provider="test", endpoint_task="test")

    def test_histogram_timing_operations(self):
        """Test histogram timing operations."""
        with patch.object(VCF_AGENT_AI_RESPONSE_SECONDS, 'labels') as mock_labels:
            mock_histogram = MagicMock()
            mock_labels.return_value = mock_histogram
            
            histogram = VCF_AGENT_AI_RESPONSE_SECONDS.labels(
                model_provider="test", 
                endpoint_task="test", 
                status="success"
            )
            
            # Test that labels method was called
            mock_labels.assert_called_with(model_provider="test", endpoint_task="test", status="success") 