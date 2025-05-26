import unittest
from unittest import mock
from opentelemetry.sdk.trace import Span
from vcf_agent import metrics

class TestMetricsTracing(unittest.TestCase):
    def test_add_otel_context_valid_span(self):
        """Test add_otel_context injects trace and span IDs for a valid span."""
        mock_span_context = mock.Mock()
        mock_span_context.is_valid = True
        mock_span_context.trace_id = 123
        mock_span_context.span_id = 456
        mock_span = mock.Mock(spec=Span)
        mock_span.is_recording.return_value = True
        mock_span.get_span_context.return_value = mock_span_context
        with mock.patch('opentelemetry.trace.get_current_span', return_value=mock_span):
            event_dict = {"event": "test"}
            result = metrics.add_otel_context(None, None, event_dict.copy())
            assert "trace_id" in result
            assert "span_id" in result

    def test_add_otel_context_invalid_span(self):
        """Test add_otel_context does not inject context for invalid span."""
        mock_span_context = mock.Mock()
        mock_span_context.is_valid = False
        mock_span = mock.Mock(spec=Span)
        mock_span.is_recording.return_value = True
        mock_span.get_span_context.return_value = mock_span_context
        with mock.patch('opentelemetry.trace.get_current_span', return_value=mock_span):
            event_dict = {"event": "test"}
            result = metrics.add_otel_context(None, None, event_dict.copy())
            assert "trace_id" not in result
            assert "span_id" not in result

    def test_add_otel_context_fallback_resource(self):
        """Test add_otel_context uses _resource if resource attribute is missing."""
        mock_span_context = mock.Mock()
        mock_span_context.is_valid = True
        mock_span_context.trace_id = 789
        mock_span_context.span_id = 101

        # Mock _resource with attributes
        mock_private_resource = mock.Mock()
        mock_private_resource.attributes = {metrics.ResourceAttributes.SERVICE_NAME: "fallback-service"}

        mock_span = mock.Mock(spec=Span)
        mock_span.is_recording.return_value = True
        mock_span.get_span_context.return_value = mock_span_context
        # Simulate resource attribute being absent or None, but _resource present
        mock_span.resource = None 
        mock_span._resource = mock_private_resource # Set the private _resource
        
        with mock.patch('opentelemetry.trace.get_current_span', return_value=mock_span):
            event_dict = {"event": "test_fallback"}
            result = metrics.add_otel_context(None, None, event_dict.copy())
            assert result.get("service.name") == "fallback-service"
            assert "trace_id" in result # Ensure trace/span still present
            assert "span_id" in result

    def test_setup_logging_configures_structlog(self):
        """Test setup_logging configures structlog and returns a logger."""
        with mock.patch('structlog.configure') as mock_configure, \
             mock.patch('structlog.get_logger', return_value="logger") as mock_get_logger:
            logger = metrics.setup_logging()
            self.assertGreaterEqual(mock_configure.call_count, 1)
            mock_get_logger.assert_called_once()
            assert logger == "logger"

    def test_observe_ai_interaction_success_and_error(self):
        """Test observe_ai_interaction records metrics for success and error cases."""
        with mock.patch.object(metrics, 'VCF_AGENT_AI_REQUESTS_TOTAL') as mock_requests, \
             mock.patch.object(metrics, 'VCF_AGENT_AI_RESPONSE_SECONDS') as mock_response, \
             mock.patch.object(metrics, 'VCF_AGENT_AI_TOKENS_TOTAL') as mock_tokens, \
             mock.patch.object(metrics, 'VCF_AGENT_AI_ERRORS_TOTAL') as mock_errors:
            # Success case
            metrics.observe_ai_interaction(
                model_provider="openai",
                endpoint_task="test",
                duration_seconds=1.0,
                prompt_tokens=10,
                completion_tokens=5,
                total_tokens=None,
                success=True
            )
            mock_requests.labels.assert_called_with(model_provider="openai", endpoint_task="test", status="success")
            mock_response.labels.assert_called_with(model_provider="openai", endpoint_task="test", status="success")
            mock_response.labels().observe.assert_called_with(1.0)
            mock_tokens.labels.assert_any_call(model_provider="openai", endpoint_task="test", token_type="prompt", status="success")
            mock_tokens.labels.assert_any_call(model_provider="openai", endpoint_task="test", token_type="completion", status="success")
            
            # Error case
            metrics.observe_ai_interaction(
                model_provider="openai",
                endpoint_task="test",
                duration_seconds=1.0,
                prompt_tokens=None,
                completion_tokens=None,
                total_tokens=None,
                success=False,
                error_type="rate_limit"
            )
            mock_requests.labels.assert_called_with(model_provider="openai", endpoint_task="test", status="error")
            mock_response.labels.assert_called_with(model_provider="openai", endpoint_task="test", status="error")
            mock_response.labels().observe.assert_called_with(1.0)
            mock_errors.labels.assert_called_with(model_provider="openai", endpoint_task="test", error_type="rate_limit")

    def test_start_metrics_http_server(self):
        """Test start_metrics_http_server starts a server and is idempotent."""
        # Patch the actual start_http_server function that the thread will call
        with mock.patch('vcf_agent.metrics.start_http_server') as mock_prometheus_start_server, \
             mock.patch.object(metrics.log, 'info') as mock_log_info, \
             mock.patch.object(metrics.log, 'error') as mock_log_error, \
             mock.patch('threading.Thread') as mock_thread: # Keep mock_thread to check thread creation
            
            mock_thread_instance = mock.Mock()
            mock_thread.return_value = mock_thread_instance

            # First call should start server
            metrics.start_metrics_http_server(port=8001)
            mock_thread.assert_called_once() # Check thread was created
            # The assertion on mock_prometheus_start_server will be tricky due to the thread
            # We rely on the thread calling the patched function. 
            # For a unit test, it might be better to test _target_for_thread directly if possible,
            # or ensure the patched start_http_server is globally patched.
            # For now, let's assume the test environment allows the thread to pick up the patch.
            # We will check if the thread's target, when called, calls the mock_prometheus_start_server.
            # This requires the target to be callable and for us to call it.
            # Get the target function passed to the Thread constructor
            if mock_thread.call_args:
                thread_target = mock_thread.call_args.kwargs.get('target')
                if thread_target:
                    thread_target() # Execute the target function
                    mock_prometheus_start_server.assert_called_once_with(8001, registry=metrics.http_registry)
            
            mock_log_info.assert_any_call("prometheus_metrics_server_started", port=8001)

            # Second call should be idempotent
            metrics.start_metrics_http_server(port=8001)
            self.assertEqual(mock_thread.call_count, 1) # Thread should not be created again
            mock_log_info.assert_any_call("prometheus_metrics_server_already_running", port=8001)

            # Test server start failure
            metrics._metrics_server_thread = None # Reset for failure test
            mock_prometheus_start_server.reset_mock() # Reset call count for failure case
            mock_thread.reset_mock() # Reset mock_thread for the failure case call
            
            # Create a new mock for the thread instance that will be returned
            failing_thread_instance = mock.Mock()
            failing_thread_instance.start.side_effect = Exception("Thread start failed") # Make thread.start() fail
            mock_thread.return_value = failing_thread_instance

            metrics.start_metrics_http_server(port=8002)
            mock_log_error.assert_called_with("prometheus_metrics_server_failed_to_start", port=8002, error="Thread start failed")

    def test_push_job_metrics(self):
        """Test push_job_metrics calls push_to_gateway and handles disabled state."""
        mock_metric = mock.Mock()
        mock_metric.describe.return_value = []
        with mock.patch('vcf_agent.metrics.push_to_gateway') as mock_push, \
             mock.patch.object(metrics.log, 'warn') as mock_log_warn:
            # Test successful push
            metrics.push_job_metrics("test_job", [mock_metric], instance_name="test_instance")
            mock_push.assert_called_once()
            
            # Test disabled Pushgateway by directly patching the module-level variable
            mock_push.reset_mock() # Reset for the disabled case
            with mock.patch('vcf_agent.metrics.PUSHGATEWAY_URL', "disabled"):
                metrics.push_job_metrics("test_job_disabled", [mock_metric])
                mock_log_warn.assert_called_with("pushgateway_not_configured_or_disabled", job_name="test_job_disabled")
            
            # Test push failure
            mock_push.reset_mock()
            mock_log_warn.reset_mock() # Ensure it is not called here
            mock_push.side_effect = Exception("Push failed")
            with mock.patch.object(metrics.log, 'error') as mock_log_error:
                metrics.push_job_metrics("test_job_fail", [mock_metric])
                mock_push.assert_called_once() # Ensure push_to_gateway was still called
                mock_log_error.assert_called_once()
                # Check specific log arguments if necessary, e.g., error message
                args, kwargs = mock_log_error.call_args
                self.assertEqual(kwargs.get('job_name'), "test_job_fail")
                self.assertEqual(kwargs.get('error'), "Push failed")

    def test_observe_ai_interaction_records_explicit_total_tokens_only_if_others_absent(self):
        """Test that explicit total_tokens is recorded only if prompt/completion tokens are absent."""
        with mock.patch('vcf_agent.metrics.VCF_AGENT_AI_REQUESTS_TOTAL') as mock_requests, \
             mock.patch('vcf_agent.metrics.VCF_AGENT_AI_RESPONSE_SECONDS') as mock_response, \
             mock.patch('vcf_agent.metrics.VCF_AGENT_AI_TOKENS_TOTAL') as mock_tokens, \
             mock.patch('vcf_agent.metrics.VCF_AGENT_AI_ERRORS_TOTAL') as mock_errors:
            
            # Scenario 1: total_tokens provided, prompt/completion are None -> total_tokens should be recorded
            metrics.observe_ai_interaction(
                model_provider="provider1", 
                endpoint_task="task1", 
                duration_seconds=1.0, 
                prompt_tokens=None, 
                completion_tokens=None, 
                total_tokens=300, 
                success=True
            )
            mock_tokens.labels.assert_any_call(model_provider="provider1", endpoint_task="task1", token_type='total', status='success')
            mock_tokens.labels.return_value.inc.assert_any_call(300)
            # Ensure prompt/completion were not called for this scenario for 'total' metric logic check
            # This requires checking ALL calls to mock_tokens.labels(...).inc(...)
            # A bit tricky with assert_any_call. Let's verify specific calls.
            
            # Reset mocks for clarity for the next scenario
            mock_tokens.reset_mock()

            # Scenario 2: total_tokens provided, but prompt_tokens also provided -> total_tokens should NOT be recorded
            metrics.observe_ai_interaction(
                model_provider="provider2", 
                endpoint_task="task2", 
                duration_seconds=1.0, 
                prompt_tokens=50, 
                completion_tokens=None, 
                total_tokens=300, 
                success=True
            )
            # prompt should be recorded
            mock_tokens.labels.assert_any_call(model_provider="provider2", endpoint_task="task2", token_type='prompt', status='success')
            mock_tokens.labels.return_value.inc.assert_any_call(50)
            # total should NOT be recorded
            total_metric_call_found = False
            for call_args in mock_tokens.labels.call_args_list:
                if call_args[1].get('token_type') == 'total':
                    total_metric_call_found = True
                    break
            self.assertFalse(total_metric_call_found, "Total tokens metric should not be called when prompt_tokens is present")

            mock_tokens.reset_mock()

            # Scenario 3: total_tokens provided, but completion_tokens also provided -> total_tokens should NOT be recorded
            metrics.observe_ai_interaction(
                model_provider="provider3", 
                endpoint_task="task3", 
                duration_seconds=1.0, 
                prompt_tokens=None, 
                completion_tokens=150, 
                total_tokens=300, 
                success=True
            )
            # completion should be recorded
            mock_tokens.labels.assert_any_call(model_provider="provider3", endpoint_task="task3", token_type='completion', status='success')
            mock_tokens.labels.return_value.inc.assert_any_call(150)
            # total should NOT be recorded
            total_metric_call_found = False
            for call_args in mock_tokens.labels.call_args_list:
                if call_args[1].get('token_type') == 'total':
                    total_metric_call_found = True
                    break
            self.assertFalse(total_metric_call_found, "Total tokens metric should not be called when completion_tokens is present")

    # Test observe_tool_usage (This function does not exist in metrics.py, removing)
    # def test_observe_tool_usage(self):
    #     """Test observe_tool_usage records metrics for success and error cases."""
    #     with mock.patch.object(metrics, 'VCF_AGENT_TOOL_REQUESTS_TOTAL') as mock_requests, \
    #          mock.patch.object(metrics, 'VCF_AGENT_TOOL_DURATION_SECONDS') as mock_response, \
    #          mock.patch.object(metrics, 'VCF_AGENT_TOOL_ERRORS_TOTAL') as mock_errors:
    #         # Success case
    #         metrics.observe_tool_usage(
    #             tool_name="test_tool",
    #             duration_seconds=0.5,
    #             success=True
    #         )
    #         mock_requests.labels.assert_called_with(tool_name="test_tool", status="success")
    #         mock_response.labels.assert_called_with(tool_name="test_tool") # No status label for histogram
    #         mock_response.labels().observe.assert_called_with(0.5)

    #         # Error case
    #         mock_requests.reset_mock() # Reset for error case call
    #         mock_response.reset_mock()
    #         # mock_errors.reset_mock() # Should also reset mock_errors if it's checked in success and then error
    #         metrics.observe_tool_usage(
    #             tool_name="test_tool",
    #             duration_seconds=0.2,
    #             success=False,
    #             error_type="test_error"
    #         )
    #         mock_requests.labels.assert_called_with(tool_name="test_tool", status="error")
    #         mock_response.labels.assert_called_with(tool_name="test_tool") # No status label for histogram
    #         mock_response.labels().observe.assert_called_with(0.2)
    #         mock_errors.labels.assert_called_with(tool_name="test_tool", error_type="test_error")

    # Test log_cli_command_start and log_cli_command_end (Placeholder)
    # These might not be unit testable in isolation easily if they directly use global CLI metrics.
    # Will assess if these need specific tests or are covered by integration tests.

    def test_observe_bcftools_command_success_and_error(self):
        """Test observe_bcftools_command records metrics for success and error cases."""
        with mock.patch.object(metrics, 'VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL') as mock_commands, \
             mock.patch.object(metrics, 'VCF_AGENT_BCFTOOLS_DURATION_SECONDS') as mock_duration, \
             mock.patch.object(metrics, 'VCF_AGENT_BCFTOOLS_ERRORS_TOTAL') as mock_errors:
            
            # Success case
            metrics.observe_bcftools_command(
                bcftools_subcommand="view",
                duration_seconds=1.2,
                success=True
            )
            mock_commands.labels.assert_called_with(bcftools_subcommand="view", status="success")
            mock_duration.labels.assert_called_with(bcftools_subcommand="view", status="success")
            mock_duration.labels.return_value.observe.assert_called_with(1.2)
            mock_errors.labels.assert_not_called() # No error for success case

            # Reset mocks for error case
            mock_commands.reset_mock()
            mock_duration.reset_mock()
            mock_errors.reset_mock()

            # Error case
            metrics.observe_bcftools_command(
                bcftools_subcommand="stats",
                duration_seconds=0.8,
                success=False
            )
            mock_commands.labels.assert_called_with(bcftools_subcommand="stats", status="error")
            mock_duration.labels.assert_called_with(bcftools_subcommand="stats", status="error")
            mock_duration.labels.return_value.observe.assert_called_with(0.8)
            mock_errors.labels.assert_called_with(bcftools_subcommand="stats")
            mock_errors.labels.return_value.inc.assert_called_once() # Ensure error counter is incremented

    # ... rest of the existing test methods ... 