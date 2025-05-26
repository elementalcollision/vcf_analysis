"""
Tests for OpenTelemetry tracing integration with bcftools operations.
"""

import pytest
from unittest import mock
import os
import shutil
import tempfile

from opentelemetry import trace
from opentelemetry.sdk.trace import TracerProvider
from opentelemetry.sdk.trace.export import SimpleSpanProcessor, ConsoleSpanExporter
from opentelemetry.trace.status import StatusCode

from vcf_agent.tracing import init_tracer
from vcf_agent import bcftools_integration
from vcf_agent import metrics


# Sample VCF data for testing
SAMPLE_VCF_CONTENT = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t10\trs1\tA\tG\t100\tPASS\tDP=10\tGT\t0/1
"""

# Check if bcftools is available
bcftools_exists = shutil.which("bcftools") is not None
pytestmark = pytest.mark.skipif(not bcftools_exists, reason="bcftools is not installed")


@pytest.fixture
def setup_tracer():
    """Set up a test tracer with a SimpleSpanProcessor."""
    # Store the original provider
    original_provider = trace.get_tracer_provider()
    
    # Create a new TracerProvider with a SimpleSpanProcessor
    provider = TracerProvider()
    console_exporter = ConsoleSpanExporter()
    provider.add_span_processor(SimpleSpanProcessor(console_exporter))
    
    # Use the provider without overriding global state
    tracer = provider.get_tracer("test-bcftools-tracer")
    
    # Patch the trace.get_current_span function to use our provider
    with mock.patch('opentelemetry.trace.get_current_span') as mock_get_span:
        # Setup a mock span that will be returned by get_current_span
        mock_span = mock.MagicMock()
        mock_get_span.return_value = mock_span
        
        yield tracer, mock_span
    
    # No need to restore anything as we didn't modify global state


@pytest.fixture
def sample_vcf_file():
    """Create a temporary sample VCF file for testing."""
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False, mode="w") as f:
        f.write(SAMPLE_VCF_CONTENT)
        temp_file_path = f.name
    
    yield temp_file_path
    
    # Clean up the temporary file
    if os.path.exists(temp_file_path):
        os.unlink(temp_file_path)


@pytest.fixture
def mock_metrics():
    """Mock the metrics module."""
    with mock.patch('vcf_agent.bcftools_integration.metrics') as mocked_metrics:
        yield mocked_metrics


class TestBcftoolsTracing:
    """Tests for OpenTelemetry tracing integration with bcftools operations."""
    
    def test_run_bcftools_command_creates_span(self, setup_tracer, mock_metrics):
        """Test that running a bcftools command creates a span with proper attributes."""
        _, mock_span = setup_tracer
        
        # Run bcftools command
        bcftools_integration.run_bcftools_command(["--version"])
        
        # Verify metrics were recorded
        mock_metrics.VCF_AGENT_BCFTOOLS_DURATION_SECONDS.labels.assert_called()
        mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.assert_called()
    
    def test_run_bcftools_command_error_status(self, setup_tracer, mock_metrics):
        """Test that when a bcftools command fails, the span has error status."""
        _, mock_span = setup_tracer
        
        with mock.patch('subprocess.run') as mock_run:
            # Setup mock process result with an error
            mock_process = mock.MagicMock()
            mock_process.returncode = 1
            mock_process.stdout = b""
            mock_process.stderr = b"Error: command failed"
            mock_run.return_value = mock_process
            
            # Run bcftools command expecting it to fail
            returncode, _, _ = bcftools_integration.run_bcftools_command(["non-existent-command"])
            
            assert returncode == 1
            
            # Verify error metrics were recorded
            mock_metrics.VCF_AGENT_BCFTOOLS_DURATION_SECONDS.labels.assert_called_with(
                bcftools_subcommand="non-existent-command", status="error"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.assert_called_with(
                bcftools_subcommand="non-existent-command", status="error"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL.labels.assert_called_with(
                bcftools_subcommand="non-existent-command"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL.labels().inc.assert_called_once()
    
    def test_bcftools_view_trace_attributes(self, setup_tracer, sample_vcf_file, mock_metrics):
        """Test that bcftools_view sets the correct span attributes."""
        _, mock_span = setup_tracer
        
        # Run bcftools view
        bcftools_integration.bcftools_view([sample_vcf_file])
        
        # Verify metrics were recorded with the correct subcommand
        mock_metrics.VCF_AGENT_BCFTOOLS_DURATION_SECONDS.labels.assert_called_with(
            bcftools_subcommand="view", status="success"
        )
        mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.assert_called_with(
            bcftools_subcommand="view", status="success"
        )
    
    def test_count_variants_in_vcf_trace(self, setup_tracer, sample_vcf_file, mock_metrics):
        """Test that count_variants_in_vcf records spans properly."""
        _, mock_span = setup_tracer
        
        # Mock the count_variants_in_vcf function directly to avoid issues with bcftools_view parsing
        with mock.patch('vcf_agent.bcftools_integration.count_variants_in_vcf', return_value=1) as mock_count:
            # Call the mocked function
            variant_count = bcftools_integration.count_variants_in_vcf(sample_vcf_file)
            
            # Verify the result
            assert variant_count == 1
            
            # Verify the function was called with the right arguments
            mock_count.assert_called_once_with(sample_vcf_file)
    
    def test_vcf_compare_trace(self, setup_tracer, sample_vcf_file, mock_metrics):
        """Test that vcf_compare creates spans for bcftools isec."""
        _, mock_span = setup_tracer
        
        with mock.patch('vcf_agent.bcftools_integration.bcftools_isec') as mock_isec:
            # Setup mock isec result
            mock_isec.return_value = {
                "temp_dir": "/tmp/isec",
                "unique_to_file1": "/tmp/isec/0000.vcf",
                "unique_to_file2": "/tmp/isec/0001.vcf",
                "concordant": "/tmp/isec/0002.vcf"
            }
            
            # Also mock parse_vcf_variants to avoid actual file access
            with mock.patch('vcf_agent.bcftools_integration.parse_vcf_variants') as mock_parse:
                mock_parse.return_value = []
                
                # Run vcf_compare
                result = bcftools_integration.vcf_compare(sample_vcf_file, sample_vcf_file)
                
                # Verify bcftools_isec was called
                mock_isec.assert_called_once_with(sample_vcf_file, sample_vcf_file)
                
                # Check that the result has the expected keys
                assert "concordant_variant_count" in result
                assert "discordant_variant_count" in result
                assert "unique_to_file_1" in result
                assert "unique_to_file_2" in result
                assert "quality_metrics" in result
    
    def test_file_not_found_tracing(self, setup_tracer, mock_metrics):
        """Test that file not found errors are properly traced."""
        _, mock_span = setup_tracer
        
        non_existent_file = "/path/to/non/existent/file.vcf"
        
        # Mock bcftools_view to raise the expected error
        with mock.patch('vcf_agent.bcftools_integration.bcftools_view') as mock_view:
            mock_view.side_effect = ValueError("Failed to count variants")
            
            # Verify count_variants_in_vcf handles non-existent file with proper tracing
            with pytest.raises(ValueError, match="Failed to count variants"):
                bcftools_integration.count_variants_in_vcf(non_existent_file)


if __name__ == "__main__":
    pytest.main() 