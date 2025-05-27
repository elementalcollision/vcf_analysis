"""
Unit tests for GATK integration module.

Tests the Python wrapper for GATK ValidateVariants command,
including subprocess execution, error handling, and edge cases.
"""

import pytest
import subprocess
from unittest.mock import patch, MagicMock
from vcf_agent.gatk_integration import run_gatk_command, gatk_validatevariants


class TestRunGatkCommand:
    """Test cases for the run_gatk_command function."""

    @patch('vcf_agent.gatk_integration.subprocess.run')
    def test_run_gatk_command_success(self, mock_run):
        """Test successful GATK command execution."""
        # Setup mock
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = b"GATK output"
        mock_result.stderr = b""
        mock_run.return_value = mock_result

        # Execute
        command = ["ValidateVariants", "-V", "test.vcf"]
        rc, stdout, stderr = run_gatk_command(command)

        # Verify
        assert rc == 0
        assert stdout == "GATK output"
        assert stderr == ""
        mock_run.assert_called_once_with(
            ["gatk", "ValidateVariants", "-V", "test.vcf"],
            input=None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )

    @patch('vcf_agent.gatk_integration.subprocess.run')
    def test_run_gatk_command_with_input_data(self, mock_run):
        """Test GATK command execution with input data."""
        # Setup mock
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = b"GATK output with input"
        mock_result.stderr = b""
        mock_run.return_value = mock_result

        # Execute
        command = ["ValidateVariants", "-V", "test.vcf"]
        input_data = b"some input data"
        rc, stdout, stderr = run_gatk_command(command, input_data)

        # Verify
        assert rc == 0
        assert stdout == "GATK output with input"
        assert stderr == ""
        mock_run.assert_called_once_with(
            ["gatk", "ValidateVariants", "-V", "test.vcf"],
            input=input_data,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )

    @patch('vcf_agent.gatk_integration.subprocess.run')
    def test_run_gatk_command_failure(self, mock_run):
        """Test GATK command execution with non-zero return code."""
        # Setup mock
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stdout = b""
        mock_result.stderr = b"GATK error message"
        mock_run.return_value = mock_result

        # Execute
        command = ["ValidateVariants", "-V", "invalid.vcf"]
        rc, stdout, stderr = run_gatk_command(command)

        # Verify
        assert rc == 1
        assert stdout == ""
        assert stderr == "GATK error message"

    @patch('vcf_agent.gatk_integration.subprocess.run')
    def test_run_gatk_command_gatk_not_found(self, mock_run):
        """Test GATK command execution when GATK is not installed."""
        # Setup mock to raise FileNotFoundError
        mock_run.side_effect = FileNotFoundError("No such file or directory: 'gatk'")

        # Execute and verify exception
        command = ["ValidateVariants", "-V", "test.vcf"]
        with pytest.raises(FileNotFoundError, match="GATK is not installed or not in PATH"):
            run_gatk_command(command)

    @patch('vcf_agent.gatk_integration.subprocess.run')
    def test_run_gatk_command_empty_command(self, mock_run):
        """Test GATK command execution with empty command list."""
        # Setup mock
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = b""
        mock_result.stderr = b""
        mock_run.return_value = mock_result

        # Execute
        command = []
        rc, stdout, stderr = run_gatk_command(command)

        # Verify
        assert rc == 0
        assert stdout == ""
        assert stderr == ""
        mock_run.assert_called_once_with(
            ["gatk"],
            input=None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )

    @patch('vcf_agent.gatk_integration.subprocess.run')
    def test_run_gatk_command_unicode_output(self, mock_run):
        """Test GATK command execution with unicode characters in output."""
        # Setup mock with unicode characters
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "GATK output with unicode: ñáéíóú".encode('utf-8')
        mock_result.stderr = "Warning with unicode: ñáéíóú".encode('utf-8')
        mock_run.return_value = mock_result

        # Execute
        command = ["ValidateVariants", "-V", "test.vcf"]
        rc, stdout, stderr = run_gatk_command(command)

        # Verify
        assert rc == 0
        assert stdout == "GATK output with unicode: ñáéíóú"
        assert stderr == "Warning with unicode: ñáéíóú"


class TestGatkValidateVariants:
    """Test cases for the gatk_validatevariants function."""

    @patch('vcf_agent.gatk_integration.run_gatk_command')
    def test_gatk_validatevariants_success(self, mock_run_gatk):
        """Test successful GATK ValidateVariants execution."""
        # Setup mock
        mock_run_gatk.return_value = (0, "Validation successful", "")

        # Execute
        args = ["-V", "test.vcf", "-R", "reference.fasta"]
        rc, stdout, stderr = gatk_validatevariants(args)

        # Verify
        assert rc == 0
        assert stdout == "Validation successful"
        assert stderr == ""
        mock_run_gatk.assert_called_once_with(["ValidateVariants", "-V", "test.vcf", "-R", "reference.fasta"], None)

    @patch('vcf_agent.gatk_integration.run_gatk_command')
    def test_gatk_validatevariants_with_input_data(self, mock_run_gatk):
        """Test GATK ValidateVariants execution with input data."""
        # Setup mock
        mock_run_gatk.return_value = (0, "Validation with input successful", "")

        # Execute
        args = ["-V", "test.vcf", "-R", "reference.fasta"]
        input_data = b"input data for validation"
        rc, stdout, stderr = gatk_validatevariants(args, input_data)

        # Verify
        assert rc == 0
        assert stdout == "Validation with input successful"
        assert stderr == ""
        mock_run_gatk.assert_called_once_with(
            ["ValidateVariants", "-V", "test.vcf", "-R", "reference.fasta"], 
            input_data
        )

    @patch('vcf_agent.gatk_integration.run_gatk_command')
    def test_gatk_validatevariants_failure(self, mock_run_gatk):
        """Test GATK ValidateVariants execution with validation errors."""
        # Setup mock
        mock_run_gatk.return_value = (1, "", "Validation failed: Invalid VCF format")

        # Execute
        args = ["-V", "invalid.vcf", "-R", "reference.fasta"]
        rc, stdout, stderr = gatk_validatevariants(args)

        # Verify
        assert rc == 1
        assert stdout == ""
        assert stderr == "Validation failed: Invalid VCF format"
        mock_run_gatk.assert_called_once_with(["ValidateVariants", "-V", "invalid.vcf", "-R", "reference.fasta"], None)

    @patch('vcf_agent.gatk_integration.run_gatk_command')
    def test_gatk_validatevariants_empty_args(self, mock_run_gatk):
        """Test GATK ValidateVariants execution with empty arguments."""
        # Setup mock
        mock_run_gatk.return_value = (1, "", "Missing required arguments")

        # Execute
        args = []
        rc, stdout, stderr = gatk_validatevariants(args)

        # Verify
        assert rc == 1
        assert stdout == ""
        assert stderr == "Missing required arguments"
        mock_run_gatk.assert_called_once_with(["ValidateVariants"], None)

    @patch('vcf_agent.gatk_integration.run_gatk_command')
    def test_gatk_validatevariants_complex_args(self, mock_run_gatk):
        """Test GATK ValidateVariants execution with complex arguments."""
        # Setup mock
        mock_run_gatk.return_value = (0, "Complex validation successful", "")

        # Execute
        args = [
            "-V", "test.vcf",
            "-R", "reference.fasta",
            "--validate-GVCF",
            "--validation-type-to-exclude", "ALL",
            "--dbsnp", "dbsnp.vcf"
        ]
        rc, stdout, stderr = gatk_validatevariants(args)

        # Verify
        assert rc == 0
        assert stdout == "Complex validation successful"
        assert stderr == ""
        expected_command = [
            "ValidateVariants",
            "-V", "test.vcf",
            "-R", "reference.fasta",
            "--validate-GVCF",
            "--validation-type-to-exclude", "ALL",
            "--dbsnp", "dbsnp.vcf"
        ]
        mock_run_gatk.assert_called_once_with(expected_command, None)

    @patch('vcf_agent.gatk_integration.run_gatk_command')
    def test_gatk_validatevariants_propagates_exception(self, mock_run_gatk):
        """Test that gatk_validatevariants propagates exceptions from run_gatk_command."""
        # Setup mock to raise exception
        mock_run_gatk.side_effect = FileNotFoundError("GATK is not installed or not in PATH")

        # Execute and verify exception propagation
        args = ["-V", "test.vcf", "-R", "reference.fasta"]
        with pytest.raises(FileNotFoundError, match="GATK is not installed or not in PATH"):
            gatk_validatevariants(args)


class TestGatkIntegrationEdgeCases:
    """Test edge cases and error conditions for GATK integration."""

    @patch('vcf_agent.gatk_integration.subprocess.run')
    def test_run_gatk_command_large_output(self, mock_run):
        """Test GATK command execution with large output."""
        # Setup mock with large output
        large_output = "A" * 10000  # 10KB of output
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = large_output.encode('utf-8')
        mock_result.stderr = b""
        mock_run.return_value = mock_result

        # Execute
        command = ["ValidateVariants", "-V", "large.vcf"]
        rc, stdout, stderr = run_gatk_command(command)

        # Verify
        assert rc == 0
        assert stdout == large_output
        assert stderr == ""

    @patch('vcf_agent.gatk_integration.subprocess.run')
    def test_run_gatk_command_special_characters(self, mock_run):
        """Test GATK command execution with special characters in arguments."""
        # Setup mock
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = b"Command with special chars executed"
        mock_result.stderr = b""
        mock_run.return_value = mock_result

        # Execute with special characters
        command = ["ValidateVariants", "-V", "file with spaces.vcf", "--arg", "value=with=equals"]
        rc, stdout, stderr = run_gatk_command(command)

        # Verify
        assert rc == 0
        assert stdout == "Command with special chars executed"
        assert stderr == ""
        mock_run.assert_called_once_with(
            ["gatk", "ValidateVariants", "-V", "file with spaces.vcf", "--arg", "value=with=equals"],
            input=None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )

    @patch('vcf_agent.gatk_integration.subprocess.run')
    def test_run_gatk_command_binary_input(self, mock_run):
        """Test GATK command execution with binary input data."""
        # Setup mock
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = b"Binary input processed"
        mock_result.stderr = b""
        mock_run.return_value = mock_result

        # Execute with binary input
        command = ["ValidateVariants", "-V", "test.vcf"]
        binary_input = bytes([0, 1, 2, 3, 255, 254, 253])
        rc, stdout, stderr = run_gatk_command(command, binary_input)

        # Verify
        assert rc == 0
        assert stdout == "Binary input processed"
        assert stderr == ""
        mock_run.assert_called_once_with(
            ["gatk", "ValidateVariants", "-V", "test.vcf"],
            input=binary_input,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        ) 