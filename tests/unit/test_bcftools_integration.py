import pytest
from unittest import mock
import subprocess
import tempfile
import os
import shutil
import gzip

from vcf_agent import bcftools_integration
from vcf_agent import metrics # To allow mocking metrics

# Assume bcftools is installed and in PATH for most tests.
# Tests for bcftools not found will mock subprocess.run.

@pytest.fixture
def mock_metrics():
    with mock.patch('vcf_agent.bcftools_integration.metrics') as mocked_metrics:
        # Configure individual metric objects if needed, e.g.:
        # mocked_metrics.VCF_AGENT_BCFTOOLS_DURATION_SECONDS = mock.Mock()
        # mocked_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL = mock.Mock()
        # mocked_metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL = mock.Mock()
        # mocked_metrics.log = mock.Mock()
        yield mocked_metrics


class TestRunBcftoolsCommand:
    def test_run_successful_command(self, mock_metrics):
        # Example: bcftools --version (should always succeed)
        # We need a real bcftools command that produces predictable output and exit code 0
        # For now, let's assume 'bcftools --version' is suitable, though it might vary slightly.
        # A better approach would be a simple bcftools view on a tiny, known VCF header.
        
        # Using a more stable command like 'echo' for initial setup if bcftools isn't guaranteed
        # For this test, we'll mock subprocess.run to simulate bcftools
        mock_process = mock.Mock()
        mock_process.returncode = 0
        mock_process.stdout = b"bcftools 1.x.y\n"
        mock_process.stderr = b""

        with mock.patch('subprocess.run', return_value=mock_process) as mock_subproc_run:
            return_code, stdout, stderr = bcftools_integration.run_bcftools_command(["--version"])
            
            assert return_code == 0
            assert "bcftools" in stdout
            assert stderr == ""
            mock_subproc_run.assert_called_once_with(
                ["bcftools", "--version"],
                input=None,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False
            )
            # Assert metrics were called
            mock_metrics.VCF_AGENT_BCFTOOLS_DURATION_SECONDS.labels.assert_called_once_with(
                bcftools_subcommand="--version", status="success"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.assert_called_once_with(
                bcftools_subcommand="--version", status="success"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL.labels().inc.assert_not_called()

    def test_run_command_with_error(self, mock_metrics):
        mock_process = mock.Mock()
        mock_process.returncode = 1
        mock_process.stdout = b""
        mock_process.stderr = b"Error: some bcftools error\n"

        with mock.patch('subprocess.run', return_value=mock_process) as mock_subproc_run:
            return_code, stdout, stderr = bcftools_integration.run_bcftools_command(["view", "non_existent_file.vcf"])
            
            assert return_code == 1
            assert stdout == ""
            assert "some bcftools error" in stderr
            mock_subproc_run.assert_called_once_with(
                ["bcftools", "view", "non_existent_file.vcf"],
                input=None,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_DURATION_SECONDS.labels.assert_called_once_with(
                bcftools_subcommand="view", status="error"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.assert_called_once_with(
                bcftools_subcommand="view", status="error"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL.labels.assert_called_once_with(
                bcftools_subcommand="view"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL.labels().inc.assert_called_once()

    def test_run_command_bcftools_not_found(self, mock_metrics):
        with mock.patch('subprocess.run', side_effect=FileNotFoundError("bcftools not found")) as mock_subproc_run:
            with pytest.raises(FileNotFoundError, match="bcftools is not installed or not in PATH"):
                bcftools_integration.run_bcftools_command(["--version"])
            
            mock_subproc_run.assert_called_once()
            # Check that specific log for bcftools_not_found was made
            mock_metrics.log.error.assert_any_call("bcftools_not_found", command_used=['bcftools', '--version'], error="bcftools not found")
            
            # General metrics for error case
            mock_metrics.VCF_AGENT_BCFTOOLS_DURATION_SECONDS.labels.assert_called_once_with(
                bcftools_subcommand="--version", status="error"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.assert_called_once_with(
                bcftools_subcommand="--version", status="error"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL.labels.assert_called_once_with(
                bcftools_subcommand="--version"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL.labels().inc.assert_called_once()

    def test_run_command_with_input_data(self, mock_metrics):
        mock_process = mock.Mock()
        mock_process.returncode = 0
        mock_process.stdout = b"Processed output"
        mock_process.stderr = b""
        input_content = b"vcf data here"

        with mock.patch('subprocess.run', return_value=mock_process) as mock_subproc_run:
            return_code, stdout, stderr = bcftools_integration.run_bcftools_command(["view", "-"], input_data=input_content)
            
            assert return_code == 0
            assert stdout == "Processed output"
            mock_subproc_run.assert_called_once_with(
                ["bcftools", "view", "-"],
                input=input_content,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.assert_called_once_with(
                bcftools_subcommand="view", status="success"
            )

    def test_run_command_unexpected_subprocess_error(self, mock_metrics):
        with mock.patch('subprocess.run', side_effect=RuntimeError("Unexpected subprocess issue")) as mock_subproc_run:
            with pytest.raises(RuntimeError, match="Unexpected subprocess issue"):
                bcftools_integration.run_bcftools_command(["--version"])
            
            mock_subproc_run.assert_called_once()
            mock_metrics.log.error.assert_any_call("bcftools_run_unexpected_error", command_used=['bcftools', '--version'], error="Unexpected subprocess issue")
            mock_metrics.VCF_AGENT_BCFTOOLS_DURATION_SECONDS.labels.assert_called_once_with(
                bcftools_subcommand="--version", status="error"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.assert_called_once_with(
                bcftools_subcommand="--version", status="error"
            )
            mock_metrics.VCF_AGENT_BCFTOOLS_ERRORS_TOTAL.labels.assert_called_once_with(
                bcftools_subcommand="--version"
            )

# Further tests would cover the wrapper functions (bcftools_view, etc.)
# primarily by mocking run_bcftools_command and asserting it's called correctly.

# Tests for bcftools_isec, count_variants_in_vcf, etc. would require actual
# sample VCF files and more complex setup/assertions. 

class TestBcftoolsWrapperFunctions:
    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_bcftools_view_wrapper(self, mock_run_command):
        mock_run_command.return_value = (0, "view_stdout", "view_stderr")
        args = ["-H", "somefile.vcf.gz"]
        input_d = b"input_data"
        
        ret_code, stdout, stderr = bcftools_integration.bcftools_view(args, input_data=input_d)
        
        mock_run_command.assert_called_once_with(["view"] + args, input_d)
        assert ret_code == 0
        assert stdout == "view_stdout"
        assert stderr == "view_stderr"

    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_bcftools_query_wrapper(self, mock_run_command):
        mock_run_command.return_value = (0, "query_stdout", "")
        args = ["-f", "%CHROM", "somefile.vcf.gz"]
        
        ret_code, stdout, stderr = bcftools_integration.bcftools_query(args)
        
        mock_run_command.assert_called_once_with(["query"] + args, None)
        assert ret_code == 0
        assert stdout == "query_stdout"

    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_bcftools_filter_wrapper(self, mock_run_command):
        args = ["-e", "QUAL>10", "somefile.vcf.gz"]
        bcftools_integration.bcftools_filter(args)
        mock_run_command.assert_called_once_with(["filter"] + args, None)

    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_bcftools_norm_wrapper(self, mock_run_command):
        args = ["-m", "-any", "somefile.vcf.gz"]
        bcftools_integration.bcftools_norm(args)
        mock_run_command.assert_called_once_with(["norm"] + args, None)

    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_bcftools_stats_wrapper(self, mock_run_command):
        args = ["somefile.vcf.gz"]
        bcftools_integration.bcftools_stats(args)
        mock_run_command.assert_called_once_with(["stats"] + args, None)

    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_bcftools_annotate_wrapper(self, mock_run_command):
        args = ["-a", "ann.vcf", "somefile.vcf.gz"]
        bcftools_integration.bcftools_annotate(args)
        mock_run_command.assert_called_once_with(["annotate"] + args, None)

# Further tests would cover the wrapper functions (bcftools_view, etc.) 

TEST_DATA_DIR = "tests/unit/test_data"
SAMPLE1_VCF = f"{TEST_DATA_DIR}/sample1.vcf"
SAMPLE2_VCF = f"{TEST_DATA_DIR}/sample2.vcf"
DETAILED_VCF = f"{TEST_DATA_DIR}/detailed.vcf"

# Helper function to bgzip a file
def bgzip_file(input_path, output_path):
    with open(input_path, "rb") as f_in, gzip.open(output_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return output_path

class TestBcftoolsFileOperations:
    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_bcftools_isec_successful(self, mock_run_bcftools_command_for_isec, mock_metrics):
        import tempfile
        import os
        import shutil

        # This mock_run_bcftools_command_for_isec is specifically for the isec test.
        # It needs to allow count_variants_in_vcf to use its own mock or real bcftools.
        # The fallback in mock_isec_effect was trying to call real bcftools, which might be the issue.
        # For this unit test, let's assume count_variants_in_vcf is tested independently
        # or also mocked if we want to purely unit test bcftools_isec structure.

        # For now, let's simplify mock_isec_effect to only handle the isec call.
        # and mock the subsequent count_variants_in_vcf calls directly if needed.
        temp_output_dir = tempfile.mkdtemp()
        def mock_isec_only_effect(command_args, input_data=None):
            if command_args[0] == "isec":
                # command_args will be ["isec", file1, file2, "-p", temp_output_dir_passed_to_isec]
                # We need to make sure the output_dir used here matches what bcftools_isec will use.
                # The actual output dir for isec files is command_args[4] if -p is command_args[3]
                isec_output_dir_from_call = command_args[4] 
                os.makedirs(isec_output_dir_from_call, exist_ok=True)
                vcf_header = (
                    "##fileformat=VCFv4.2\n"
                    "##contig=<ID=chr1>\n"
                    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                )
                # These are the files bcftools_isec expects based on its return dict
                with open(os.path.join(isec_output_dir_from_call, "0000.vcf"), "w") as f: # unique to file1
                    f.write(f"{vcf_header}chr1\t30\trs3\tG\tA\t100\tPASS\tDP=30\n") # 1 variant
                with open(os.path.join(isec_output_dir_from_call, "0001.vcf"), "w") as f: # unique to file2
                    f.write(f"{vcf_header}chr1\t40\trs4\tC\tG\t100\tPASS\tDP=40\n") # 1 variant
                with open(os.path.join(isec_output_dir_from_call, "0002.vcf"), "w") as f: # concordant
                    f.write(f"{vcf_header}chr1\t10\trs1\tA\tG\t100\tPASS\tDP=10\n") # 1 variant
                return (0, "mock isec stdout", "mock isec stderr")
            # For this specific test of bcftools_isec, we don't want a fallback to real bcftools for other commands.
            # Let calls to count_variants_in_vcf be handled by a separate mock if needed.
            pytest.fail(f"run_bcftools_command called with unexpected command in mock_isec_only_effect: {command_args}")
    
        mock_run_bcftools_command_for_isec.side_effect = mock_isec_only_effect
    
        isec_result_dir_from_func = None
        try:
            # We mock run_bcftools_command, so the actual content of SAMPLE1_VCF/SAMPLE2_VCF doesn't matter for the isec call itself.
            # The important part is that bcftools_isec passes the temp_dir to the mocked command.
            with mock.patch('vcf_agent.bcftools_integration.count_variants_in_vcf') as mock_count_variants:
                # Configure mock_count_variants to return expected counts for the files created by mock_isec_only_effect
                def count_side_effect(filepath):
                    if "0000.vcf" in filepath: return 1 # unique1
                    if "0001.vcf" in filepath: return 1 # unique2
                    if "0002.vcf" in filepath: return 1 # concordant
                    return 0
                mock_count_variants.side_effect = count_side_effect

                result = bcftools_integration.bcftools_isec(SAMPLE1_VCF, SAMPLE2_VCF)
                isec_result_dir_from_func = result.get("temp_dir") # This is the dir created by mkdtemp in bcftools_isec
        
                # Assert that the mocked run_bcftools_command was called for isec
                # The path passed to -p will be the one from result.get("temp_dir")
                mock_run_bcftools_command_for_isec.assert_called_once_with([
                    "isec", SAMPLE1_VCF, SAMPLE2_VCF, "-p", isec_result_dir_from_func
                ])

                assert "temp_dir" in result
                assert os.path.exists(result["unique_to_file1"]) # File path should exist due to mock
                assert os.path.exists(result["unique_to_file2"])
                assert os.path.exists(result["concordant"])
        
                # These calls will now use mock_count_variants
                common_count = bcftools_integration.count_variants_in_vcf(result["concordant"])
                unique1_count = bcftools_integration.count_variants_in_vcf(result["unique_to_file1"])
                unique2_count = bcftools_integration.count_variants_in_vcf(result["unique_to_file2"])
        
                assert common_count == 1
                assert unique1_count == 1 
                assert unique2_count == 1
        finally:
            if isec_result_dir_from_func and os.path.exists(isec_result_dir_from_func):
                shutil.rmtree(isec_result_dir_from_func)
            # Clean up the temp_output_dir created at the start of the test too, if it was used by the mock.
            if os.path.exists(temp_output_dir):
                 shutil.rmtree(temp_output_dir) 

    def test_bcftools_isec_file_not_found(self, mock_metrics):
        with pytest.raises(RuntimeError, match="bcftools isec failed"):
            bcftools_integration.bcftools_isec("non_existent1.vcf", "non_existent2.vcf")
        # Check that metrics reflect the error from run_bcftools_command for 'isec'
        # This test assumes run_bcftools_command correctly reports errors if bcftools fails (e.g. file not found for bcftools itself)
        # The RuntimeError is raised by bcftools_isec if run_bcftools_command returns non-zero for isec
        
        # Find the call to VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL for 'isec'
        # This is a bit tricky as metrics are global. We expect at least one 'isec' error.
        isec_error_call_found = False
        for call_args, call_kwargs in mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.call_args_list:
            if call_kwargs.get('bcftools_subcommand') == 'isec' and call_kwargs.get('status') == 'error':
                isec_error_call_found = True
                break
        assert isec_error_call_found, "Metrics for 'isec' error not found"


    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_count_variants_in_vcf(self, mock_run_bcftools_cmd, mock_metrics):
        print(f"DEBUG: Type of run_bcftools_command in test_count_variants_in_vcf: {type(bcftools_integration.run_bcftools_command)}")
        # Test with 2 variant lines
        mock_stdout_two_lines = "variant_line1\nvariant_line2"
        mock_run_bcftools_cmd.return_value = (0, mock_stdout_two_lines, "")
        assert bcftools_integration.count_variants_in_vcf("dummy_path1.vcf") == 2
        mock_run_bcftools_cmd.assert_called_with(["view", "-I", "dummy_path1.vcf"])

        # Test with 1 variant line
        mock_run_bcftools_cmd.reset_mock()
        mock_stdout_one_line = "variant_line_only"
        mock_run_bcftools_cmd.return_value = (0, mock_stdout_one_line, "")
        assert bcftools_integration.count_variants_in_vcf("dummy_path2.vcf") == 1
        mock_run_bcftools_cmd.assert_called_with(["view", "-I", "dummy_path2.vcf"])

    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_count_variants_empty_vcf(self, mock_run_bcftools_cmd, mock_metrics):
        print(f"DEBUG: Type of run_bcftools_command in test_count_variants_empty_vcf: {type(bcftools_integration.run_bcftools_command)}")
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".vcf") as tmp_file:
            tmp_file.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            tmp_file_path = tmp_file.name
        
        mock_run_bcftools_cmd.return_value = (0, "", "") # Empty stdout
        try:
            assert bcftools_integration.count_variants_in_vcf(tmp_file_path) == 0
            mock_run_bcftools_cmd.assert_called_once_with(["view", "-I", tmp_file_path])
        finally:
            os.remove(tmp_file_path)

    @mock.patch('vcf_agent.bcftools_integration.run_bcftools_command')
    def test_count_variants_no_records_vcf(self, mock_run_bcftools_cmd, mock_metrics):
        print(f"DEBUG: Type of run_bcftools_command in test_count_variants_no_records_vcf: {type(bcftools_integration.run_bcftools_command)}")
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".vcf") as tmp_file:
            tmp_file.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            tmp_file_path = tmp_file.name

        mock_run_bcftools_cmd.return_value = (0, "   \n   ", "") # Stdout with only whitespace and actual newline
        try:
            assert bcftools_integration.count_variants_in_vcf(tmp_file_path) == 0
            mock_run_bcftools_cmd.assert_called_once_with(["view", "-I", tmp_file_path])
        finally:
            os.remove(tmp_file_path)

    def test_count_variants_file_not_found(self, mock_metrics):
        # run_bcftools_command internally will try to run `bcftools view -H non_existent.vcf`
        # This should lead to an error status in metrics for the 'view' subcommand.
        with pytest.raises(ValueError, match="Failed to count variants"):
             bcftools_integration.count_variants_in_vcf("non_existent.vcf")
        
        # Check metrics for 'view' error from run_bcftools_command
        view_error_call_found = False
        # mock_metrics might be reset if fixture scope is function. Let's be careful.
        # Need to ensure this check is specific enough or the mock_metrics fixture is class/module scoped if needed.
        for call_args, call_kwargs in mock_metrics.VCF_AGENT_BCFTOOLS_COMMANDS_TOTAL.labels.call_args_list:
            if call_kwargs.get('bcftools_subcommand') == 'view' and call_kwargs.get('status') == 'error':
                # This is a weak check as other tests might also call view with error.
                # A more robust way would be to clear mocks or use a spy on run_bcftools_command.
                view_error_call_found = True
                break
        assert view_error_call_found, "Metrics for 'view' error not found in count_variants_file_not_found"


# Tests for parse_vcf_variants and vcf_compare would follow a similar pattern,
# requiring sample VCFs and asserting their output.

# Further tests would cover the wrapper functions (bcftools_view, etc.) 