"""
Unit tests for VCF utility functions (src/vcf_agent/vcf_utils.py).
"""

import pytest
from unittest.mock import patch, MagicMock, call, mock_open

# Import the module to be tested
from vcf_agent import vcf_utils
from vcf_agent import graph_integration # Used for type hinting if needed, and for mocking its methods

# Mock VCF record object that VCF reader would yield
class MockVCFRecord:
    def __init__(self, CHROM, POS, REF, ALT, gt_types, ID=None):
        self.CHROM = CHROM
        self.POS = POS
        self.REF = REF
        self.ALT = ALT # List of strings
        self.ID = ID # Optional, can be None
        self.gt_types = gt_types # Array-like, e.g., [0] for het for the first sample

class TestVCFUtilsUnit:
    """Unit tests for vcf_utils.py"""

    def test_example_placeholder_vcf_utils(self):
        """A placeholder test for vcf_utils."""
        assert True

    @patch('vcf_agent.vcf_utils.os.path.isfile', return_value=True)
    @patch('vcf_agent.vcf_utils.graph_integration')
    @patch('vcf_agent.vcf_utils.VCF')
    def test_populate_kuzu_from_vcf_basic_processing(self, mock_vcf_constructor, mock_graph_integration, mock_isfile):
        """Test basic VCF processing with one het and one hom_alt variant."""
        mock_kuzu_conn = MagicMock()
        vcf_file_path = "dummy.vcf"

        # Setup mock VCF reader and records
        mock_vcf_reader = MagicMock()
        mock_vcf_reader.samples = ["test_sample_from_vcf"]
        
        record1_het = MockVCFRecord(CHROM="chr1", POS=100, REF="A", ALT=["T"], gt_types=[0]) # 0 for HET
        record2_hom = MockVCFRecord(CHROM="chr1", POS=200, REF="G", ALT=["C"], gt_types=[1]) # 1 for HOM_ALT
        mock_vcf_reader.__iter__.return_value = iter([record1_het, record2_hom])
        mock_vcf_constructor.return_value = mock_vcf_reader

        # Call the function
        counts = vcf_utils.populate_kuzu_from_vcf(mock_kuzu_conn, vcf_file_path)

        # Assertions
        # 1. Sample addition
        mock_graph_integration.add_sample.assert_called_once_with(
            mock_kuzu_conn, 
            {"sample_id": "test_sample_from_vcf"}
        )

        # 2. Variant additions
        expected_variant_calls = [
            call(mock_kuzu_conn, {"variant_id": "chr1_100_A_T", "chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}),
            call(mock_kuzu_conn, {"variant_id": "chr1_200_G_C", "chrom": "chr1", "pos": 200, "ref": "G", "alt": "C"})
        ]
        mock_graph_integration.add_variant.assert_has_calls(expected_variant_calls)
        assert mock_graph_integration.add_variant.call_count == 2

        # 3. Link additions
        expected_link_calls = [
            call(mock_kuzu_conn, "test_sample_from_vcf", "chr1_100_A_T", {"zygosity": "HET"}),
            call(mock_kuzu_conn, "test_sample_from_vcf", "chr1_200_G_C", {"zygosity": "HOM_ALT"})
        ]
        mock_graph_integration.link_variant_to_sample.assert_has_calls(expected_link_calls)
        assert mock_graph_integration.link_variant_to_sample.call_count == 2
        
        # 4. Return counts
        assert counts == {"variants": 2, "samples": 1, "links": 2}
        
        # Ensure VCF was 'opened'
        mock_vcf_constructor.assert_called_once_with(vcf_file_path)

    @patch('vcf_agent.vcf_utils.os.path.isfile', return_value=True)
    @patch('vcf_agent.vcf_utils.graph_integration')
    @patch('vcf_agent.vcf_utils.VCF')
    def test_populate_kuzu_from_vcf_sample_override(self, mock_vcf_constructor, mock_graph_integration, mock_isfile):
        """Test VCF processing when sample_name_override is provided."""
        mock_kuzu_conn = MagicMock()
        vcf_file_path = "dummy_override.vcf"
        override_sample_name = "override_sample_001"

        mock_vcf_reader = MagicMock()
        # samples in VCF will be ignored due to override
        mock_vcf_reader.samples = ["sample_in_vcf_ignored"]
        record1 = MockVCFRecord(CHROM="chrX", POS=500, REF="C", ALT=["G"], gt_types=[0])
        mock_vcf_reader.__iter__.return_value = iter([record1])
        mock_vcf_constructor.return_value = mock_vcf_reader

        counts = vcf_utils.populate_kuzu_from_vcf(mock_kuzu_conn, vcf_file_path, sample_name_override=override_sample_name)

        mock_graph_integration.add_sample.assert_called_once_with(
            mock_kuzu_conn, 
            {"sample_id": override_sample_name}
        )
        mock_graph_integration.add_variant.assert_called_once_with(
            mock_kuzu_conn, 
            {"variant_id": "chrX_500_C_G", "chrom": "chrX", "pos": 500, "ref": "C", "alt": "G"}
        )
        mock_graph_integration.link_variant_to_sample.assert_called_once_with(
            mock_kuzu_conn, 
            override_sample_name, 
            "chrX_500_C_G", 
            {"zygosity": "HET"}
        )
        assert counts == {"variants": 1, "samples": 1, "links": 1}
        mock_vcf_constructor.assert_called_once_with(vcf_file_path)

    @patch('vcf_agent.vcf_utils.os.path.isfile', return_value=False)
    @patch('vcf_agent.vcf_utils.VCF')
    def test_populate_kuzu_from_vcf_file_not_found(self, mock_vcf_constructor, mock_isfile):
        """Test VCF processing when the VCF file is not found."""
        mock_kuzu_conn = MagicMock()
        vcf_file_path = "non_existent_file.vcf"

        with pytest.raises(FileNotFoundError) as excinfo:
            vcf_utils.populate_kuzu_from_vcf(mock_kuzu_conn, vcf_file_path)
        
        assert vcf_file_path in str(excinfo.value)
        mock_vcf_constructor.assert_not_called()
        mock_isfile.assert_called_once_with(vcf_file_path)

    @patch('vcf_agent.vcf_utils.os.path.isfile', return_value=True)
    @patch('vcf_agent.vcf_utils.VCF')
    def test_populate_kuzu_from_vcf_no_samples_in_vcf(self, mock_vcf_constructor, mock_isfile):
        """Test VCF processing when VCF has no samples and no override is given."""
        mock_kuzu_conn = MagicMock()
        vcf_file_path = "no_samples.vcf"

        mock_vcf_reader = MagicMock()
        mock_vcf_reader.samples = [] # Empty samples list
        # No records needed as it should fail before processing them
        mock_vcf_reader.__iter__.return_value = iter([]) 
        mock_vcf_constructor.return_value = mock_vcf_reader

        with pytest.raises(IndexError) as excinfo: # Expect IndexError from vcf_reader.samples[0]
            vcf_utils.populate_kuzu_from_vcf(mock_kuzu_conn, vcf_file_path)
        
        # Check that the error message is what you'd expect from trying to access samples[0]
        # This is a bit dependent on Python's default IndexError message for empty lists.
        # A more robust check might be to ensure no graph operations were attempted.
        # For now, checking the type is the primary goal.
        assert "list index out of range" in str(excinfo.value).lower()
        mock_vcf_constructor.assert_called_once_with(vcf_file_path)

    @patch('vcf_agent.vcf_utils.os.path.isfile', return_value=True)
    @patch('vcf_agent.vcf_utils.graph_integration')
    @patch('vcf_agent.vcf_utils.VCF')
    def test_populate_kuzu_from_vcf_multiple_alt_alleles(self, mock_vcf_constructor, mock_graph_integration, mock_isfile):
        """Test VCF processing ensuring only the first ALT allele is used."""
        mock_kuzu_conn = MagicMock()
        vcf_file_path = "multi_alt.vcf"
        sample_name = "test_sample_multi_alt"

        mock_vcf_reader = MagicMock()
        mock_vcf_reader.samples = [sample_name]
        
        # Record with two ALT alleles, only the first ('T') should be used
        record_multi_alt = MockVCFRecord(CHROM="chr2", POS=300, REF="C", ALT=["T", "G"], gt_types=[0])
        mock_vcf_reader.__iter__.return_value = iter([record_multi_alt])
        mock_vcf_constructor.return_value = mock_vcf_reader

        counts = vcf_utils.populate_kuzu_from_vcf(mock_kuzu_conn, vcf_file_path)

        mock_graph_integration.add_sample.assert_called_once_with(
            mock_kuzu_conn, 
            {"sample_id": sample_name}
        )
        # Ensure variant_id and alt field in add_variant use only the first ALT allele ('T')
        mock_graph_integration.add_variant.assert_called_once_with(
            mock_kuzu_conn, 
            {"variant_id": "chr2_300_C_T", "chrom": "chr2", "pos": 300, "ref": "C", "alt": "T"}
        )
        mock_graph_integration.link_variant_to_sample.assert_called_once_with(
            mock_kuzu_conn, 
            sample_name, 
            "chr2_300_C_T", # Uses variant_id based on first ALT
            {"zygosity": "HET"}
        )
        assert counts == {"variants": 1, "samples": 1, "links": 1}
        mock_vcf_constructor.assert_called_once_with(vcf_file_path)

    # TODO: Add tests for populate_kuzu_from_vcf
    #   - Basic processing (het, hom)
    #   - sample_name_override
    #   - FileNotFoundError
    #   - VCF with no samples
    #   - VCF with multiple ALT alleles 