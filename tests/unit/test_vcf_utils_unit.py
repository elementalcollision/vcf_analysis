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
        self.genotypes = gt_types # RENAMED from gt_types to genotypes

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
        
        record1_het = MockVCFRecord(CHROM="chr1", POS=100, REF="A", ALT=["T"], gt_types=[[0,1,True]]) # Het
        record2_hom = MockVCFRecord(CHROM="chr1", POS=200, REF="G", ALT=["C"], gt_types=[[1,1,True]]) # Hom Alt
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
            call(mock_kuzu_conn, {"variant_id": "chr1-100-A-T", "chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "rs_id": None}),
            call(mock_kuzu_conn, {"variant_id": "chr1-200-G-C", "chrom": "chr1", "pos": 200, "ref": "G", "alt": "C", "rs_id": None})
        ]
        mock_graph_integration.add_variant.assert_has_calls(expected_variant_calls)
        assert mock_graph_integration.add_variant.call_count == 2

        # 3. Link additions
        expected_link_calls = [
            call(mock_kuzu_conn, "test_sample_from_vcf", "chr1-100-A-T", {"zygosity": "HET"}),
            call(mock_kuzu_conn, "test_sample_from_vcf", "chr1-200-G-C", {"zygosity": "HOM_ALT"})
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
        # For HET, gt_types could be e.g. [[0, 1, True]] (phased) or [[0, 1, False]] (unphased)
        # Let's use [[0,1,True]] for consistency
        record1 = MockVCFRecord(CHROM="chrX", POS=500, REF="C", ALT=["G"], gt_types=[[0,1,True]])
        mock_vcf_reader.__iter__.return_value = iter([record1])
        mock_vcf_constructor.return_value = mock_vcf_reader

        counts = vcf_utils.populate_kuzu_from_vcf(mock_kuzu_conn, vcf_file_path, sample_name_override=override_sample_name)

        mock_graph_integration.add_sample.assert_called_once_with(
            mock_kuzu_conn, 
            {"sample_id": override_sample_name}
        )
        mock_graph_integration.add_variant.assert_called_once_with(
            mock_kuzu_conn, 
            {"variant_id": "chrX-500-C-G", "chrom": "chrX", "pos": 500, "ref": "C", "alt": "G", "rs_id": None}
        )
        mock_graph_integration.link_variant_to_sample.assert_called_once_with(
            mock_kuzu_conn,
            override_sample_name,
            "chrX-500-C-G",
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
    @patch('vcf_agent.vcf_utils.graph_integration')
    @patch('vcf_agent.vcf_utils.VCF')
    def test_populate_kuzu_from_vcf_no_samples_in_vcf(self, mock_vcf_constructor, mock_graph_integration, mock_isfile):
        """Test VCF processing when VCF has no samples and no override is given."""
        mock_kuzu_conn = MagicMock()
        vcf_file_path = "no_samples.vcf"

        mock_vcf_reader = MagicMock()
        mock_vcf_reader.samples = [] # Empty samples list
        # No records needed as it should fail before processing them
        mock_vcf_reader.__iter__.return_value = iter([]) 
        mock_vcf_constructor.return_value = mock_vcf_reader

        # With no samples in VCF and no override, no samples should be processed,
        # no records iterated, and thus no variants/links. No IndexError expected here.
        counts = vcf_utils.populate_kuzu_from_vcf(mock_kuzu_conn, vcf_file_path)
        
        assert counts == {"variants": 0, "samples": 0, "links": 0}
        # Ensure graph_integration calls were not made if no samples/records processed
        mock_vcf_constructor.assert_called_once_with(vcf_file_path)
        mock_graph_integration.add_sample.assert_not_called()
        mock_graph_integration.add_variant.assert_not_called()
        mock_graph_integration.link_variant_to_sample.assert_not_called()
        mock_isfile.assert_called_once_with(vcf_file_path)

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
        
        # Record with two ALT alleles, only the first ('T') should be used. Genotype is HET [0,1,True]
        record_multi_alt = MockVCFRecord(CHROM="chr2", POS=300, REF="C", ALT=["T", "G"], gt_types=[[0,1,True]])
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
            {"variant_id": "chr2-300-C-T", "chrom": "chr2", "pos": 300, "ref": "C", "alt": "T", "rs_id": None}
        )
        mock_graph_integration.link_variant_to_sample.assert_called_once_with(
            mock_kuzu_conn,
            sample_name,
            "chr2-300-C-T", # Uses variant_id based on first ALT
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