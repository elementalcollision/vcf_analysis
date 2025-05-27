"""
Tests for SAMspec compliance validation module.

Tests the SAMspecValidator class and related functionality for validating
VCF files against the SAM/VCF specification standards.
"""

import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import patch, mock_open

from vcf_agent.samspec_compliance import (
    SAMspecValidator,
    validate_vcf_samspec_compliance,
    ComplianceLevel,
    ComplianceViolation,
    ComplianceReport
)


class TestSAMspecValidator:
    """Test the SAMspecValidator class."""
    
    def test_validator_initialization(self):
        """Test validator initializes correctly."""
        validator = SAMspecValidator()
        assert validator.violations == []
        assert validator.vcf_version is None
        assert validator.header_lines == []
        assert validator.data_lines == []
        assert validator.info_fields == {}
        assert validator.format_fields == {}
        assert validator.filter_fields == {}
        assert validator.contig_fields == {}
        assert validator.samples == []
    
    def test_validate_valid_vcf_file(self):
        """Test validation of a valid VCF file."""
        valid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	12345	rs123	A	G	30.0	PASS	DP=10	GT	0/1
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(valid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.file_path == f.name
                assert report.vcf_version == "4.2"
                assert report.is_compliant == True
                assert report.critical_count == 0
                assert len(report.violations) == 0
                
            finally:
                os.unlink(f.name)
    
    def test_validate_missing_fileformat(self):
        """Test validation detects missing ##fileformat line."""
        invalid_vcf_content = """##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	A	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Check for missing fileformat violation
                missing_format_violations = [
                    v for v in report.violations 
                    if v.rule_id == "MISSING_FILEFORMAT"
                ]
                assert len(missing_format_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_invalid_fileformat(self):
        """Test validation detects invalid ##fileformat line."""
        invalid_vcf_content = """##fileformat=InvalidFormat
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	A	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Check for invalid fileformat violation
                invalid_format_violations = [
                    v for v in report.violations 
                    if v.rule_id == "INVALID_FILEFORMAT"
                ]
                assert len(invalid_format_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_missing_chrom_line(self):
        """Test validation detects missing #CHROM header line."""
        invalid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
1	12345	rs123	A	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Check for missing CHROM line violation
                missing_chrom_violations = [
                    v for v in report.violations 
                    if v.rule_id == "MISSING_CHROM_LINE"
                ]
                assert len(missing_chrom_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_invalid_chrom_fields(self):
        """Test validation detects invalid #CHROM line fields."""
        invalid_vcf_content = """##fileformat=VCFv4.2
#CHROM	POSITION	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	A	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Check for invalid CHROM field violation
                invalid_field_violations = [
                    v for v in report.violations 
                    if v.rule_id == "INVALID_CHROM_FIELD"
                ]
                assert len(invalid_field_violations) > 0
                
            finally:
                os.unlink(f.name)
    
    def test_validate_invalid_info_definition(self):
        """Test validation detects invalid INFO definitions."""
        invalid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=InvalidType,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	A	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Check for invalid INFO type violation
                invalid_type_violations = [
                    v for v in report.violations 
                    if v.rule_id == "INVALID_INFO_TYPE"
                ]
                assert len(invalid_type_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_missing_info_fields(self):
        """Test validation detects missing required INFO fields."""
        invalid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	A	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Check for missing INFO field violation
                missing_field_violations = [
                    v for v in report.violations 
                    if v.rule_id == "MISSING_INFO_FIELD"
                ]
                assert len(missing_field_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_invalid_data_record_fields(self):
        """Test validation detects invalid data record fields."""
        invalid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	invalid_pos	rs123	A	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Check for invalid POS format violation
                invalid_pos_violations = [
                    v for v in report.violations 
                    if v.rule_id == "INVALID_POS_FORMAT"
                ]
                assert len(invalid_pos_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_invalid_ref_field(self):
        """Test validation detects invalid REF field."""
        invalid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	.	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Check for invalid REF violation
                invalid_ref_violations = [
                    v for v in report.violations 
                    if v.rule_id == "INVALID_REF"
                ]
                assert len(invalid_ref_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_invalid_alt_bases(self):
        """Test validation detects invalid ALT bases."""
        invalid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	A	X	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Check for invalid ALT bases violation
                invalid_alt_violations = [
                    v for v in report.violations 
                    if v.rule_id == "INVALID_ALT_BASES"
                ]
                assert len(invalid_alt_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_undefined_filter(self):
        """Test validation detects undefined FILTER values."""
        invalid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	A	G	30.0	UNDEFINED_FILTER	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                # This should be a warning, not critical
                assert report.warning_count > 0
                
                # Check for undefined filter violation
                undefined_filter_violations = [
                    v for v in report.violations 
                    if v.rule_id == "UNDEFINED_FILTER"
                ]
                assert len(undefined_filter_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_compressed_file(self):
        """Test validation works with gzip-compressed VCF files."""
        import gzip
        
        valid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	A	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False) as f:
            with gzip.open(f.name, 'wt') as gz_file:
                gz_file.write(valid_vcf_content)
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.vcf_version == "4.2"
                assert report.is_compliant == True
                assert report.critical_count == 0
                
            finally:
                os.unlink(f.name)
    
    def test_validate_file_not_found(self):
        """Test validation handles file not found error."""
        report = validate_vcf_samspec_compliance("nonexistent_file.vcf")
        
        assert report.is_compliant == False
        assert report.critical_count > 0
        
        # Check for file read error
        file_error_violations = [
            v for v in report.violations 
            if v.rule_id == "FILE_READ_ERROR"
        ]
        assert len(file_error_violations) == 1
    
    def test_compliance_report_to_dict(self):
        """Test ComplianceReport to_dict method."""
        violation = ComplianceViolation(
            level=ComplianceLevel.CRITICAL,
            rule_id="TEST_RULE",
            message="Test message",
            line_number=1,
            field="TEST_FIELD",
            value="test_value",
            suggestion="Test suggestion"
        )
        
        report = ComplianceReport(
            file_path="test.vcf",
            vcf_version="4.2",
            total_violations=1,
            critical_count=1,
            warning_count=0,
            info_count=0,
            violations=[violation],
            is_compliant=False
        )
        
        report_dict = report.to_dict()
        
        assert report_dict["file_path"] == "test.vcf"
        assert report_dict["vcf_version"] == "4.2"
        assert report_dict["total_violations"] == 1
        assert report_dict["critical_count"] == 1
        assert report_dict["warning_count"] == 0
        assert report_dict["info_count"] == 0
        assert report_dict["is_compliant"] == False
        assert len(report_dict["violations"]) == 1
        
        violation_dict = report_dict["violations"][0]
        assert violation_dict["level"] == "critical"
        assert violation_dict["rule_id"] == "TEST_RULE"
        assert violation_dict["message"] == "Test message"
        assert violation_dict["line_number"] == 1
        assert violation_dict["field"] == "TEST_FIELD"
        assert violation_dict["value"] == "test_value"
        assert violation_dict["suggestion"] == "Test suggestion"


class TestSAMspecValidatorEdgeCases:
    """Test edge cases and complex scenarios."""
    
    def test_validate_format_sample_mismatch(self):
        """Test validation detects FORMAT/sample data mismatch."""
        invalid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	12345	rs123	A	G	30.0	PASS	DP=10	GT:DP	0/1
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                # This should generate a warning about format/sample mismatch
                assert report.warning_count > 0
                
                # Check for format sample mismatch violation
                mismatch_violations = [
                    v for v in report.violations 
                    if v.rule_id == "FORMAT_SAMPLE_MISMATCH"
                ]
                assert len(mismatch_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_multiple_violations(self):
        """Test validation detects multiple violations in one file."""
        invalid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Type=Integer,Description="Total Depth">
#CHROM	POSITION	ID	REF	ALT	QUAL	FILTER	INFO
1	invalid_pos	rs123	.	X	-10.0	UNDEFINED_FILTER	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(invalid_vcf_content)
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.total_violations > 3  # Multiple violations expected
                
                # Should have critical violations
                assert report.critical_count > 0
                
                # Should have warnings
                assert report.warning_count > 0
                
            finally:
                os.unlink(f.name)
    
    def test_validate_empty_file(self):
        """Test validation handles empty VCF file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write("")  # Empty file
            f.flush()
            
            try:
                report = validate_vcf_samspec_compliance(f.name)
                
                assert report.is_compliant == False
                assert report.critical_count > 0
                
                # Should detect missing fileformat
                missing_format_violations = [
                    v for v in report.violations 
                    if v.rule_id == "MISSING_FILEFORMAT"
                ]
                assert len(missing_format_violations) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_validate_header_content_parsing(self):
        """Test parsing of complex header content."""
        validator = SAMspecValidator()
        
        # Test parsing with quoted values
        content = 'ID=TEST,Number=1,Type=String,Description="Test with spaces and, commas"'
        result = validator._parse_header_content(content)
        
        assert result["ID"] == "TEST"
        assert result["Number"] == "1"
        assert result["Type"] == "String"
        assert result["Description"] == "Test with spaces and, commas"
        
        # Test parsing without quotes
        content = 'ID=SIMPLE,Number=A,Type=Integer,Description=SimpleDescription'
        result = validator._parse_header_content(content)
        
        assert result["ID"] == "SIMPLE"
        assert result["Number"] == "A"
        assert result["Type"] == "Integer"
        assert result["Description"] == "SimpleDescription"


class TestSAMspecValidatorPerformance:
    """Test performance aspects of the validator."""
    
    def test_validate_large_file_simulation(self):
        """Test validation performance with simulated large file."""
        # Create a VCF with many data lines
        header = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
"""
        
        # Generate 1000 data lines
        data_lines = []
        for i in range(1000):
            pos = 10000 + i
            data_lines.append(f"1	{pos}	rs{i}	A	G	30.0	PASS	DP=10	GT	0/1")
        
        vcf_content = header + "\n".join(data_lines)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(vcf_content)
            f.flush()
            
            try:
                import time
                start_time = time.time()
                report = validate_vcf_samspec_compliance(f.name)
                end_time = time.time()
                
                # Validation should complete in reasonable time (< 5 seconds)
                assert (end_time - start_time) < 5.0
                
                # File should be compliant
                assert report.is_compliant == True
                assert report.critical_count == 0
                
            finally:
                os.unlink(f.name)


class TestConvenienceFunction:
    """Test the convenience function."""
    
    def test_validate_vcf_samspec_compliance_function(self):
        """Test the convenience function works correctly."""
        valid_vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs123	A	G	30.0	PASS	DP=10
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(valid_vcf_content)
            f.flush()
            
            try:
                # Test the convenience function
                report = validate_vcf_samspec_compliance(f.name)
                
                assert isinstance(report, ComplianceReport)
                assert report.file_path == f.name
                assert report.vcf_version == "4.2"
                assert report.is_compliant == True
                
            finally:
                os.unlink(f.name) 