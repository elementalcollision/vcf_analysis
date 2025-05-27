#!/usr/bin/env python3
"""
Golden File Generation Script

Generates golden files by running validated VCF operations on known-good test data.
This script should be run manually after verifying that the operations produce
correct outputs.
"""

import sys
import os
import json
from pathlib import Path

# Add the src and project root to the path so we can import our modules
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))
sys.path.insert(0, str(project_root))

from vcf_agent.bcftools_integration import (
    bcftools_view,
    bcftools_stats, 
    bcftools_query,
    bcftools_filter,
    bcftools_norm
)
from vcf_agent.validation import validate_vcf_file
from tests.utils.golden_file_utils import GoldenFileManager, normalize_output


def generate_bcftools_golden_files():
    """Generate golden files for bcftools operations."""
    print("Generating bcftools golden files...")
    
    # Test files to use
    test_files = [
        "sample_test_data/small_valid.vcf",
        "sample_test_data/sample1.vcf", 
        "sample_test_data/sample2.vcf"
    ]
    
    golden_manager = GoldenFileManager()
    
    for test_file in test_files:
        if not Path(test_file).exists():
            print(f"Skipping {test_file} - file not found")
            continue
            
        print(f"Processing {test_file}...")
        file_base = Path(test_file).stem
        
        try:
            # bcftools view (header only)
            print(f"  Generating view header for {file_base}")
            rc, stdout, stderr = bcftools_view(["-h", test_file])
            if rc == 0:
                content = normalize_output(stdout)
                golden_path = f"bcftools/view/view_header_{file_base}.txt"
                golden_manager.save_golden_file(golden_path, content)
                print(f"    Saved: {golden_path}")
            else:
                print(f"    Error: {stderr}")
            
            # bcftools stats
            print(f"  Generating stats for {file_base}")
            rc, stdout, stderr = bcftools_stats([test_file])
            if rc == 0:
                content = normalize_output(stdout)
                golden_path = f"bcftools/stats/stats_{file_base}.txt"
                golden_manager.save_golden_file(golden_path, content)
                print(f"    Saved: {golden_path}")
            else:
                print(f"    Error: {stderr}")
            
            # bcftools query (basic format)
            print(f"  Generating query for {file_base}")
            rc, stdout, stderr = bcftools_query(["-f", "%CHROM\t%POS\t%REF\t%ALT\n", test_file])
            if rc == 0:
                content = normalize_output(stdout)
                golden_path = f"bcftools/query/query_basic_{file_base}.txt"
                golden_manager.save_golden_file(golden_path, content)
                print(f"    Saved: {golden_path}")
            else:
                print(f"    Error: {stderr}")
                
        except Exception as e:
            print(f"  Error processing {test_file}: {e}")


def generate_validation_golden_files():
    """Generate golden files for VCF validation operations."""
    print("Generating validation golden files...")
    
    golden_manager = GoldenFileManager()
    
    # Valid files
    valid_files = [
        "sample_test_data/small_valid.vcf",
        "sample_test_data/sample1.vcf"
    ]
    
    for test_file in valid_files:
        if not Path(test_file).exists():
            print(f"Skipping {test_file} - file not found")
            continue
            
        print(f"Validating {test_file}...")
        file_base = Path(test_file).stem
        
        try:
            validation_result = validate_vcf_file(test_file)
            
            # Convert result to JSON for consistent storage
            result_json = json.dumps(validation_result, indent=2, sort_keys=True)
            golden_path = f"validation/valid_files/validation_{file_base}.json"
            golden_manager.save_golden_file(golden_path, result_json)
            print(f"  Saved: {golden_path}")
            
        except Exception as e:
            print(f"  Error validating {test_file}: {e}")
    
    # Test invalid file validation (if we have any)
    invalid_test_file = "sample_test_data/invalid_test.vcf"
    if Path(invalid_test_file).exists():
        print(f"Validating invalid file {invalid_test_file}...")
        try:
            validation_result = validate_vcf_file(invalid_test_file)
            result_json = json.dumps(validation_result, indent=2, sort_keys=True)
            golden_path = "validation/invalid_files/validation_invalid_test.json"
            golden_manager.save_golden_file(golden_path, result_json)
            print(f"  Saved: {golden_path}")
        except Exception as e:
            print(f"  Error validating invalid file: {e}")


def generate_comparison_golden_files():
    """Generate golden files for VCF comparison operations."""
    print("Generating comparison golden files...")
    
    # This would involve running VCF comparison operations
    # For now, we'll create placeholder structure
    golden_manager = GoldenFileManager()
    
    # Placeholder for comparison results
    comparison_result = {
        "comparison_type": "basic",
        "file1": "sample1.vcf",
        "file2": "sample2.vcf", 
        "matches": 0,
        "differences": 0,
        "summary": "Comparison completed"
    }
    
    result_json = json.dumps(comparison_result, indent=2, sort_keys=True)
    golden_path = "comparison/basic/comparison_sample1_sample2.json"
    golden_manager.save_golden_file(golden_path, result_json)
    print(f"  Saved placeholder: {golden_path}")


def generate_analysis_golden_files():
    """Generate golden files for agent analysis operations."""
    print("Generating analysis golden files...")
    
    golden_manager = GoldenFileManager()
    
    # Placeholder for analysis results
    analysis_result = """# VCF Analysis Summary

## File Information
- File: sample1.vcf
- Format: VCF 4.2
- Samples: 1

## Variant Summary
- Total variants: 5
- SNPs: 3
- Indels: 2

## Quality Metrics
- Average quality: 30.5
- Variants passing filters: 5

## Recommendations
- File appears to be valid and well-formatted
- All variants pass quality thresholds
"""
    
    golden_path = "analysis/summary/analysis_summary_sample1.md"
    golden_manager.save_golden_file(golden_path, analysis_result)
    print(f"  Saved placeholder: {golden_path}")


def main():
    """Main function to generate all golden files."""
    print("=== Golden File Generation Script ===")
    print("This script generates golden files from validated VCF operations.")
    print("Make sure to review all outputs before committing to version control.")
    print()
    
    # Check if golden directory exists
    golden_dir = Path("golden")
    if not golden_dir.exists():
        print("Error: golden/ directory not found. Please run from project root.")
        sys.exit(1)
    
    try:
        generate_bcftools_golden_files()
        print()
        
        generate_validation_golden_files()
        print()
        
        generate_comparison_golden_files()
        print()
        
        generate_analysis_golden_files()
        print()
        
        print("=== Golden file generation completed ===")
        print("Please review all generated files before committing:")
        print("  1. Check that outputs are correct and expected")
        print("  2. Verify that sensitive data is not included")
        print("  3. Test the golden files with the test suite")
        print("  4. Commit to version control after validation")
        
    except Exception as e:
        print(f"Error during golden file generation: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 