"""
Test utilities package for VCF Analysis Agent.

Provides common utilities for testing including golden file management,
test data handling, and assertion helpers.
"""

from .golden_file_utils import (
    GoldenFileManager,
    GoldenFileError,
    assert_golden_match,
    generate_golden_file,
    load_golden_file,
    save_golden_file,
    normalize_output,
    compare_text_content,
    compare_json_content,
)

__all__ = [
    'GoldenFileManager',
    'GoldenFileError', 
    'assert_golden_match',
    'generate_golden_file',
    'load_golden_file',
    'save_golden_file',
    'normalize_output',
    'compare_text_content',
    'compare_json_content',
] 