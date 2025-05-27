"""
Golden File Testing Utilities

Provides utilities for loading, comparing, and validating golden files
in the VCF Analysis Agent test suite.
"""

import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from difflib import unified_diff
import tempfile


class GoldenFileError(Exception):
    """Exception raised for golden file operations."""
    pass


class GoldenFileManager:
    """Manager for golden file operations."""
    
    def __init__(self, golden_dir: Optional[Path] = None):
        """Initialize with golden file directory."""
        if golden_dir is None:
            # Default to golden/ directory in project root
            project_root = Path(__file__).parent.parent.parent
            golden_dir = project_root / "golden"
        
        self.golden_dir = Path(golden_dir)
        if not self.golden_dir.exists():
            raise GoldenFileError(f"Golden directory not found: {self.golden_dir}")
    
    def load_golden_file(self, relative_path: str) -> str:
        """Load a golden file and return its contents."""
        golden_path = self.golden_dir / relative_path
        
        if not golden_path.exists():
            raise GoldenFileError(f"Golden file not found: {golden_path}")
        
        try:
            with open(golden_path, 'r', encoding='utf-8') as f:
                return f.read()
        except Exception as e:
            raise GoldenFileError(f"Failed to read golden file {golden_path}: {e}")
    
    def save_golden_file(self, relative_path: str, content: str) -> None:
        """Save content to a golden file."""
        golden_path = self.golden_dir / relative_path
        
        # Create parent directories if they don't exist
        golden_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            with open(golden_path, 'w', encoding='utf-8') as f:
                f.write(content)
        except Exception as e:
            raise GoldenFileError(f"Failed to write golden file {golden_path}: {e}")
    
    def golden_file_exists(self, relative_path: str) -> bool:
        """Check if a golden file exists."""
        golden_path = self.golden_dir / relative_path
        return golden_path.exists()


def normalize_output(content: str, 
                    remove_timestamps: bool = True,
                    remove_paths: bool = True,
                    normalize_whitespace: bool = True) -> str:
    """
    Normalize output content for consistent comparison.
    
    Args:
        content: Raw output content
        remove_timestamps: Remove timestamp patterns
        remove_paths: Remove absolute file paths
        normalize_whitespace: Normalize whitespace and line endings
    
    Returns:
        Normalized content string
    """
    normalized = content
    
    if remove_timestamps:
        # Remove common timestamp patterns
        timestamp_patterns = [
            r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}',  # ISO format
            r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}',   # Standard format
            r'\w{3} \w{3} \d{1,2} \d{2}:\d{2}:\d{2} \d{4}',  # Date command format
        ]
        for pattern in timestamp_patterns:
            normalized = re.sub(pattern, '[TIMESTAMP]', normalized)
    
    if remove_paths:
        # Remove absolute paths, keep relative paths
        normalized = re.sub(r'/[^\s]+/', '[PATH]/', normalized)
        # Remove Windows paths
        normalized = re.sub(r'[A-Z]:\\\\[^\s]+\\\\', '[PATH]\\\\', normalized)
    
    if normalize_whitespace:
        # Normalize line endings and trailing whitespace
        lines = normalized.splitlines()
        lines = [line.rstrip() for line in lines]
        normalized = '\n'.join(lines)
        # Remove multiple consecutive empty lines
        normalized = re.sub(r'\n\s*\n\s*\n', '\n\n', normalized)
    
    return normalized.strip()


def compare_text_content(actual: str, 
                        expected: str,
                        tolerance: float = 0.0,
                        normalize: bool = True) -> tuple[bool, str]:
    """
    Compare text content with optional normalization and tolerance.
    
    Args:
        actual: Actual output content
        expected: Expected golden file content
        tolerance: Tolerance for floating-point comparisons (0.0 = exact match)
        normalize: Whether to normalize content before comparison
    
    Returns:
        Tuple of (is_match, diff_message)
    """
    if normalize:
        actual_norm = normalize_output(actual)
        expected_norm = normalize_output(expected)
    else:
        actual_norm = actual
        expected_norm = expected
    
    # Handle floating-point tolerance if specified
    if tolerance > 0.0:
        actual_norm = _apply_float_tolerance(actual_norm, expected_norm, tolerance)
    
    if actual_norm == expected_norm:
        return True, ""
    
    # Generate diff for debugging
    diff_lines = list(unified_diff(
        expected_norm.splitlines(keepends=True),
        actual_norm.splitlines(keepends=True),
        fromfile='expected (golden)',
        tofile='actual (current)',
        lineterm=''
    ))
    
    diff_message = ''.join(diff_lines)
    return False, diff_message


def compare_json_content(actual: Union[str, Dict], 
                        expected: Union[str, Dict],
                        tolerance: float = 0.0) -> tuple[bool, str]:
    """
    Compare JSON content with optional floating-point tolerance.
    
    Args:
        actual: Actual JSON content (string or dict)
        expected: Expected JSON content (string or dict)
        tolerance: Tolerance for floating-point comparisons
    
    Returns:
        Tuple of (is_match, diff_message)
    """
    try:
        if isinstance(actual, str):
            actual_data = json.loads(actual)
        else:
            actual_data = actual
        
        if isinstance(expected, str):
            expected_data = json.loads(expected)
        else:
            expected_data = expected
        
        if tolerance > 0.0:
            actual_data = _apply_json_float_tolerance(actual_data, expected_data, tolerance)
        
        if actual_data == expected_data:
            return True, ""
        
        # Generate readable diff
        actual_str = json.dumps(actual_data, indent=2, sort_keys=True)
        expected_str = json.dumps(expected_data, indent=2, sort_keys=True)
        
        diff_lines = list(unified_diff(
            expected_str.splitlines(keepends=True),
            actual_str.splitlines(keepends=True),
            fromfile='expected (golden)',
            tofile='actual (current)',
            lineterm=''
        ))
        
        diff_message = ''.join(diff_lines)
        return False, diff_message
        
    except json.JSONDecodeError as e:
        return False, f"JSON parsing error: {e}"


def assert_golden_match(actual: str,
                       golden_path: str,
                       tolerance: float = 0.0,
                       normalize: bool = True,
                       content_type: str = 'text',
                       golden_manager: Optional[GoldenFileManager] = None) -> None:
    """
    Assert that actual content matches golden file content.
    
    Args:
        actual: Actual output content
        golden_path: Relative path to golden file
        tolerance: Tolerance for floating-point comparisons
        normalize: Whether to normalize content before comparison
        content_type: Type of content ('text' or 'json')
        golden_manager: Optional golden file manager instance
    
    Raises:
        AssertionError: If content doesn't match
        GoldenFileError: If golden file operations fail
    """
    if golden_manager is None:
        golden_manager = GoldenFileManager()
    
    if not golden_manager.golden_file_exists(golden_path):
        raise GoldenFileError(f"Golden file not found: {golden_path}")
    
    expected = golden_manager.load_golden_file(golden_path)
    
    if content_type == 'json':
        is_match, diff_message = compare_json_content(actual, expected, tolerance)
    else:
        is_match, diff_message = compare_text_content(actual, expected, tolerance, normalize)
    
    if not is_match:
        raise AssertionError(f"Content doesn't match golden file {golden_path}:\n{diff_message}")


def generate_golden_file(operation_func,
                        golden_path: str,
                        *args,
                        golden_manager: Optional[GoldenFileManager] = None,
                        **kwargs) -> str:
    """
    Generate a golden file by running an operation and saving the output.
    
    Args:
        operation_func: Function to run to generate output
        golden_path: Relative path where to save golden file
        *args: Arguments to pass to operation_func
        golden_manager: Optional golden file manager instance
        **kwargs: Keyword arguments to pass to operation_func
    
    Returns:
        Generated content string
    """
    if golden_manager is None:
        golden_manager = GoldenFileManager()
    
    # Run the operation
    content = operation_func(*args, **kwargs)
    
    # Normalize the content
    normalized_content = normalize_output(content)
    
    # Save to golden file
    golden_manager.save_golden_file(golden_path, normalized_content)
    
    return normalized_content


def _apply_float_tolerance(actual: str, expected: str, tolerance: float) -> str:
    """Apply floating-point tolerance to text content."""
    # This is a simplified implementation
    # In practice, you might want more sophisticated float matching
    import re
    
    def replace_float(match):
        try:
            actual_val = float(match.group())
            # Find corresponding value in expected
            # This is a simplified approach
            return match.group()
        except ValueError:
            return match.group()
    
    # Find floating-point numbers and apply tolerance
    float_pattern = r'-?\d+\.\d+'
    return re.sub(float_pattern, replace_float, actual)


def _apply_json_float_tolerance(actual_data: Any, expected_data: Any, tolerance: float) -> Any:
    """Apply floating-point tolerance to JSON data recursively."""
    if isinstance(actual_data, dict) and isinstance(expected_data, dict):
        result = {}
        for key in actual_data:
            if key in expected_data:
                result[key] = _apply_json_float_tolerance(actual_data[key], expected_data[key], tolerance)
            else:
                result[key] = actual_data[key]
        return result
    elif isinstance(actual_data, list) and isinstance(expected_data, list):
        return [_apply_json_float_tolerance(a, e, tolerance) 
                for a, e in zip(actual_data, expected_data)]
    elif isinstance(actual_data, float) and isinstance(expected_data, float):
        if abs(actual_data - expected_data) <= tolerance:
            return expected_data  # Use expected value if within tolerance
        return actual_data
    else:
        return actual_data


# Convenience functions for common operations
def load_golden_file(relative_path: str) -> str:
    """Load a golden file (convenience function)."""
    manager = GoldenFileManager()
    return manager.load_golden_file(relative_path)


def save_golden_file(relative_path: str, content: str) -> None:
    """Save a golden file (convenience function)."""
    manager = GoldenFileManager()
    manager.save_golden_file(relative_path, content) 