"""
pytest configuration and shared fixtures for CLI Enhanced Validation Engine tests.

This file contains shared pytest fixtures and configuration based on best practices
from pytest-mock research and CLI testing patterns.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch
from typing import Generator, Dict, Any

# Import the components we're testing
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.cli_enhanced_validation import (
    CacheManager,
    GitIntegration,
    ValidationConfig,
    ASTAnalyzer,
    MultiFormatDocstringParser,
    PerformanceOptimizer,
    EnhancedCLIValidator
)


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test files."""
    temp_path = Path(tempfile.mkdtemp())
    try:
        yield temp_path
    finally:
        shutil.rmtree(temp_path, ignore_errors=True)


@pytest.fixture
def sample_python_file(temp_dir: Path) -> Path:
    """Create a sample Python file with various docstring formats for testing."""
    content = '''"""Module docstring for testing."""

def function_with_google_docstring(param1: str, param2: int = 10) -> str:
    """Function with Google-style docstring.
    
    Args:
        param1: First parameter description
        param2: Second parameter with default value
        
    Returns:
        A string result
        
    Raises:
        ValueError: If param1 is empty
    """
    if not param1:
        raise ValueError("param1 cannot be empty")
    return f"{param1}_{param2}"


def function_with_numpy_docstring(data: list) -> float:
    """Function with NumPy-style docstring.
    
    Parameters
    ----------
    data : list
        Input data list
        
    Returns
    -------
    float
        Mean of the data
    """
    return sum(data) / len(data) if data else 0.0


def function_without_docstring():
    return "no docs"


class TestClass:
    """Test class with methods."""
    
    def method_with_docstring(self, value: int) -> int:
        """Method with docstring.
        
        Args:
            value: Input value
            
        Returns:
            Modified value
        """
        return value * 2
        
    def method_without_docstring(self):
        return True
'''
    
    file_path = temp_dir / "sample_module.py"
    file_path.write_text(content)
    return file_path


@pytest.fixture
def sample_config_yaml(temp_dir: Path) -> Path:
    """Create a sample YAML configuration file."""
    content = '''validation:
  level: comprehensive
  formats: [google, numpy, rest]
  cache: true
  parallel: true

rules:
  coverage_threshold: 95
  require_examples: true
  validate_types: true

ci:
  format: json
  github_annotations: true
  fail_on_warnings: false

cache:
  ttl_hours: 24
  max_size_mb: 100
  enabled: true

performance:
  parallel_workers: 4
  timeout_seconds: 300
'''
    
    file_path = temp_dir / "test_config.yml"
    file_path.write_text(content)
    return file_path


@pytest.fixture
def sample_config_toml(temp_dir: Path) -> Path:
    """Create a sample TOML configuration file."""
    content = '''[validation]
level = "comprehensive"
formats = ["google", "numpy", "rest"]
cache = true
parallel = true

[rules]
coverage_threshold = 95
require_examples = true
validate_types = true

[ci]
format = "json"
github_annotations = true
fail_on_warnings = false

[cache]
ttl_hours = 24
max_size_mb = 100
enabled = true

[performance]
parallel_workers = 4
timeout_seconds = 300
'''
    
    file_path = temp_dir / "test_config.toml"
    file_path.write_text(content)
    return file_path


@pytest.fixture
def mock_git_repo(temp_dir: Path) -> Path:
    """Create a mock Git repository for testing Git integration."""
    git_dir = temp_dir / ".git"
    git_dir.mkdir()
    
    # Create some sample files
    (temp_dir / "file1.py").write_text("# Sample file 1")
    (temp_dir / "file2.py").write_text("# Sample file 2")
    (temp_dir / "README.md").write_text("# Test repo")
    
    return temp_dir


@pytest.fixture
def cache_manager(temp_dir: Path) -> CacheManager:
    """Create a CacheManager instance with temporary cache directory."""
    cache_dir = temp_dir / "test_cache"
    return CacheManager(cache_dir=str(cache_dir), ttl_hours=1, max_size_mb=10)


@pytest.fixture
def git_integration() -> GitIntegration:
    """Create a GitIntegration instance."""
    return GitIntegration()


@pytest.fixture
def validation_config() -> ValidationConfig:
    """Create a ValidationConfig instance with default settings."""
    return ValidationConfig()


@pytest.fixture
def ast_analyzer() -> ASTAnalyzer:
    """Create an ASTAnalyzer instance."""
    return ASTAnalyzer()


@pytest.fixture
def docstring_parser() -> MultiFormatDocstringParser:
    """Create a MultiFormatDocstringParser instance."""
    return MultiFormatDocstringParser()


@pytest.fixture
def performance_optimizer() -> PerformanceOptimizer:
    """Create a PerformanceOptimizer instance."""
    config = ValidationConfig()
    return PerformanceOptimizer(config)


@pytest.fixture
def enhanced_validator(temp_dir: Path) -> EnhancedCLIValidator:
    """Create an EnhancedCLIValidator instance for testing."""
    # Create a temporary config file
    config_file = temp_dir / "test_config.yml"
    config_content = '''validation:
  level: comprehensive
cache:
  enabled: true
  ttl_hours: 1
performance:
  parallel_workers: 2
'''
    config_file.write_text(config_content)
    return EnhancedCLIValidator(config_file=str(config_file))


@pytest.fixture
def mock_subprocess():
    """Mock subprocess operations for Git commands."""
    with patch('subprocess.run') as mock_run:
        # Default successful Git command
        mock_run.return_value = Mock(
            returncode=0,
            stdout="",
            stderr=""
        )
        yield mock_run


@pytest.fixture
def performance_baseline() -> Dict[str, Any]:
    """Baseline performance metrics for regression testing."""
    return {
        'cache_hit_rate_threshold': 0.8,  # 80% minimum cache hit rate
        'validation_speed_threshold': 1000,  # Max 1000ms for 100 definitions
        'memory_usage_threshold': 100,  # Max 100MB memory usage
        'parallel_speedup_threshold': 1.5,  # Min 1.5x speedup with parallelization
    }


# Pytest configuration
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "performance: marks tests as performance benchmarks"
    )
    config.addinivalue_line(
        "markers", "unit: marks tests as unit tests"
    )


def pytest_collection_modifyitems(config, items):
    """Automatically mark test types based on directory structure."""
    for item in items:
        # Mark tests based on directory
        if "unit" in str(item.fspath):
            item.add_marker(pytest.mark.unit)
        elif "integration" in str(item.fspath):
            item.add_marker(pytest.mark.integration)
        elif "performance" in str(item.fspath):
            item.add_marker(pytest.mark.performance)
            item.add_marker(pytest.mark.slow) 