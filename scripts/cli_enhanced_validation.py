#!/usr/bin/env python3
"""
Enhanced CLI Documentation Validation Engine

Comprehensive validation engine that integrates existing CLI documentation tools 
with advanced AST parsing, multi-format docstring support, and CI/CD integration.

This engine serves as the unified entry point for all CLI documentation validation,
providing:
- Integration with existing validation tools
- Advanced AST-based analysis
- Multi-format docstring parsing
- Performance optimization with caching
- CI/CD integration with proper exit codes and reporting

Usage:
    python scripts/cli_enhanced_validation.py [--mode=comprehensive] [--format=console] [--verbose]
    
Exit codes:
    0: All validation passed
    1: Validation failures found
    2: Script execution error
"""

import ast
import sys
import os
import time
import json
import logging
import hashlib
import argparse
import pickle
import threading
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Set, Callable, TYPE_CHECKING
from dataclasses import dataclass, field
from enum import Enum
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import existing validation tools
from scripts.validate_cli_documentation import CLIDocumentationValidator, ValidationResult as ModuleValidationResult
from scripts.cli_documentation_style_guide import CLIDocstringValidator, CLIDocumentationStandards

# Import docstring_parser for structured parsing
try:
    from docstring_parser import parse as parse_docstring
    DOCSTRING_PARSER_AVAILABLE = True
except ImportError:
    DOCSTRING_PARSER_AVAILABLE = False
    parse_docstring = None

# Optional imports for configuration management
if TYPE_CHECKING:
    # Type checking imports - these are for IDE/linter support only
    from pydantic_settings import BaseSettings, SettingsConfigDict
    from pydantic import Field
    import tomli
    import tomllib

# Runtime imports with graceful fallback
PYDANTIC_SETTINGS_AVAILABLE = False
try:
    from pydantic import BaseSettings, Field
    from pydantic.v1.config import Extra
    PYDANTIC_SETTINGS_AVAILABLE = True
except ImportError:
    try:
        # Try newer pydantic-settings package
        from pydantic_settings import BaseSettings, SettingsConfigDict  # type: ignore
        from pydantic import Field  # type: ignore
        PYDANTIC_SETTINGS_AVAILABLE = True
    except ImportError:
        # No pydantic support available
        PYDANTIC_SETTINGS_AVAILABLE = False

# TOML support detection
TOML_SUPPORT = False
try:
    # Python 3.11+ built-in support
    import tomllib  # type: ignore
    TOML_SUPPORT = True
    TOML_LOADER = tomllib
except ImportError:
    try:
        # Python < 3.11 fallback
        import tomli as TOML_LOADER  # type: ignore
        TOML_SUPPORT = True
    except ImportError:
        TOML_SUPPORT = False
        TOML_LOADER = None


# Validation modes
class ValidationMode(Enum):
    QUICK = "quick"
    COMPREHENSIVE = "comprehensive" 
    STRICT = "strict"


class OutputFormat(Enum):
    CONSOLE = "console"
    JSON = "json"
    GITHUB = "github"


@dataclass
class CodeDefinition:
    """Represents a code definition (function, method, class) found in source code."""
    name: str
    type: str  # 'function', 'method', 'class'
    file_path: str
    line_number: int
    docstring: Optional[str]
    parent_class: Optional[str] = None
    args: List[str] = field(default_factory=list)
    is_cli_handler: bool = False


@dataclass
class CoverageReport:
    """Docstring coverage analysis report."""
    total_definitions: int
    documented_definitions: int
    undocumented_definitions: int
    coverage_percentage: float
    undocumented_items: List[CodeDefinition] = field(default_factory=list)


@dataclass
class ValidationLocation:
    """Represents the location of a validation issue."""
    file_path: str
    line_number: int
    column_number: Optional[int] = None


@dataclass 
class ValidationIssue:
    """Represents a validation issue found during analysis (updated structure)."""
    type: str  # Issue type (e.g., 'missing_docstring', 'parameter_mismatch')
    severity: str  # 'error', 'warning', 'info'
    message: str
    location: ValidationLocation
    suggestions: List[str] = field(default_factory=list)


@dataclass
class FileValidationResult:
    """Results of validating a single file."""
    file_path: str
    definitions: List[CodeDefinition]
    issues: List[ValidationIssue]
    total_definitions: int
    documented_count: int
    coverage_percentage: float
    execution_time: float
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ValidationReport:
    """Comprehensive validation report."""
    total_files: int
    total_definitions: int
    coverage_percentage: float
    execution_time: float
    issues: List[ValidationIssue]
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Parameter:
    """Represents a function parameter from docstring or signature."""
    name: str
    type_hint: Optional[str] = None
    description: Optional[str] = None
    is_optional: bool = False
    default_value: Optional[str] = None


@dataclass
class ReturnInfo:
    """Represents return value information from docstring."""
    type_hint: Optional[str] = None
    description: Optional[str] = None


@dataclass
class ExceptionInfo:
    """Represents exception information from docstring."""
    type_name: str
    description: Optional[str] = None


@dataclass
class ParsedDocstring:
    """Represents a parsed docstring with structured components."""
    raw_text: str
    detected_format: str
    summary: str
    description: str
    parameters: List[Parameter] = field(default_factory=list)
    returns: Optional[ReturnInfo] = None
    raises: List[ExceptionInfo] = field(default_factory=list)
    examples: List[str] = field(default_factory=list)
    parsing_successful: bool = True
    parsing_errors: List[str] = field(default_factory=list)


@dataclass
class FunctionSignature:
    """Represents a function signature extracted from AST."""
    parameters: List[Parameter]
    return_type: Optional[str] = None


@dataclass
class StructureValidation:
    """Results of structured docstring validation."""
    is_valid: bool
    missing_parameters: List[str] = field(default_factory=list)
    extra_parameters: List[str] = field(default_factory=list)
    missing_return_doc: bool = False
    undocumented_exceptions: List[str] = field(default_factory=list)
    type_mismatches: List[str] = field(default_factory=list)
    format_issues: List[str] = field(default_factory=list)
    suggestions: List[str] = field(default_factory=list)


class CacheManager:
    """
    File-based caching system for validation results with content-hash invalidation.
    
    Provides efficient caching of AST parsing results and validation outcomes
    with automatic cache invalidation when file content changes.
    """
    
    def __init__(self, cache_dir: str = ".validation-cache", ttl_hours: int = 24, max_size_mb: int = 100):
        self.cache_dir = Path(cache_dir)
        self.ttl_seconds = ttl_hours * 3600
        self.max_size_bytes = max_size_mb * 1024 * 1024
        self.lock = threading.Lock()
        self.logger = logging.getLogger(__name__ + '.CacheManager')
        
        # Ensure cache directory exists
        self.cache_dir.mkdir(exist_ok=True)
        
        # Cache hit/miss statistics
        self.hits = 0
        self.misses = 0
        
        self.logger.debug(f"CacheManager initialized: dir={cache_dir}, ttl={ttl_hours}h, max_size={max_size_mb}MB")
    
    def _get_file_hash(self, file_path: str) -> str:
        """Generate content hash for a file."""
        try:
            with open(file_path, 'rb') as f:
                content = f.read()
            return hashlib.sha256(content).hexdigest()[:16]  # Use first 16 chars for performance
        except Exception as e:
            self.logger.debug(f"Failed to hash {file_path}: {e}")
            return "invalid"
    
    def _get_cache_key(self, file_path: str, validation_type: str = "ast") -> str:
        """Generate cache key for a file and validation type."""
        file_hash = self._get_file_hash(file_path)
        path_normalized = str(Path(file_path).resolve())
        key_content = f"{path_normalized}:{validation_type}:{file_hash}"
        return hashlib.md5(key_content.encode()).hexdigest()
    
    def _get_cache_file_path(self, cache_key: str) -> Path:
        """Get the cache file path for a given cache key."""
        return self.cache_dir / f"{cache_key}.cache"
    
    def get(self, file_path: str, validation_type: str = "ast") -> Optional[Any]:
        """
        Retrieve cached validation result for a file.
        
        Args:
            file_path: Path to the file being validated
            validation_type: Type of validation (ast, module, style)
            
        Returns:
            Cached result if valid, None if cache miss
        """
        with self.lock:
            try:
                cache_key = self._get_cache_key(file_path, validation_type)
                cache_file = self._get_cache_file_path(cache_key)
                
                if not cache_file.exists():
                    self.misses += 1
                    return None
                
                # Check if cache entry has expired
                cache_age = time.time() - cache_file.stat().st_mtime
                if cache_age > self.ttl_seconds:
                    cache_file.unlink()  # Remove expired cache
                    self.misses += 1
                    return None
                
                # Load cached result
                with open(cache_file, 'rb') as f:
                    cached_result = pickle.load(f)
                
                self.hits += 1
                self.logger.debug(f"Cache HIT: {file_path} ({validation_type})")
                return cached_result
                
            except Exception as e:
                self.logger.debug(f"Cache retrieval failed for {file_path}: {e}")
                self.misses += 1
                return None
    
    def set(self, file_path: str, result: Any, validation_type: str = "ast") -> bool:
        """
        Store validation result in cache.
        
        Args:
            file_path: Path to the file being validated
            result: Validation result to cache
            validation_type: Type of validation (ast, module, style)
            
        Returns:
            True if successfully cached, False otherwise
        """
        with self.lock:
            try:
                cache_key = self._get_cache_key(file_path, validation_type)
                cache_file = self._get_cache_file_path(cache_key)
                
                # Store result in cache
                with open(cache_file, 'wb') as f:
                    pickle.dump(result, f)
                
                self.logger.debug(f"Cache SET: {file_path} ({validation_type})")
                
                # Trigger cleanup if cache is getting large
                self._cleanup_if_needed()
                
                return True
                
            except Exception as e:
                self.logger.error(f"Failed to cache result for {file_path}: {e}")
                return False
    
    def invalidate(self, file_path: str = None) -> int:
        """
        Invalidate cache entries.
        
        Args:
            file_path: Specific file to invalidate, or None to clear all
            
        Returns:
            Number of cache entries removed
        """
        with self.lock:
            removed_count = 0
            
            if file_path is None:
                # Clear entire cache
                for cache_file in self.cache_dir.glob("*.cache"):
                    cache_file.unlink()
                    removed_count += 1
                self.logger.info(f"Cleared entire cache: {removed_count} entries")
            else:
                # Invalidate specific file (all validation types)
                for validation_type in ["ast", "module", "style", "structured"]:
                    cache_key = self._get_cache_key(file_path, validation_type)
                    cache_file = self._get_cache_file_path(cache_key)
                    if cache_file.exists():
                        cache_file.unlink()
                        removed_count += 1
                
                self.logger.debug(f"Invalidated cache for {file_path}: {removed_count} entries")
            
            return removed_count
    
    def _cleanup_if_needed(self) -> None:
        """Clean up cache if it exceeds size limits."""
        try:
            total_size = sum(f.stat().st_size for f in self.cache_dir.glob("*.cache"))
            
            if total_size <= self.max_size_bytes:
                return
            
            # Remove oldest cache files until under size limit
            cache_files = list(self.cache_dir.glob("*.cache"))
            cache_files.sort(key=lambda f: f.stat().st_mtime)  # Oldest first
            
            removed_count = 0
            for cache_file in cache_files:
                cache_file.unlink()
                removed_count += 1
                total_size -= cache_file.stat().st_size
                
                if total_size <= self.max_size_bytes * 0.8:  # Remove extra 20% for buffer
                    break
            
            self.logger.info(f"Cache cleanup: removed {removed_count} old entries")
            
        except Exception as e:
            self.logger.error(f"Cache cleanup failed: {e}")
    
    def get_stats(self) -> Dict[str, Any]:
        """Get cache performance statistics."""
        total_requests = self.hits + self.misses
        hit_rate = (self.hits / total_requests * 100) if total_requests > 0 else 0
        
        cache_files = list(self.cache_dir.glob("*.cache"))
        total_size = sum(f.stat().st_size for f in cache_files)
        
        return {
            "hits": self.hits,
            "misses": self.misses,
            "hit_rate_percent": round(hit_rate, 1),
            "total_entries": len(cache_files),
            "total_size_mb": round(total_size / (1024 * 1024), 2),
            "max_size_mb": round(self.max_size_bytes / (1024 * 1024), 2)
        }


class GitIntegration:
    """
    Git integration for incremental validation and change detection.
    
    Provides functionality to detect changed files, staged files, and compare
    against different branches for efficient incremental validation.
    """
    
    def __init__(self):
        self.logger = logging.getLogger(__name__ + '.GitIntegration')
        self._git_available = self._check_git_available()
    
    def _check_git_available(self) -> bool:
        """Check if git is available and we're in a git repository."""
        try:
            result = subprocess.run(['git', 'rev-parse', '--git-dir'], 
                                  capture_output=True, text=True, timeout=5)
            return result.returncode == 0
        except Exception:
            return False
    
    def get_changed_files(self, base_branch: str = 'main', include_untracked: bool = False) -> List[str]:
        """
        Get list of Python files changed since base branch.
        
        Args:
            base_branch: Base branch to compare against
            include_untracked: Whether to include untracked files
            
        Returns:
            List of changed Python file paths
        """
        if not self._git_available:
            self.logger.warning("Git not available, cannot detect changed files")
            return []
        
        changed_files = []
        
        try:
            # Get files changed compared to base branch
            result = subprocess.run(['git', 'diff', '--name-only', f'{base_branch}...HEAD'], 
                                  capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                changed_files.extend(result.stdout.strip().split('\n'))
            
            # Get unstaged changes
            result = subprocess.run(['git', 'diff', '--name-only'], 
                                  capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                changed_files.extend(result.stdout.strip().split('\n'))
            
            # Get staged changes
            result = subprocess.run(['git', 'diff', '--cached', '--name-only'], 
                                  capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                changed_files.extend(result.stdout.strip().split('\n'))
            
            # Get untracked files if requested
            if include_untracked:
                result = subprocess.run(['git', 'ls-files', '--others', '--exclude-standard'], 
                                      capture_output=True, text=True, timeout=10)
                
                if result.returncode == 0:
                    changed_files.extend(result.stdout.strip().split('\n'))
            
            # Filter to Python files and remove duplicates
            python_files = list(set([
                f for f in changed_files 
                if f.endswith('.py') and f.strip() and Path(f).exists()
            ]))
            
            self.logger.info(f"Found {len(python_files)} changed Python files")
            return python_files
            
        except Exception as e:
            self.logger.error(f"Failed to get changed files: {e}")
            return []
    
    def get_staged_files(self) -> List[str]:
        """
        Get list of staged Python files for pre-commit validation.
        
        Returns:
            List of staged Python file paths
        """
        if not self._git_available:
            return []
        
        try:
            result = subprocess.run(['git', 'diff', '--cached', '--name-only', '--diff-filter=ACMR'], 
                                  capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                staged_files = [
                    f for f in result.stdout.strip().split('\n')
                    if f.endswith('.py') and f.strip() and Path(f).exists()
                ]
                
                self.logger.info(f"Found {len(staged_files)} staged Python files")
                return staged_files
            
        except Exception as e:
            self.logger.error(f"Failed to get staged files: {e}")
        
        return []
    
    def is_file_ignored(self, file_path: str) -> bool:
        """
        Check if a file is ignored by git.
        
        Args:
            file_path: Path to check
            
        Returns:
            True if file is ignored by git
        """
        if not self._git_available:
            return False
        
        try:
            result = subprocess.run(['git', 'check-ignore', file_path], 
                                  capture_output=True, text=True, timeout=5)
            return result.returncode == 0
        except Exception:
            return False


class ValidationConfig:
    """
    Configuration management with multiple source support.
    
    Supports loading configuration from:
    - .cli-validation.yml
    - pyproject.toml 
    - Environment variables
    - Direct initialization
    """
    
    def __init__(self, config_file: Optional[str] = None):
        self.logger = logging.getLogger(__name__ + '.ValidationConfig')
        
        # Default configuration
        self.cache_enabled = True
        self.cache_ttl_hours = 24
        self.cache_max_size_mb = 100
        self.cache_dir = ".validation-cache"
        
        # Performance settings
        self.parallel_workers = 4
        self.incremental_mode = True
        self.timeout_seconds = 300
        
        # Validation rules
        self.coverage_threshold = 95.0
        self.docstring_formats = ['google', 'numpy', 'sphinx']
        self.require_examples = True
        self.validate_parameter_types = True
        self.fail_on_warnings = False
        
        # CI/CD settings
        self.github_annotations = True
        self.pre_commit_quick_mode = True
        self.staged_files_only = True
        self.json_output = False
        
        # Load configuration from various sources
        self._load_configuration(config_file)
    
    def _load_configuration(self, config_file: Optional[str] = None) -> None:
        """Load configuration from multiple sources with priority order."""
        
        # 1. Try to load from specified config file
        if config_file and Path(config_file).exists():
            self._load_from_yaml(config_file)
            return
        
        # 2. Try default config file locations
        default_configs = [
            '.cli-validation.yml',
            '.cli-validation.yaml',
            'cli-validation.yml',
            'cli-validation.yaml'
        ]
        
        for config_path in default_configs:
            if Path(config_path).exists():
                self._load_from_yaml(config_path)
                break
        
        # 3. Try pyproject.toml
        self._load_from_pyproject()
        
        # 4. Load environment variables
        self._load_from_env()
        
        self.logger.debug("Configuration loaded from multiple sources")
    
    def _load_from_yaml(self, config_path: str) -> None:
        """Load configuration from YAML file."""
        try:
            import yaml
            
            with open(config_path) as f:
                config_data = yaml.safe_load(f)
            
            if not config_data:
                return
            
            # Update configuration from YAML
            validation_config = config_data.get('validation', {})
            performance_config = config_data.get('performance', {})
            rules_config = config_data.get('rules', {})
            ci_cd_config = config_data.get('ci_cd', {})
            
            # Cache settings
            cache_config = validation_config.get('cache', {})
            if cache_config:
                self.cache_enabled = cache_config.get('enabled', self.cache_enabled)
                self.cache_ttl_hours = cache_config.get('ttl_hours', self.cache_ttl_hours)
                self.cache_max_size_mb = cache_config.get('max_size_mb', self.cache_max_size_mb)
                self.cache_dir = cache_config.get('dir', self.cache_dir)
            
            # Performance settings
            if performance_config:
                self.parallel_workers = performance_config.get('parallel_workers', self.parallel_workers)
                self.incremental_mode = performance_config.get('incremental', self.incremental_mode)
                self.timeout_seconds = performance_config.get('timeout_seconds', self.timeout_seconds)
            
            # Validation rules
            if rules_config:
                self.coverage_threshold = rules_config.get('coverage_threshold', self.coverage_threshold)
                self.docstring_formats = rules_config.get('docstring_formats', self.docstring_formats)
                self.require_examples = rules_config.get('require_examples', self.require_examples)
                self.validate_parameter_types = rules_config.get('validate_parameter_types', self.validate_parameter_types)
                self.fail_on_warnings = rules_config.get('fail_on_warnings', self.fail_on_warnings)
            
            # CI/CD settings
            if ci_cd_config:
                github_config = ci_cd_config.get('github_actions', {})
                if github_config:
                    self.github_annotations = github_config.get('annotations', self.github_annotations)
                
                precommit_config = ci_cd_config.get('pre_commit', {})
                if precommit_config:
                    self.pre_commit_quick_mode = precommit_config.get('quick_mode', self.pre_commit_quick_mode)
                    self.staged_files_only = precommit_config.get('staged_files_only', self.staged_files_only)
            
            self.logger.info(f"Configuration loaded from {config_path}")
            
        except Exception as e:
            self.logger.warning(f"Failed to load YAML config from {config_path}: {e}")
    
    def _load_from_pyproject(self) -> None:
        """Load configuration from pyproject.toml."""
        if not TOML_SUPPORT:
            self.logger.debug("TOML support not available (tomllib/tomli not found)")
            return
            
        try:
            pyproject_path = Path('pyproject.toml')
            if not pyproject_path.exists():
                return
            
            with open(pyproject_path, 'rb') as f:
                pyproject_data = TOML_LOADER.load(f)
            
            # Look for tool.cli-validation section
            tool_config = pyproject_data.get('tool', {})
            cli_validation_config = tool_config.get('cli-validation', {})
            
            if cli_validation_config:
                # Apply same logic as YAML loading
                self._apply_config_dict(cli_validation_config)
                self.logger.info("Configuration loaded from pyproject.toml")
                
        except Exception as e:
            self.logger.warning(f"Failed to load pyproject.toml config: {e}")
    
    def _load_from_env(self) -> None:
        """Load configuration from environment variables."""
        env_prefix = 'CLI_VALIDATION_'
        
        env_mappings = {
            'CACHE_ENABLED': ('cache_enabled', bool),
            'CACHE_TTL_HOURS': ('cache_ttl_hours', int),
            'CACHE_MAX_SIZE_MB': ('cache_max_size_mb', int),
            'CACHE_DIR': ('cache_dir', str),
            'PARALLEL_WORKERS': ('parallel_workers', int),
            'INCREMENTAL_MODE': ('incremental_mode', bool),
            'TIMEOUT_SECONDS': ('timeout_seconds', int),
            'COVERAGE_THRESHOLD': ('coverage_threshold', float),
            'REQUIRE_EXAMPLES': ('require_examples', bool),
            'VALIDATE_PARAMETER_TYPES': ('validate_parameter_types', bool),
            'FAIL_ON_WARNINGS': ('fail_on_warnings', bool),
            'GITHUB_ANNOTATIONS': ('github_annotations', bool),
            'JSON_OUTPUT': ('json_output', bool)
        }
        
        for env_key, (attr_name, type_func) in env_mappings.items():
            env_value = os.environ.get(env_prefix + env_key)
            if env_value is not None:
                try:
                    if type_func == bool:
                        converted_value = env_value.lower() in ('true', '1', 'yes', 'on')
                    else:
                        converted_value = type_func(env_value)
                    setattr(self, attr_name, converted_value)
                    self.logger.debug(f"Environment variable {env_prefix + env_key} set {attr_name} = {converted_value}")
                except ValueError as e:
                    self.logger.warning(f"Invalid environment variable {env_prefix + env_key}: {e}")
    
    def _apply_config_dict(self, config_dict: Dict[str, Any]) -> None:
        """Apply configuration from a dictionary (helper for multiple sources)."""
        # This method can be used to apply configuration from any dictionary source
        # Implementation similar to _load_from_yaml but more generic
        pass
    
    def get_config_summary(self) -> Dict[str, Any]:
        """Get a summary of current configuration."""
        return {
            'cache': {
                'enabled': self.cache_enabled,
                'ttl_hours': self.cache_ttl_hours,
                'max_size_mb': self.cache_max_size_mb,
                'dir': self.cache_dir
            },
            'performance': {
                'parallel_workers': self.parallel_workers,
                'incremental_mode': self.incremental_mode,
                'timeout_seconds': self.timeout_seconds
            },
            'rules': {
                'coverage_threshold': self.coverage_threshold,
                'docstring_formats': self.docstring_formats,
                'require_examples': self.require_examples,
                'validate_parameter_types': self.validate_parameter_types,
                'fail_on_warnings': self.fail_on_warnings
            },
            'ci_cd': {
                'github_annotations': self.github_annotations,
                'pre_commit_quick_mode': self.pre_commit_quick_mode,
                'staged_files_only': self.staged_files_only,
                'json_output': self.json_output
            }
        }


class PerformanceOptimizer:
    """
    Performance optimization features including parallel processing and progress reporting.
    
    Provides multi-threaded validation, memory optimization, and progress tracking
    for efficient processing of large codebases.
    """
    
    def __init__(self, config: ValidationConfig):
        self.config = config
        self.logger = logging.getLogger(__name__ + '.PerformanceOptimizer')
        
    def validate_files_parallel(self, file_paths: List[str], validation_func: Callable[[str], Any]) -> List[Any]:
        """
        Validate multiple files in parallel using thread pool.
        
        Args:
            file_paths: List of file paths to validate
            validation_func: Function to call for each file
            
        Returns:
            List of validation results
        """
        if len(file_paths) <= 1 or self.config.parallel_workers <= 1:
            # Use sequential processing for single files or when parallel is disabled
            return [validation_func(file_path) for file_path in file_paths]
        
        results = []
        
        with ThreadPoolExecutor(max_workers=self.config.parallel_workers) as executor:
            # Submit all tasks
            future_to_file = {
                executor.submit(validation_func, file_path): file_path 
                for file_path in file_paths
            }
            
            # Collect results with progress reporting
            completed = 0
            total = len(file_paths)
            
            for future in future_to_file:
                try:
                    result = future.result(timeout=self.config.timeout_seconds)
                    results.append(result)
                    completed += 1
                    
                    if completed % 10 == 0 or completed == total:
                        self.logger.info(f"Validation progress: {completed}/{total} files completed")
                        
                except Exception as e:
                    file_path = future_to_file[future]
                    self.logger.error(f"Validation failed for {file_path}: {e}")
                    # Add empty result to maintain order
                    results.append(None)
        
        return results
    
    def optimize_memory_usage(self) -> None:
        """Optimize memory usage by cleaning up caches and forcing garbage collection."""
        import gc
        
        # Force garbage collection
        collected = gc.collect()
        
        self.logger.debug(f"Memory optimization: collected {collected} objects")


class ASTAnalyzer:
    """Advanced AST analyzer for comprehensive code analysis."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__ + '.ASTAnalyzer')
    
    def extract_all_definitions(self, file_path: str) -> List[CodeDefinition]:
        """Extract all functions, methods, classes from a Python file."""
        definitions = []
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            tree = ast.parse(content, filename=file_path)
            
            # Extract classes and their methods
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef):
                    # Add class definition
                    class_docstring = ast.get_docstring(node)
                    definitions.append(CodeDefinition(
                        name=node.name,
                        type='class',
                        file_path=file_path,
                        line_number=node.lineno,
                        docstring=class_docstring
                    ))
                    
                    # Add methods within the class
                    for item in node.body:
                        if isinstance(item, ast.FunctionDef):
                            method_docstring = ast.get_docstring(item)
                            definitions.append(CodeDefinition(
                                name=item.name,
                                type='method',
                                file_path=file_path,
                                line_number=item.lineno,
                                docstring=method_docstring,
                                parent_class=node.name,
                                args=[arg.arg for arg in item.args.args],
                                is_cli_handler=self._is_cli_handler(item, content)
                            ))
                
                elif isinstance(node, ast.FunctionDef) and not self._is_nested_function(node, tree):
                    # Add top-level function definitions
                    func_docstring = ast.get_docstring(node)
                    definitions.append(CodeDefinition(
                        name=node.name,
                        type='function',
                        file_path=file_path,
                        line_number=node.lineno,
                        docstring=func_docstring,
                        args=[arg.arg for arg in node.args.args],
                        is_cli_handler=self._is_cli_handler(node, content)
                    ))
        
        except Exception as e:
            self.logger.error(f"Failed to parse {file_path}: {e}")
        
        return definitions
    
    def _is_cli_handler(self, node: ast.FunctionDef, content: str) -> bool:
        """Detect if a function is likely a CLI command handler."""
        # Check for common CLI patterns
        cli_indicators = [
            'subparsers.add_parser',
            'add_parser',
            'set_defaults',
            'argparse',
            'parser.add_argument'
        ]
        
        # Get function source
        try:
            lines = content.split('\n')
            func_lines = lines[node.lineno-1:node.end_lineno if hasattr(node, 'end_lineno') else node.lineno+20]
            func_source = '\n'.join(func_lines)
            
            return any(indicator in func_source for indicator in cli_indicators)
        except:
            return False
    
    def _is_nested_function(self, node: ast.FunctionDef, tree: ast.AST) -> bool:
        """Check if a function is nested inside another function."""
        for parent in ast.walk(tree):
            if isinstance(parent, (ast.FunctionDef, ast.AsyncFunctionDef)) and parent != node:
                for child in ast.walk(parent):
                    if child == node:
                        return True
        return False
    
    def analyze_docstring_coverage(self, definitions: List[CodeDefinition]) -> CoverageReport:
        """Calculate comprehensive docstring coverage metrics."""
        total = len(definitions)
        documented = sum(1 for d in definitions if d.docstring and d.docstring.strip())
        undocumented = total - documented
        
        coverage_percentage = (documented / total * 100) if total > 0 else 100.0
        
        undocumented_items = [d for d in definitions if not d.docstring or not d.docstring.strip()]
        
        return CoverageReport(
            total_definitions=total,
            documented_definitions=documented,
            undocumented_definitions=undocumented,
            coverage_percentage=coverage_percentage,
            undocumented_items=undocumented_items
        )
    
    def extract_function_signature(self, node: ast.FunctionDef) -> FunctionSignature:
        """
        Extract complete function signature including parameter and return types.
        
        Args:
            node: AST FunctionDef node
            
        Returns:
            FunctionSignature: Complete signature information
        """
        parameters = []
        
        # Extract parameters with type annotations
        for arg in node.args.args:
            param_type = None
            default_value = None
            
            # Get type annotation if present
            if arg.annotation:
                try:
                    param_type = ast.unparse(arg.annotation)
                except AttributeError:
                    # For Python < 3.9, use a simpler approach
                    param_type = self._annotation_to_string(arg.annotation)
            
            # Check for default values
            num_defaults = len(node.args.defaults)
            num_params = len(node.args.args)
            default_offset = num_params - num_defaults
            current_param_index = node.args.args.index(arg)
            
            if current_param_index >= default_offset:
                default_index = current_param_index - default_offset
                if default_index < len(node.args.defaults):
                    try:
                        default_value = ast.unparse(node.args.defaults[default_index])
                    except AttributeError:
                        default_value = self._node_to_string(node.args.defaults[default_index])
            
            parameters.append(Parameter(
                name=arg.arg,
                type_hint=param_type,
                is_optional=default_value is not None,
                default_value=default_value
            ))
        
        # Extract return type annotation
        return_type = None
        if node.returns:
            try:
                return_type = ast.unparse(node.returns)
            except AttributeError:
                return_type = self._annotation_to_string(node.returns)
        
        return FunctionSignature(parameters=parameters, return_type=return_type)
    
    def _annotation_to_string(self, annotation) -> str:
        """Convert AST annotation to string (fallback for older Python versions)."""
        if isinstance(annotation, ast.Name):
            return annotation.id
        elif isinstance(annotation, ast.Constant):
            return str(annotation.value)
        elif isinstance(annotation, ast.Attribute):
            return f"{self._annotation_to_string(annotation.value)}.{annotation.attr}"
        elif isinstance(annotation, ast.Subscript):
            value = self._annotation_to_string(annotation.value)
            slice_val = self._annotation_to_string(annotation.slice)
            return f"{value}[{slice_val}]"
        else:
            return str(annotation)
    
    def _node_to_string(self, node) -> str:
        """Convert AST node to string representation."""
        if isinstance(node, ast.Constant):
            return repr(node.value)
        elif isinstance(node, ast.Name):
            return node.id
        elif isinstance(node, ast.Attribute):
            return f"{self._node_to_string(node.value)}.{node.attr}"
        else:
            return str(node)


class MultiFormatDocstringParser:
    """
    Multi-format docstring parser with structured validation capabilities.
    
    Supports Google, NumPy, ReST, and Epydoc docstring formats using the
    docstring_parser library with graceful fallback for unsupported formats.
    """
    
    def __init__(self):
        self.logger = logging.getLogger(__name__ + '.MultiFormatDocstringParser')
        self.supported_formats = ['google', 'numpy', 'sphinx', 'epydoc', 'auto']
        
        if not DOCSTRING_PARSER_AVAILABLE:
            self.logger.warning("docstring_parser library not available. Structured parsing disabled.")
    
    def parse_docstring(self, docstring_text: str, format_hint: str = 'auto') -> ParsedDocstring:
        """
        Parse a docstring into structured components.
        
        Args:
            docstring_text: Raw docstring text to parse
            format_hint: Format hint ('google', 'numpy', 'sphinx', 'epydoc', 'auto')
            
        Returns:
            ParsedDocstring: Structured representation of the docstring
        """
        if not docstring_text or not docstring_text.strip():
            return ParsedDocstring(
                raw_text="",
                detected_format="none",
                summary="",
                description="",
                parsing_successful=False,
                parsing_errors=["Empty docstring"]
            )
        
        # If docstring_parser is not available, return basic parsing
        if not DOCSTRING_PARSER_AVAILABLE:
            return self._basic_parse_docstring(docstring_text)
        
        try:
            # Parse using docstring_parser
            parsed = parse_docstring(docstring_text)
            
            # Detect format
            detected_format = self._detect_format(docstring_text)
            
            # Extract components
            parameters = []
            for param in parsed.params:
                parameters.append(Parameter(
                    name=param.arg_name,
                    type_hint=param.type_name,
                    description=param.description,
                    is_optional=param.is_optional
                ))
            
            # Extract return information
            returns = None
            if parsed.returns:
                returns = ReturnInfo(
                    type_hint=parsed.returns.type_name,
                    description=parsed.returns.description
                )
            
            # Extract exception information
            raises = []
            for exc in parsed.raises:
                raises.append(ExceptionInfo(
                    type_name=exc.type_name,
                    description=exc.description
                ))
            
            # Extract examples
            examples = []
            if parsed.examples:
                for example in parsed.examples:
                    examples.append(example.description)
            
            return ParsedDocstring(
                raw_text=docstring_text,
                detected_format=detected_format,
                summary=parsed.short_description or "",
                description=parsed.long_description or "",
                parameters=parameters,
                returns=returns,
                raises=raises,
                examples=examples,
                parsing_successful=True
            )
        
        except Exception as e:
            self.logger.debug(f"Failed to parse docstring with docstring_parser: {e}")
            # Fallback to basic parsing
            basic_result = self._basic_parse_docstring(docstring_text)
            basic_result.parsing_errors.append(f"Structured parsing failed: {e}")
            return basic_result
    
    def _detect_format(self, docstring_text: str) -> str:
        """
        Detect the docstring format based on content patterns.
        
        Args:
            docstring_text: Raw docstring text
            
        Returns:
            str: Detected format ('google', 'numpy', 'sphinx', 'unknown')
        """
        import re
        
        # Google style patterns
        google_patterns = [
            r'^\s*Args:\s*$',
            r'^\s*Arguments:\s*$',
            r'^\s*Returns:\s*$',
            r'^\s*Yields:\s*$',
            r'^\s*Raises:\s*$'
        ]
        
        # NumPy style patterns  
        numpy_patterns = [
            r'^\s*Parameters\s*\n\s*-+',
            r'^\s*Returns\s*\n\s*-+',
            r'^\s*Yields\s*\n\s*-+',
            r'^\s*Raises\s*\n\s*-+'
        ]
        
        # Sphinx/ReST style patterns
        sphinx_patterns = [
            r':param\s+\w+:',
            r':type\s+\w+:',
            r':returns?:',
            r':rtype:',
            r':raises?\s+\w+:'
        ]
        
        # Check for each format
        for pattern in google_patterns:
            if re.search(pattern, docstring_text, re.MULTILINE):
                return 'google'
        
        for pattern in numpy_patterns:
            if re.search(pattern, docstring_text, re.MULTILINE):
                return 'numpy'
        
        for pattern in sphinx_patterns:
            if re.search(pattern, docstring_text, re.MULTILINE):
                return 'sphinx'
        
        return 'unknown'
    
    def _basic_parse_docstring(self, docstring_text: str) -> ParsedDocstring:
        """
        Basic docstring parsing fallback when structured parsing fails.
        
        Args:
            docstring_text: Raw docstring text
            
        Returns:
            ParsedDocstring: Basic parsed representation
        """
        lines = docstring_text.strip().split('\n')
        summary = lines[0].strip() if lines else ""
        
        # Find description (lines after summary but before sections)
        description_lines = []
        for i, line in enumerate(lines[1:], 1):
            if line.strip() and not line.startswith(' '):
                # Potential section header
                if any(section in line.lower() for section in ['args:', 'parameters', 'returns:', 'raises:']):
                    break
            description_lines.append(line)
        
        description = '\n'.join(description_lines).strip()
        detected_format = self._detect_format(docstring_text)
        
        return ParsedDocstring(
            raw_text=docstring_text,
            detected_format=detected_format,
            summary=summary,
            description=description,
            parsing_successful=False,
            parsing_errors=["Using basic parsing fallback"]
        )
    
    def validate_structure(self, parsed_docstring: ParsedDocstring, function_signature: FunctionSignature) -> StructureValidation:
        """
        Validate parsed docstring structure against function signature.
        
        Args:
            parsed_docstring: Parsed docstring components
            function_signature: Function signature from AST analysis
            
        Returns:
            StructureValidation: Validation results with suggestions
        """
        missing_parameters = []
        extra_parameters = []
        type_mismatches = []
        format_issues = []
        suggestions = []
        
        # Get parameter names from docstring and signature
        docstring_params = {p.name for p in parsed_docstring.parameters}
        signature_params = {p.name for p in function_signature.parameters}
        
        # Check for missing parameters in docstring
        missing_parameters = list(signature_params - docstring_params)
        
        # Check for extra parameters in docstring
        extra_parameters = list(docstring_params - signature_params)
        
        # Check return documentation
        missing_return_doc = (function_signature.return_type is not None and 
                            function_signature.return_type != 'None' and 
                            parsed_docstring.returns is None)
        
        # Generate suggestions based on validation results
        if missing_parameters:
            suggestions.append(f"Document missing parameters: {', '.join(missing_parameters)}")
        
        if extra_parameters:
            suggestions.append(f"Remove or correct extra parameters: {', '.join(extra_parameters)}")
        
        if missing_return_doc:
            suggestions.append("Add return value documentation")
        
        if not parsed_docstring.parsing_successful:
            suggestions.append("Fix docstring format for better structured parsing")
        
        # Check format consistency
        if parsed_docstring.detected_format == 'unknown':
            format_issues.append("Docstring format not recognized")
            suggestions.append("Use consistent docstring format (Google, NumPy, or Sphinx)")
        
        is_valid = (not missing_parameters and 
                   not extra_parameters and 
                   not missing_return_doc and
                   not format_issues)
        
        return StructureValidation(
            is_valid=is_valid,
            missing_parameters=missing_parameters,
            extra_parameters=extra_parameters,
            missing_return_doc=missing_return_doc,
            type_mismatches=type_mismatches,
            format_issues=format_issues,
            suggestions=suggestions
        )


class EnhancedCLIValidator:
    """Enhanced CLI documentation validation engine with AST analysis, multi-format parsing, and CI/CD integration."""
    
    def __init__(self, config_file: Optional[str] = None):
        """Initialize the enhanced validator with all components."""
        # Load configuration
        self.config = ValidationConfig(config_file)
        
        # Initialize core components
        self.ast_analyzer = ASTAnalyzer()
        self.docstring_parser = MultiFormatDocstringParser()
        self.cache_manager = CacheManager(
            cache_dir=self.config.cache_dir,
            ttl_hours=self.config.cache_ttl_hours,
            max_size_mb=self.config.cache_max_size_mb
        ) if self.config.cache_enabled else None
        self.git_integration = GitIntegration()
        self.performance_optimizer = PerformanceOptimizer(self.config)
        
        # Initialize existing tools
        self.module_validator = CLIDocumentationValidator()
        self.style_validator = CLIDocstringValidator()
        
        self.logger = logging.getLogger(__name__)
    
    def validate_comprehensive(self, target_paths: List[str], mode: str = "comprehensive") -> ValidationReport:
        """
        Complete validation using all available tools with caching and performance optimization.
        
        Args:
            target_paths: List of paths to validate (files or directories)
            mode: Validation mode (quick, comprehensive, strict)
            
        Returns:
            ValidationReport with comprehensive results
        """
        start_time = time.time()
        
        self.logger.info(f"Starting comprehensive validation in {mode} mode")
        
        # Collect all Python files
        all_files = []
        for target_path in target_paths:
            if Path(target_path).is_file():
                all_files.append(target_path)
            else:
                all_files.extend(self._collect_python_files(target_path))
        
        # Remove duplicates and filter out ignored files
        unique_files = list(set(all_files))
        filtered_files = [f for f in unique_files if not self.git_integration.is_file_ignored(f)]
        
        self.logger.info(f"Found {len(filtered_files)} Python files to validate")
        
        # Check for incremental validation opportunity
        if self.config.incremental_mode and mode == "quick":
            changed_files = self.git_integration.get_changed_files()
            if changed_files:
                # Only validate changed files that are in our target set
                filtered_files = [f for f in filtered_files if f in changed_files]
                self.logger.info(f"Incremental mode: validating {len(filtered_files)} changed files")
        
        # Perform validation with parallel processing and caching
        validation_results = self._validate_files_with_cache(filtered_files, mode)
        
        # Aggregate results
        report = self._aggregate_validation_results(validation_results, start_time)
        
        # Add cache statistics if available
        if self.cache_manager:
            cache_stats = self.cache_manager.get_stats()
            report.metadata["cache_stats"] = cache_stats
            self.logger.info(f"Cache performance: {cache_stats['hit_rate_percent']}% hit rate")
        
        return report
    
    def validate_incremental(self, base_branch: str = "main") -> ValidationReport:
        """
        Performance-optimized validation for CI/CD workflows using Git integration.
        
        Args:
            base_branch: Base branch to compare against for changed files
            
        Returns:
            ValidationReport with results for changed files only
        """
        start_time = time.time()
        
        self.logger.info(f"Starting incremental validation against {base_branch}")
        
        # Get changed files
        changed_files = self.git_integration.get_changed_files(base_branch, include_untracked=True)
        
        if not changed_files:
            self.logger.info("No changed Python files found")
            return ValidationReport(
                total_files=0,
                total_definitions=0,
                coverage_percentage=100.0,
                execution_time=time.time() - start_time,
                issues=[],
                metadata={"incremental": True, "base_branch": base_branch}
            )
        
        # Validate changed files with full caching support
        validation_results = self._validate_files_with_cache(changed_files, "comprehensive")
        
        # Aggregate results
        report = self._aggregate_validation_results(validation_results, start_time)
        report.metadata["incremental"] = True
        report.metadata["base_branch"] = base_branch
        report.metadata["changed_files"] = changed_files
        
        return report
    
    def validate_staged_files(self) -> ValidationReport:
        """
        Validate only staged files for pre-commit hook integration.
        
        Returns:
            ValidationReport with results for staged files only
        """
        start_time = time.time()
        
        self.logger.info("Starting staged files validation for pre-commit")
        
        # Get staged files
        staged_files = self.git_integration.get_staged_files()
        
        if not staged_files:
            self.logger.info("No staged Python files found")
            return ValidationReport(
                total_files=0,
                total_definitions=0,
                coverage_percentage=100.0,
                execution_time=time.time() - start_time,
                issues=[],
                metadata={"pre_commit": True}
            )
        
        # Use quick mode for pre-commit to maintain fast feedback
        mode = "quick" if self.config.pre_commit_quick_mode else "comprehensive"
        validation_results = self._validate_files_with_cache(staged_files, mode)
        
        # Aggregate results
        report = self._aggregate_validation_results(validation_results, start_time)
        report.metadata["pre_commit"] = True
        report.metadata["staged_files"] = staged_files
        
        return report
    
    def _validate_files_with_cache(self, file_paths: List[str], mode: str) -> List[FileValidationResult]:
        """
        Validate files with caching support and parallel processing.
        
        Args:
            file_paths: List of file paths to validate
            mode: Validation mode (quick, comprehensive, strict)
            
        Returns:
            List of FileValidationResult objects
        """
        def validate_single_file(file_path: str) -> FileValidationResult:
            """Validate a single file with caching."""
            
            # Try cache first
            if self.cache_manager:
                cached_result = self.cache_manager.get(file_path, f"validation_{mode}")
                if cached_result:
                    return cached_result
            
            # Perform validation
            try:
                result = self._validate_file_comprehensive(file_path, mode)
                
                # Cache the result
                if self.cache_manager:
                    self.cache_manager.set(file_path, result, f"validation_{mode}")
                
                return result
                
            except Exception as e:
                self.logger.error(f"Validation failed for {file_path}: {e}")
                return FileValidationResult(
                    file_path=file_path,
                    definitions=[],
                    issues=[ValidationIssue(
                        type="error",
                        severity="error",
                        message=f"Validation failed: {str(e)}",
                        location=ValidationLocation(file_path=file_path, line_number=1),
                        suggestions=[]
                    )],
                    total_definitions=0,
                    documented_count=0,
                    coverage_percentage=0.0,
                    execution_time=0.0,
                    metadata={"error": True}
                )
        
        # Use parallel processing if enabled
        if self.config.parallel_workers > 1 and len(file_paths) > 1:
            return self.performance_optimizer.validate_files_parallel(file_paths, validate_single_file)
        else:
            return [validate_single_file(file_path) for file_path in file_paths]
    
    def generate_ci_report(self, results: ValidationReport) -> Dict:
        """
        Generate JSON output for CI/CD integration with GitHub Actions annotations.
        
        Args:
            results: ValidationReport to convert
            
        Returns:
            Dictionary suitable for JSON serialization
        """
        # Convert ValidationReport to serializable format
        ci_report = {
            "validation_results": {
                "total_files": results.total_files,
                "total_definitions": results.total_definitions,
                "coverage_percentage": results.coverage_percentage,
                "execution_time": results.execution_time,
                "success": results.coverage_percentage >= self.config.coverage_threshold
            },
            "issues": [
                {
                    "type": issue.type,
                    "severity": issue.severity,
                    "message": issue.message,
                    "file": issue.location.file_path,
                    "line": issue.location.line_number,
                    "suggestions": issue.suggestions
                }
                for issue in results.issues
            ],
            "metadata": results.metadata
        }
        
        # Add GitHub Actions annotations if enabled
        if self.config.github_annotations:
            annotations = []
            for issue in results.issues:
                annotation_level = "error" if issue.severity == "error" else "warning"
                annotations.append({
                    "path": issue.location.file_path,
                    "start_line": issue.location.line_number,
                    "end_line": issue.location.line_number,
                    "annotation_level": annotation_level,
                    "message": issue.message,
                    "title": f"CLI Documentation {issue.type}"
                })
            
            ci_report["github_annotations"] = annotations
        
        return ci_report
    
    def clear_cache(self) -> Dict[str, Any]:
        """
        Clear validation cache and return statistics.
        
        Returns:
            Dictionary with cache clearing statistics
        """
        if not self.cache_manager:
            return {"cache_enabled": False, "entries_removed": 0}
        
        entries_removed = self.cache_manager.invalidate()
        return {
            "cache_enabled": True,
            "entries_removed": entries_removed,
            "cache_dir": self.config.cache_dir
        }
    
    def get_performance_stats(self) -> Dict[str, Any]:
        """
        Get comprehensive performance statistics.
        
        Returns:
            Dictionary with performance metrics
        """
        stats = {
            "configuration": self.config.get_config_summary(),
            "git_integration": {
                "available": self.git_integration._git_available
            }
        }
        
        if self.cache_manager:
            stats["cache"] = self.cache_manager.get_stats()
        
        return stats
    
    def _collect_python_files(self, target_path: str) -> List[str]:
        """Collect all Python files from a target path."""
        target = Path(target_path)
        
        if target.is_file():
            return [str(target)] if target.suffix == '.py' else []
        elif target.is_dir():
            python_files = []
            for py_file in target.rglob('*.py'):
                if not self.git_integration.is_file_ignored(str(py_file)):
                    python_files.append(str(py_file))
            return python_files
        else:
            return []
    
    def _validate_file_comprehensive(self, file_path: str, mode: str) -> FileValidationResult:
        """
        Perform comprehensive validation on a single file.
        
        Args:
            file_path: Path to the file to validate
            mode: Validation mode (quick, comprehensive, strict)
            
        Returns:
            FileValidationResult with all validation results
        """
        start_time = time.time()
        
        try:
            # Extract all definitions using AST
            definitions = self.ast_analyzer.extract_all_definitions(file_path)
            
            # Calculate coverage
            documented_count = sum(1 for d in definitions if d.docstring and d.docstring.strip())
            total_count = len(definitions)
            coverage = (documented_count / total_count * 100) if total_count > 0 else 100.0
            
            # Collect validation issues
            issues = []
            
            # Check for undocumented definitions
            for definition in definitions:
                if not definition.docstring or not definition.docstring.strip():
                    issues.append(ValidationIssue(
                        type="missing_docstring",
                        severity="warning",
                        message=f"{definition.type.title()} '{definition.name}' is missing docstring",
                        location=ValidationLocation(file_path=file_path, line_number=definition.line_number),
                        suggestions=[f"Add a docstring to document the {definition.type}"]
                    ))
            
            # Perform structured docstring validation if in comprehensive mode
            if mode in ["comprehensive", "strict"] and DOCSTRING_PARSER_AVAILABLE:
                for definition in definitions:
                    if definition.docstring and definition.type in ['function', 'method']:
                        try:
                            parsed = self.docstring_parser.parse_docstring(definition.docstring)
                            signature = self.ast_analyzer.extract_function_signature(definition)
                            
                            if signature:
                                structure_validation = self.docstring_parser.validate_structure(parsed, signature)
                                
                                # Add structured validation issues
                                for missing_param in structure_validation.missing_parameters:
                                    issues.append(ValidationIssue(
                                        type="missing_parameter_doc",
                                        severity="warning",
                                        message=f"Parameter '{missing_param}' not documented",
                                        location=ValidationLocation(file_path=file_path, line_number=definition.line_number),
                                        suggestions=[f"Add '{missing_param}' to docstring parameters"]
                                    ))
                                
                                for extra_param in structure_validation.extra_parameters:
                                    issues.append(ValidationIssue(
                                        type="extra_parameter_doc",
                                        severity="warning", 
                                        message=f"Parameter '{extra_param}' documented but not in signature",
                                        location=ValidationLocation(file_path=file_path, line_number=definition.line_number),
                                        suggestions=[f"Remove '{extra_param}' from docstring or add to function"]
                                    ))
                        except Exception as e:
                            self.logger.debug(f"Structured validation failed for {definition.name}: {e}")
            
            execution_time = time.time() - start_time
            
            return FileValidationResult(
                file_path=file_path,
                definitions=definitions,
                issues=issues,
                total_definitions=total_count,
                documented_count=documented_count,
                coverage_percentage=coverage,
                execution_time=execution_time,
                metadata={"mode": mode}
            )
            
        except Exception as e:
            self.logger.error(f"File validation failed for {file_path}: {e}")
            return FileValidationResult(
                file_path=file_path,
                definitions=[],
                issues=[ValidationIssue(
                    type="validation_error",
                    severity="error",
                    message=f"Validation failed: {str(e)}",
                    location=ValidationLocation(file_path=file_path, line_number=1),
                    suggestions=[]
                )],
                total_definitions=0,
                documented_count=0,
                coverage_percentage=0.0,
                execution_time=time.time() - start_time,
                metadata={"error": True}
            )
    
    def _aggregate_validation_results(self, file_results: List[FileValidationResult], start_time: float) -> ValidationReport:
        """
        Aggregate individual file validation results into a comprehensive report.
        
        Args:
            file_results: List of FileValidationResult objects
            start_time: Start time of the validation process
            
        Returns:
            ValidationReport with aggregated results
        """
        # Filter out None results
        valid_results = [r for r in file_results if r is not None]
        
        # Aggregate statistics
        total_files = len(valid_results)
        total_definitions = sum(r.total_definitions for r in valid_results)
        total_documented = sum(r.documented_count for r in valid_results)
        
        # Calculate overall coverage
        overall_coverage = (total_documented / total_definitions * 100) if total_definitions > 0 else 100.0
        
        # Collect all issues
        all_issues = []
        for result in valid_results:
            all_issues.extend(result.issues)
        
        execution_time = time.time() - start_time
        
        return ValidationReport(
            total_files=total_files,
            total_definitions=total_definitions,
            coverage_percentage=overall_coverage,
            execution_time=execution_time,
            issues=all_issues,
            metadata={}
        )
    
    def _extract_function_signature_for_definition(self, definition: CodeDefinition) -> Optional[FunctionSignature]:
        """
        Extract function signature for a given code definition.
        
        Args:
            definition: Code definition to extract signature for
            
        Returns:
            Optional[FunctionSignature]: Function signature if extractable
        """
        try:
            # Re-parse the file to get the AST node
            with open(definition.file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            tree = ast.parse(content, filename=definition.file_path)
            
            # Find the specific function node
            for node in ast.walk(tree):
                if isinstance(node, ast.FunctionDef):
                    if (node.name == definition.name and 
                        node.lineno == definition.line_number):
                        return self.ast_analyzer.extract_function_signature(node)
            
            return None
        
        except Exception as e:
            self.logger.debug(f"Failed to extract signature for {definition.name}: {e}")
            return None


def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler()
        ]
    )


def main():
    """Main CLI interface for enhanced CLI documentation validation."""
    parser = argparse.ArgumentParser(
        description="Enhanced CLI Documentation Validation Engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic validation
  python scripts/cli_enhanced_validation.py src/vcf_agent/cli/

  # Quick validation with cache
  python scripts/cli_enhanced_validation.py src/vcf_agent/ --mode quick --cache

  # Incremental validation for CI/CD
  python scripts/cli_enhanced_validation.py --incremental --base-branch main

  # Pre-commit hook validation
  python scripts/cli_enhanced_validation.py --staged-files

  # JSON output for CI integration
  python scripts/cli_enhanced_validation.py src/vcf_agent/ --json

  # Clear cache
  python scripts/cli_enhanced_validation.py --clear-cache

  # Performance statistics
  python scripts/cli_enhanced_validation.py --stats
        """
    )
    
    # Positional arguments
    parser.add_argument(
        'paths',
        nargs='*',
        help='File or directory paths to validate (default: current directory)'
    )
    
    # Validation mode options
    parser.add_argument(
        '--mode',
        choices=['quick', 'comprehensive', 'strict'],
        default='comprehensive',
        help='Validation mode (default: comprehensive)'
    )
    
    # Configuration options
    parser.add_argument(
        '--config',
        type=str,
        help='Path to configuration file (.cli-validation.yml)'
    )
    
    parser.add_argument(
        '--no-cache',
        action='store_true',
        help='Disable caching for this run'
    )
    
    # Git integration options
    parser.add_argument(
        '--incremental',
        action='store_true',
        help='Only validate files changed since base branch'
    )
    
    parser.add_argument(
        '--base-branch',
        type=str,
        default='main',
        help='Base branch for incremental validation (default: main)'
    )
    
    parser.add_argument(
        '--staged-files',
        action='store_true',
        help='Only validate staged files (for pre-commit hooks)'
    )
    
    # Output options
    parser.add_argument(
        '--json',
        action='store_true',
        help='Output results in JSON format for CI/CD integration'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--github-annotations',
        action='store_true',
        help='Include GitHub Actions annotations in JSON output'
    )
    
    # Utility options
    parser.add_argument(
        '--clear-cache',
        action='store_true',
        help='Clear validation cache and exit'
    )
    
    parser.add_argument(
        '--stats',
        action='store_true',
        help='Show performance statistics and exit'
    )
    
    # Performance options
    parser.add_argument(
        '--parallel-workers',
        type=int,
        help='Number of parallel workers for validation'
    )
    
    parser.add_argument(
        '--timeout',
        type=int,
        help='Timeout in seconds for validation operations'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    try:
        # Initialize validator with configuration
        validator = EnhancedCLIValidator(config_file=args.config)
        
        # Override configuration with CLI arguments
        if args.no_cache:
            validator.config.cache_enabled = False
            validator.cache_manager = None
        
        if args.parallel_workers:
            validator.config.parallel_workers = args.parallel_workers
        
        if args.timeout:
            validator.config.timeout_seconds = args.timeout
        
        if args.github_annotations:
            validator.config.github_annotations = True
        
        if args.json:
            validator.config.json_output = True
        
        # Handle utility commands
        if args.clear_cache:
            cache_stats = validator.clear_cache()
            if args.json:
                print(json.dumps(cache_stats, indent=2))
            else:
                if cache_stats["cache_enabled"]:
                    print(f"✅ Cache cleared: {cache_stats['entries_removed']} entries removed")
                    print(f"Cache directory: {cache_stats['cache_dir']}")
                else:
                    print("❌ Cache is not enabled")
            return 0
        
        if args.stats:
            perf_stats = validator.get_performance_stats()
            if args.json:
                print(json.dumps(perf_stats, indent=2))
            else:
                print("📊 Performance Statistics:")
                config = perf_stats["configuration"]
                print(f"Cache: {'enabled' if config['cache']['enabled'] else 'disabled'}")
                print(f"Parallel workers: {config['performance']['parallel_workers']}")
                print(f"Git integration: {'available' if perf_stats['git_integration']['available'] else 'unavailable'}")
                
                if "cache" in perf_stats:
                    cache = perf_stats["cache"]
                    print(f"Cache hit rate: {cache['hit_rate_percent']}%")
                    print(f"Cache entries: {cache['total_entries']}")
                    print(f"Cache size: {cache['total_size_mb']} MB")
            return 0
        
        # Determine validation approach
        if args.staged_files:
            # Pre-commit hook validation
            report = validator.validate_staged_files()
        elif args.incremental:
            # Incremental validation for CI/CD
            report = validator.validate_incremental(args.base_branch)
        else:
            # Standard comprehensive validation
            target_paths = args.paths if args.paths else ['.']
            report = validator.validate_comprehensive(target_paths, args.mode)
        
        # Output results
        if args.json:
            ci_report = validator.generate_ci_report(report)
            print(json.dumps(ci_report, indent=2, default=str))
        else:
            # Console output
            print(f"\n{'='*60}")
            print(f"CLI Documentation Validation Report")
            print(f"{'='*60}")
            
            if hasattr(report, 'metadata') and report.metadata.get('incremental'):
                print(f"🔄 Incremental validation against '{report.metadata.get('base_branch', 'main')}'")
                if report.metadata.get('changed_files'):
                    print(f"📁 Changed files: {len(report.metadata['changed_files'])}")
            elif hasattr(report, 'metadata') and report.metadata.get('pre_commit'):
                print(f"🚀 Pre-commit validation")
                if report.metadata.get('staged_files'):
                    print(f"📁 Staged files: {len(report.metadata['staged_files'])}")
            else:
                print(f"📊 Comprehensive validation")
            
            print(f"📁 Files analyzed: {report.total_files}")
            print(f"🔍 Definitions found: {report.total_definitions}")
            print(f"📝 Coverage: {report.coverage_percentage:.1f}%")
            print(f"⏱️  Execution time: {report.execution_time:.2f}s")
            
            # Show cache statistics if available
            if hasattr(report, 'metadata') and 'cache_stats' in report.metadata:
                cache_stats = report.metadata['cache_stats']
                print(f"💾 Cache hit rate: {cache_stats['hit_rate_percent']}%")
            
            # Show issues summary
            if report.issues:
                errors = [i for i in report.issues if i.severity == 'error']
                warnings = [i for i in report.issues if i.severity == 'warning']
                info = [i for i in report.issues if i.severity == 'info']
                
                print(f"\n📋 Issues found:")
                if errors:
                    print(f"   ❌ Errors: {len(errors)}")
                if warnings:
                    print(f"   ⚠️  Warnings: {len(warnings)}")
                if info:
                    print(f"   ℹ️  Info: {len(info)}")
                
                # Show detailed issues in verbose mode
                if args.verbose:
                    print(f"\n{'='*60}")
                    print(f"Detailed Issues:")
                    print(f"{'='*60}")
                    
                    for issue in report.issues[:20]:  # Limit to first 20 issues
                        severity_icon = {'error': '❌', 'warning': '⚠️', 'info': 'ℹ️'}.get(issue.severity, '•')
                        print(f"\n{severity_icon} {issue.type.upper()}")
                        print(f"   📄 {issue.location.file_path}:{issue.location.line_number}")
                        print(f"   💬 {issue.message}")
                        if issue.suggestions:
                            print(f"   💡 {issue.suggestions[0]}")
                    
                    if len(report.issues) > 20:
                        print(f"\n... and {len(report.issues) - 20} more issues")
            else:
                print(f"\n✅ No issues found")
            
            # Show success/failure status
            print(f"\n{'='*60}")
            success = report.coverage_percentage >= validator.config.coverage_threshold
            if success:
                print(f"✅ VALIDATION PASSED")
                if not validator.config.fail_on_warnings or not any(i.severity == 'warning' for i in report.issues):
                    exit_code = 0
                else:
                    print(f"⚠️  (but warnings found and fail_on_warnings=True)")
                    exit_code = 1
            else:
                print(f"❌ VALIDATION FAILED")
                print(f"   Coverage {report.coverage_percentage:.1f}% below threshold {validator.config.coverage_threshold}%")
                exit_code = 1
        
        # Determine exit code
        if args.json:
            # For JSON output, success is determined by configuration
            success = report.coverage_percentage >= validator.config.coverage_threshold
            if not success:
                exit_code = 1
            elif validator.config.fail_on_warnings and any(i.severity == 'warning' for i in report.issues):
                exit_code = 1
            else:
                exit_code = 0
        
        return exit_code
        
    except Exception as e:
        if args.json:
            error_report = {
                "error": True,
                "message": str(e),
                "validation_results": {
                    "success": False,
                    "total_files": 0,
                    "total_definitions": 0,
                    "coverage_percentage": 0.0,
                    "execution_time": 0.0
                }
            }
            print(json.dumps(error_report, indent=2))
        else:
            logging.error(f"Validation failed: {e}")
            print(f"❌ ERROR: {e}")
        
        return 2  # Script execution error


if __name__ == "__main__":
    main() 