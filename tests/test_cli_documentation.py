"""
Test CLI Documentation Completeness

Tests to ensure the CLI module docstring accurately documents all implemented commands.
This prevents documentation drift and ensures the module documentation is trustworthy.
"""

import pytest
import re
import ast
from pathlib import Path
from typing import Set, List
from unittest.mock import patch


@pytest.fixture
def cli_module_path():
    """Path to the CLI module."""
    return Path("src/vcf_agent/cli.py")

@pytest.fixture
def cli_module_content(cli_module_path):
    """Content of the CLI module."""
    if not cli_module_path.exists():
        pytest.skip(f"CLI module not found at {cli_module_path}")
    
    with open(cli_module_path, 'r', encoding='utf-8') as f:
        return f.read()

@pytest.fixture
def cli_module_docstring(cli_module_content):
    """Extract the CLI module docstring."""
    tree = ast.parse(cli_module_content)
    
    if not (tree.body and isinstance(tree.body[0], ast.Expr) 
            and isinstance(tree.body[0].value, ast.Constant)):
        pytest.fail("CLI module has no docstring")
    
    return tree.body[0].value.value


class TestCLIDocumentationCompleteness:
    """Test suite for CLI documentation completeness validation."""
    
    def extract_documented_commands(self, docstring: str) -> Set[str]:
        """Extract command names from the module docstring."""
        commands = set()
        
        # Pattern to match command documentation lines like "• command-name: Description"
        command_pattern = r'^[•*-]\s+([a-z][a-z0-9-]*):?\s+'
        
        for line in docstring.split('\n'):
            line = line.strip()
            match = re.match(command_pattern, line)
            if match:
                command_name = match.group(1)
                # Filter out non-command entries and subcommands
                if not command_name.startswith('--') and command_name not in [
                    'validate', 'batch-validate', 'explain'  # samspec subcommands
                ]:
                    commands.add(command_name)
        
        return commands
    
    def extract_implemented_commands(self, cli_content: str) -> Set[str]:
        """Extract command names from argparse implementation."""
        commands = set()
        
        # Pattern to match ONLY main subparser.add_parser calls (not nested ones)
        # Look for lines that add parsers to the main subparsers object
        # Exclude lines with prefixes like "samspec_subparsers"
        parser_pattern = r'^[^#]*\bsubparsers\.add_parser\(\s*["\']([^"\']+)["\']'
        
        for line in cli_content.split('\n'):
            line = line.strip()
            match = re.match(parser_pattern, line)
            if match and 'samspec_subparsers' not in line:
                command_name = match.group(1)
                commands.add(command_name)
        
        return commands
    
    def extract_samspec_subcommands(self, cli_content: str) -> Set[str]:
        """Extract samspec subcommands for validation."""
        subcommands = set()
        
        # Look for samspec subparser creation
        samspec_pattern = r'samspec_subparsers\.add_parser\(\s*["\']([^"\']+)["\']'
        
        for match in re.finditer(samspec_pattern, cli_content):
            subcommand_name = match.group(1)
            subcommands.add(subcommand_name)
        
        return subcommands
    
    def test_cli_module_has_docstring(self, cli_module_docstring):
        """Test that CLI module has a docstring."""
        assert cli_module_docstring is not None
        assert len(cli_module_docstring.strip()) > 0
        assert "VCF Analysis Agent CLI" in cli_module_docstring
    
    def test_docstring_documents_all_implemented_commands(
        self, cli_module_content, cli_module_docstring
    ):
        """Test that all implemented commands are documented in the module docstring."""
        documented_commands = self.extract_documented_commands(cli_module_docstring)
        implemented_commands = self.extract_implemented_commands(cli_module_content)
        
        missing_from_docs = implemented_commands - documented_commands
        
        if missing_from_docs:
            pytest.fail(
                f"Commands implemented but not documented in module docstring: "
                f"{sorted(missing_from_docs)}. Update the module docstring to include these commands."
            )
    
    def test_no_phantom_commands_in_documentation(
        self, cli_module_content, cli_module_docstring
    ):
        """Test that documented commands are actually implemented."""
        documented_commands = self.extract_documented_commands(cli_module_docstring)
        implemented_commands = self.extract_implemented_commands(cli_module_content)
        
        documented_but_not_implemented = documented_commands - implemented_commands
        
        if documented_but_not_implemented:
            pytest.fail(
                f"Commands documented but not implemented: "
                f"{sorted(documented_but_not_implemented)}. Remove these from the docstring or implement them."
            )
    
    def test_command_examples_use_valid_commands(
        self, cli_module_content, cli_module_docstring
    ):
        """Test that command examples in docstring use actual implemented commands."""
        implemented_commands = self.extract_implemented_commands(cli_module_content)
        samspec_subcommands = self.extract_samspec_subcommands(cli_module_content)
        all_valid_commands = implemented_commands | samspec_subcommands
        
        # Pattern to match example command lines
        example_pattern = r'\$\s+python\s+-m\s+vcf_agent\.cli\s+([a-z][a-z0-9-]*)'
        
        invalid_examples = []
        for match in re.finditer(example_pattern, cli_module_docstring):
            command_used = match.group(1)
            
            if command_used not in all_valid_commands:
                invalid_examples.append(command_used)
        
        if invalid_examples:
            pytest.fail(
                f"Examples use undefined commands: {sorted(invalid_examples)}. "
                f"Valid commands are: {sorted(all_valid_commands)}"
            )
    
    def test_documentation_coverage_is_100_percent(
        self, cli_module_content, cli_module_docstring
    ):
        """Test that documentation coverage is 100%."""
        documented_commands = self.extract_documented_commands(cli_module_docstring)
        implemented_commands = self.extract_implemented_commands(cli_module_content)
        
        if not implemented_commands:
            pytest.skip("No implemented commands found")
        
        coverage = len(documented_commands) / len(implemented_commands) * 100
        
        assert coverage == 100.0, (
            f"Documentation coverage is {coverage:.1f}%, expected 100%. "
            f"Documented: {len(documented_commands)}, Implemented: {len(implemented_commands)}. "
            f"Missing: {sorted(implemented_commands - documented_commands)}"
        )
    
    def test_docstring_has_required_sections(self, cli_module_docstring):
        """Test that docstring has all required sections."""
        required_sections = [
            "Interactive Commands",
            "LanceDB Vector Database Operations", 
            "VCF Processing",
            "Compliance & Quality Validation",
            "Common Usage Examples"
        ]
        
        missing_sections = []
        for section in required_sections:
            if section not in cli_module_docstring:
                missing_sections.append(section)
        
        if missing_sections:
            pytest.fail(
                f"Docstring missing required sections: {missing_sections}. "
                f"Ensure the module docstring includes all standard sections."
            )
    
    def test_docstring_has_usage_examples(self, cli_module_docstring):
        """Test that docstring contains practical usage examples."""
        # Check for example patterns
        example_indicators = [
            "$ python -m vcf_agent.cli",
            "Usage Examples:",
            "Common Usage Examples:"
        ]
        
        has_examples = any(indicator in cli_module_docstring for indicator in example_indicators)
        
        assert has_examples, (
            "Docstring should contain usage examples. "
            "Add practical examples showing how to use the CLI commands."
        )
    
    def test_samspec_subcommands_are_documented(
        self, cli_module_content, cli_module_docstring
    ):
        """Test that samspec subcommands are mentioned in documentation."""
        samspec_subcommands = self.extract_samspec_subcommands(cli_module_content)
        
        if not samspec_subcommands:
            pytest.skip("No samspec subcommands found")
        
        # Check that samspec section exists and mentions subcommands
        samspec_section_exists = "samspec:" in cli_module_docstring.lower()
        
        assert samspec_section_exists, (
            "Docstring should document the samspec command and its subcommands. "
            f"Found subcommands: {sorted(samspec_subcommands)}"
        )
        
        # Check that each subcommand is mentioned
        for subcommand in samspec_subcommands:
            subcommand_mentioned = subcommand in cli_module_docstring
            assert subcommand_mentioned, (
                f"Samspec subcommand '{subcommand}' should be mentioned in docstring"
            )


class TestCLIDocumentationValidationScript:
    """Test the CLI documentation validation script itself."""
    
    def test_validation_script_exists(self):
        """Test that the validation script exists."""
        script_path = Path("scripts/validate_cli_documentation.py")
        assert script_path.exists(), "CLI documentation validation script should exist"
    
    def test_validation_script_is_executable(self):
        """Test that the validation script can be imported and run."""
        import sys
        import subprocess
        
        # Test that the script can be executed
        result = subprocess.run([
            sys.executable, 
            "scripts/validate_cli_documentation.py", 
            "--help"
        ], capture_output=True, text=True)
        
        assert result.returncode == 0, (
            f"Validation script should be executable. Error: {result.stderr}"
        )
        assert "CLI module documentation completeness" in result.stdout
    
    @pytest.mark.slow
    def test_validation_script_passes_on_current_cli(self):
        """Test that the validation script passes on the current CLI module."""
        import sys
        import subprocess
        
        # Run the validation script with completeness check
        result = subprocess.run([
            sys.executable,
            "scripts/validate_cli_documentation.py", 
            "--check-completeness"
        ], capture_output=True, text=True)
        
        if result.returncode != 0:
            pytest.fail(
                f"CLI documentation validation failed. "
                f"Output: {result.stdout}\nError: {result.stderr}"
            )


class TestCLIDocumentationMetrics:
    """Test CLI documentation metrics and quality."""
    
    def test_docstring_length_is_adequate(self, cli_module_docstring):
        """Test that docstring is comprehensive (adequate length)."""
        lines = [line.strip() for line in cli_module_docstring.split('\n') if line.strip()]
        
        # Should have substantial documentation
        assert len(lines) >= 30, (
            f"Docstring should be comprehensive. Found {len(lines)} lines, expected at least 30"
        )
    
    def test_docstring_has_multiple_examples(self, cli_module_docstring):
        """Test that docstring contains multiple usage examples."""
        example_pattern = r'\$\s+python\s+-m\s+vcf_agent\.cli'
        examples = re.findall(example_pattern, cli_module_docstring)
        
        assert len(examples) >= 5, (
            f"Docstring should contain multiple examples. Found {len(examples)}, expected at least 5"
        )
    
    def test_commands_are_categorized(self, cli_module_docstring):
        """Test that commands are organized into logical categories."""
        # Check for category headers
        category_indicators = [
            "Commands:",
            "Operations:",
            "Interactive",
            "LanceDB",
            "VCF Processing", 
            "Compliance"
        ]
        
        categories_found = sum(
            1 for indicator in category_indicators 
            if indicator in cli_module_docstring
        )
        
        assert categories_found >= 3, (
            f"Commands should be organized into categories. "
            f"Found {categories_found} category indicators, expected at least 3"
        )


# Integration test to ensure this works with existing test infrastructure
def test_cli_documentation_validation_integration():
    """Integration test to ensure CLI documentation validation works with pytest."""
    # This test validates that our validation logic works correctly
    from scripts.validate_cli_documentation import CLIDocumentationValidator
    
    try:
        validator = CLIDocumentationValidator()
        result = validator.validate()
        
        # The result should be valid (since we fixed the documentation)
        assert result.is_valid, (
            f"CLI documentation should be valid. Errors: {result.errors}"
        )
        
        # Should have documented commands
        assert len(result.documented_commands) > 0, "Should have documented commands"
        assert len(result.implemented_commands) > 0, "Should have implemented commands"
        
        # Coverage should be 100%
        coverage = len(result.documented_commands) / len(result.implemented_commands) * 100
        assert coverage == 100.0, f"Documentation coverage should be 100%, got {coverage:.1f}%"
        
    except FileNotFoundError:
        pytest.skip("CLI module not found - test environment issue")
    except Exception as e:
        pytest.fail(f"CLI documentation validation failed: {e}")


if __name__ == "__main__":
    # Allow running this test file directly
    pytest.main([__file__, "-v"]) 