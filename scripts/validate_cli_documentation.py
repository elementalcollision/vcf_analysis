#!/usr/bin/env python3
"""
CLI Documentation Validation Script

Validates that the CLI module docstring accurately documents all implemented commands.
This prevents documentation drift and ensures developers can trust the module documentation.

Usage:
    python scripts/validate_cli_documentation.py [--check-completeness] [--verbose]
    
Exit codes:
    0: All documentation is complete and accurate
    1: Documentation discrepancies found
    2: Script execution error
"""

import sys
import re
import ast
import argparse
from pathlib import Path
from typing import Set, List, Dict, Any, Optional
from dataclasses import dataclass


@dataclass
class Command:
    """Represents a CLI command with its metadata."""
    name: str
    description: str
    category: Optional[str] = None
    is_subcommand: bool = False
    parent_command: Optional[str] = None


@dataclass
class ValidationResult:
    """Results of CLI documentation validation."""
    is_valid: bool
    documented_commands: Set[str]
    implemented_commands: Set[str]
    missing_from_docs: Set[str]
    documented_but_not_implemented: Set[str]
    errors: List[str]
    warnings: List[str]


class CLIDocumentationValidator:
    """Validates CLI module documentation against implementation."""
    
    def __init__(self, cli_module_path: str = "src/vcf_agent/cli.py"):
        self.cli_module_path = Path(cli_module_path)
        if not self.cli_module_path.exists():
            raise FileNotFoundError(f"CLI module not found at {cli_module_path}")
    
    def extract_commands_from_docstring(self) -> Set[str]:
        """Extract command names from the module docstring."""
        with open(self.cli_module_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Extract the module docstring
        tree = ast.parse(content)
        if not (tree.body and isinstance(tree.body[0], ast.Expr) 
                and isinstance(tree.body[0].value, ast.Constant)):
            raise ValueError("No module docstring found")
        
        docstring = tree.body[0].value.value
        commands = set()
        
        # Pattern to match command documentation lines like "• command-name: Description"
        command_pattern = r'^[•*-]\s+([a-z][a-z0-9-]*):?\s+'
        
        for line in docstring.split('\n'):
            line = line.strip()
            match = re.match(command_pattern, line)
            if match:
                command_name = match.group(1)
                # Filter out non-command entries
                if not command_name.startswith('--') and command_name not in [
                    'validate', 'batch-validate', 'explain'  # These are samspec subcommands
                ]:
                    commands.add(command_name)
        
        return commands
    
    def extract_commands_from_implementation(self) -> Set[str]:
        """Extract top-level command names from argparse implementation."""
        with open(self.cli_module_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        commands = set()
        
        # Pattern to match ONLY main subparser.add_parser calls (not nested ones)
        # Look for lines that add parsers to the main subparsers object
        # Exclude lines with prefixes like "samspec_subparsers"
        parser_pattern = r'^[^#]*\bsubparsers\.add_parser\(\s*["\']([^"\']+)["\']'
        
        for line in content.split('\n'):
            line = line.strip()
            match = re.match(parser_pattern, line)
            if match and 'samspec_subparsers' not in line:
                command_name = match.group(1)
                commands.add(command_name)
        
        return commands
    
    def extract_samspec_subcommands(self) -> Set[str]:
        """Extract samspec subcommands for special handling."""
        with open(self.cli_module_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        subcommands = set()
        
        # Look for samspec subparser creation
        samspec_pattern = r'samspec_subparsers\.add_parser\(\s*["\']([^"\']+)["\']'
        
        for match in re.finditer(samspec_pattern, content):
            subcommand_name = match.group(1)
            subcommands.add(subcommand_name)
        
        return subcommands
    
    def validate_command_examples(self) -> List[str]:
        """Validate that command examples in docstring use actual commands."""
        with open(self.cli_module_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Extract docstring
        tree = ast.parse(content)
        docstring = tree.body[0].value.value
        
        implemented_commands = self.extract_commands_from_implementation()
        samspec_subcommands = self.extract_samspec_subcommands()
        
        warnings = []
        
        # Pattern to match example command lines
        example_pattern = r'\$\s+python\s+-m\s+vcf_agent\.cli\s+([a-z][a-z0-9-]*)'
        
        for match in re.finditer(example_pattern, docstring):
            command_used = match.group(1)
            
            # Check if command exists (top-level or samspec subcommand)
            if (command_used not in implemented_commands and 
                command_used not in samspec_subcommands):
                warnings.append(f"Example uses undefined command: {command_used}")
        
        return warnings
    
    def validate_samspec_subcommands_documented(self) -> List[str]:
        """Validate that samspec subcommands are mentioned in docstring."""
        with open(self.cli_module_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Extract docstring
        tree = ast.parse(content)
        docstring = tree.body[0].value.value
        
        samspec_subcommands = self.extract_samspec_subcommands()
        warnings = []
        
        if samspec_subcommands:
            # Check that samspec section mentions the subcommands
            for subcommand in samspec_subcommands:
                if subcommand not in docstring:
                    warnings.append(f"Samspec subcommand '{subcommand}' not mentioned in docstring")
        
        return warnings
    
    def validate(self) -> ValidationResult:
        """Perform complete validation of CLI documentation."""
        errors = []
        warnings = []
        
        try:
            documented_commands = self.extract_commands_from_docstring()
        except Exception as e:
            errors.append(f"Failed to extract commands from docstring: {e}")
            documented_commands = set()
        
        try:
            implemented_commands = self.extract_commands_from_implementation()
        except Exception as e:
            errors.append(f"Failed to extract commands from implementation: {e}")
            implemented_commands = set()
        
        # Find discrepancies (only for top-level commands)
        missing_from_docs = implemented_commands - documented_commands
        documented_but_not_implemented = documented_commands - implemented_commands
        
        # Validate examples
        try:
            example_warnings = self.validate_command_examples()
            warnings.extend(example_warnings)
        except Exception as e:
            warnings.append(f"Failed to validate examples: {e}")
        
        # Validate samspec subcommands are mentioned
        try:
            samspec_warnings = self.validate_samspec_subcommands_documented()
            warnings.extend(samspec_warnings)
        except Exception as e:
            warnings.append(f"Failed to validate samspec subcommands: {e}")
        
        # Check for critical issues
        if missing_from_docs:
            errors.append(f"Commands implemented but not documented: {sorted(missing_from_docs)}")
        
        if documented_but_not_implemented:
            errors.append(f"Commands documented but not implemented: {sorted(documented_but_not_implemented)}")
        
        is_valid = len(errors) == 0
        
        return ValidationResult(
            is_valid=is_valid,
            documented_commands=documented_commands,
            implemented_commands=implemented_commands,
            missing_from_docs=missing_from_docs,
            documented_but_not_implemented=documented_but_not_implemented,
            errors=errors,
            warnings=warnings
        )


def print_validation_report(result: ValidationResult, verbose: bool = False) -> None:
    """Print a formatted validation report."""
    
    print("=" * 60)
    print("CLI DOCUMENTATION VALIDATION REPORT")
    print("=" * 60)
    
    if result.is_valid:
        print("✅ SUCCESS: CLI documentation is complete and accurate!")
        print(f"\n📊 Summary:")
        print(f"   • Total commands documented: {len(result.documented_commands)}")
        print(f"   • Total commands implemented: {len(result.implemented_commands)}")
        print(f"   • Documentation coverage: 100%")
    else:
        print("❌ FAILURE: CLI documentation has discrepancies!")
        print(f"\n📊 Summary:")
        print(f"   • Commands documented: {len(result.documented_commands)}")
        print(f"   • Commands implemented: {len(result.implemented_commands)}")
        print(f"   • Missing from docs: {len(result.missing_from_docs)}")
        print(f"   • Documented but not implemented: {len(result.documented_but_not_implemented)}")
        
        coverage = (len(result.documented_commands) / len(result.implemented_commands) * 100 
                   if result.implemented_commands else 0)
        print(f"   • Documentation coverage: {coverage:.1f}%")
    
    # Print errors
    if result.errors:
        print(f"\n🚨 ERRORS ({len(result.errors)}):")
        for i, error in enumerate(result.errors, 1):
            print(f"   {i}. {error}")
    
    # Print warnings
    if result.warnings:
        print(f"\n⚠️  WARNINGS ({len(result.warnings)}):")
        for i, warning in enumerate(result.warnings, 1):
            print(f"   {i}. {warning}")
    
    # Verbose output
    if verbose:
        print(f"\n📝 DETAILED ANALYSIS:")
        
        if result.documented_commands:
            print(f"\nDocumented Commands ({len(result.documented_commands)}):")
            for cmd in sorted(result.documented_commands):
                print(f"   • {cmd}")
        
        if result.implemented_commands:
            print(f"\nImplemented Commands ({len(result.implemented_commands)}):")
            for cmd in sorted(result.implemented_commands):
                print(f"   • {cmd}")
        
        if result.missing_from_docs:
            print(f"\nMissing from Documentation ({len(result.missing_from_docs)}):")
            for cmd in sorted(result.missing_from_docs):
                print(f"   ❌ {cmd}")
        
        if result.documented_but_not_implemented:
            print(f"\nDocumented but Not Implemented ({len(result.documented_but_not_implemented)}):")
            for cmd in sorted(result.documented_but_not_implemented):
                print(f"   ❌ {cmd}")
    
    print("=" * 60)


def main():
    """Main function for CLI documentation validation."""
    parser = argparse.ArgumentParser(
        description="Validate CLI module documentation completeness",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic validation
    python scripts/validate_cli_documentation.py
    
    # Detailed validation with verbose output
    python scripts/validate_cli_documentation.py --verbose
    
    # Validation for CI/CD (exits with non-zero on failure)
    python scripts/validate_cli_documentation.py --check-completeness
        """
    )
    
    parser.add_argument(
        "--check-completeness", 
        action="store_true",
        help="Exit with non-zero code if documentation is incomplete (for CI/CD)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true", 
        help="Show detailed validation information"
    )
    
    parser.add_argument(
        "--cli-module",
        default="src/vcf_agent/cli.py",
        help="Path to CLI module file (default: src/vcf_agent/cli.py)"
    )
    
    args = parser.parse_args()
    
    try:
        # Perform validation
        validator = CLIDocumentationValidator(args.cli_module)
        result = validator.validate()
        
        # Print report
        print_validation_report(result, verbose=args.verbose)
        
        # Exit with appropriate code
        if args.check_completeness and not result.is_valid:
            print(f"\n💥 VALIDATION FAILED: CLI documentation is incomplete!")
            print(f"   Fix the documentation discrepancies before deployment.")
            sys.exit(1)
        elif result.warnings and args.check_completeness:
            print(f"\n⚠️  VALIDATION PASSED WITH WARNINGS")
            print(f"   Consider addressing warnings for optimal documentation quality.")
            sys.exit(0)
        else:
            sys.exit(0)
            
    except FileNotFoundError as e:
        print(f"❌ ERROR: {e}", file=sys.stderr)
        sys.exit(2)
    except Exception as e:
        print(f"❌ UNEXPECTED ERROR: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main() 