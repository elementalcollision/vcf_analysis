"""
VCF Analysis Agent: Strands-based agent for VCF/BCF analysis, validation, and processing.

- Integrates bcftools and AI models for genomics workflows
- Provides CLI and API tools for validation and analysis
- Extensible with additional tools and models
"""

from strands import Agent, tools
from strands.tools import tool
from strands.models.ollama import OllamaModel
from typing import Any
from .validation import validate_vcf_file

SYSTEM_PROMPT = "You are the VCF Analysis Agent. You analyze, validate, and process VCF files for genomics workflows."

# Placeholder for tool imports and configuration
# from strands_tools import calculator, file_read, shell

# Placeholder echo tool using @tools.tool decorator
@tools.tool
def echo(text: str) -> str:
    """
    Echoes the input text back to the user.

    Args:
        text (str): Text to echo back.

    Returns:
        str: The echoed text.

    Example:
        >>> echo("Hello, world!")
        'Echo: Hello, world!'
    """
    return f"Echo: {text}"

# VCF validation tool
@tools.tool
def validate_vcf(filepath: str) -> str:
    """
    Validates a VCF/BCF file for existence, format, index, and bcftools stats.

    Args:
        filepath (str): Path to the VCF/BCF file.

    Returns:
        str: Validation result as a string (valid/invalid and error message if any).

    Example:
        >>> validate_vcf('sample_data/HG00098.vcf.gz')
        'VALID: sample_data/HG00098.vcf.gz passed all checks.'
    """
    is_valid, error = validate_vcf_file(filepath)
    if is_valid:
        return f"VALID: {filepath} passed all checks."
    else:
        return f"INVALID: {filepath} failed validation. Reason: {error}"

# Register the decorated tool directly
tools_list = [echo, validate_vcf]

# Configure Ollama model
ollama_model = OllamaModel(
    host="http://192.168.0.225:11434",
    model_id="qwen3:latest"
)

agent = Agent(
    system_prompt=SYSTEM_PROMPT,
    tools=tools_list,  # type: ignore # Suppress linter warning as runtime is OK
    model=ollama_model,
    callback_handler=None
)

__all__ = ["agent", "SYSTEM_PROMPT"] 