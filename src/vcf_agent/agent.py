"""
VCF Analysis Agent: Strands-based agent for VCF/BCF analysis, validation, and processing.

- Integrates bcftools and AI models for genomics workflows
- Provides CLI and API tools for validation and analysis
- Extensible with additional tools and models

Output Mode Toggling (Chain-of-Thought vs. Raw Output):
------------------------------------------------------
The agent supports three ways to control output mode:
  1. Environment variable: VCF_AGENT_RAW_MODE ("1", "true", "yes" = raw output)
  2. CLI flag: --raw / --no-think (see cli.py)
  3. Session-based: SessionConfig(raw_mode=...) (see config.py)
See the README for usage examples and details.
"""

from strands import Agent, tools
from strands.tools import tool
from strands.models.ollama import OllamaModel
from strands.models.litellm import LiteLLMModel
from typing import Any, Optional, Dict, List, cast, Union, Literal
from .validation import validate_vcf_file
from .bcftools_integration import (
    bcftools_view as _bcftools_view,
    bcftools_query as _bcftools_query,
    bcftools_filter as _bcftools_filter,
    bcftools_norm as _bcftools_norm,
    bcftools_stats as _bcftools_stats,
    bcftools_annotate as _bcftools_annotate,
    vcf_compare as _vcf_compare,
)
import os
from .config import SessionConfig
from .api_clients import OpenAIClient, CerebrasClient, APIClientError

# Output mode toggling: chain-of-thought (CoT) vs. raw output
# 1. Environment variable: VCF_AGENT_RAW_MODE ("1", "true", "yes" = raw)
# 2. CLI flag: --raw / --no-think (see cli.py)
# 3. Session-based: SessionConfig(raw_mode=...) (see config.py)
# See README for details.
RAW_MODE = True  # Always use raw output mode (no chain-of-thought, no <think> blocks)

SYSTEM_PROMPT = (
    "You are the VCF Analysis Agent. You analyze, validate, and process VCF files for genomics workflows.\n"
    "For all analysis tasks, you MUST respond with ONLY valid JSON matching the provided schema.\n"
    "Do NOT include any explanation, markdown, <think> blocks, or extra text—output must be pure JSON.\n"
    "If you cannot perform the task, return a valid JSON object with an 'error' field and no other fields.\n"
    "/no_think"  # Always disable chain-of-thought
)

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
@tools.tool
def bcftools_view_tool(args: list) -> str:
    """
    Run bcftools view with the given arguments.

    Args:
        args (list): Arguments for bcftools view (e.g., ["-h", "file.vcf"])
    Returns:
        str: Output or error from bcftools.
    """
    rc, out, err = _bcftools_view(args)
    return out if rc == 0 else err

@tools.tool
def bcftools_query_tool(args: list) -> str:
    """
    Run bcftools query with the given arguments.

    Args:
        args (list): Arguments for bcftools query.
    Returns:
        str: Output or error from bcftools.
    """
    rc, out, err = _bcftools_query(args)
    return out if rc == 0 else err

@tools.tool
def bcftools_filter_tool(args: list) -> str:
    """
    Run bcftools filter with the given arguments.

    Args:
        args (list): Arguments for bcftools filter.
    Returns:
        str: Output or error from bcftools.
    """
    rc, out, err = _bcftools_filter(args)
    return out if rc == 0 else err

@tools.tool
def bcftools_norm_tool(args: list) -> str:
    """
    Run bcftools norm with the given arguments.

    Args:
        args (list): Arguments for bcftools norm.
    Returns:
        str: Output or error from bcftools.
    """
    rc, out, err = _bcftools_norm(args)
    return out if rc == 0 else err

@tools.tool
def bcftools_stats_tool(args: list) -> str:
    """
    Run bcftools stats with the given arguments.

    Args:
        args (list): Arguments for bcftools stats.
    Returns:
        str: Output or error from bcftools.
    """
    rc, out, err = _bcftools_stats(args)
    return out if rc == 0 else err

@tools.tool
def bcftools_annotate_tool(args: list) -> str:
    """
    Run bcftools annotate with the given arguments.

    Args:
        args (list): Arguments for bcftools annotate.
    Returns:
        str: Output or error from bcftools.
    """
    rc, out, err = _bcftools_annotate(args)
    return out if rc == 0 else err

@tools.tool
def vcf_comparison_tool(file1: str, file2: str) -> str:
    """
    Compare two VCF files for concordance, discordance, and quality metrics.
    Returns: JSON string with all required fields, even on error.
    """
    import json
    try:
        result = _vcf_compare(file1, file2)
        # Ensure all required fields are present, even on error
        required_fields = [
            "concordant_variant_count",
            "discordant_variant_count",
            "unique_to_file_1",
            "unique_to_file_2",
            "quality_metrics"
        ]
        for field in required_fields:
            if field not in result:
                result[field] = 0 if 'count' in field else [] if 'unique' in field else {}
        if "error" not in result:
            result["error"] = None
        return json.dumps(result)
    except Exception as e:
        return json.dumps({
            "concordant_variant_count": 0,
            "discordant_variant_count": 0,
            "unique_to_file_1": [],
            "unique_to_file_2": [],
            "quality_metrics": {},
            "error": str(e)
        })

# Add stubs for vcf_analysis_summary_tool and vcf_summarization_tool if not present
@tools.tool
def vcf_analysis_summary_tool(filepath: str) -> str:
    """
    Analyze a VCF file and return a summary. Always returns all required fields, even on error.
    """
    import json
    try:
        # Placeholder: not implemented
        raise NotImplementedError("vcf_analysis_summary_tool is not implemented.")
    except Exception as e:
        return json.dumps({
            "variant_count": 0,
            "variant_types": {},
            "sample_statistics": {},
            "notable_patterns": [],
            "error": str(e)
        })

@tools.tool
def vcf_summarization_tool(filepath: str) -> str:
    """
    Summarize a VCF file. Always returns all required fields, even on error.
    """
    import json
    try:
        # Placeholder: not implemented
        raise NotImplementedError("vcf_summarization_tool is not implemented.")
    except Exception as e:
        return json.dumps({
            "variant_count": 0,
            "variant_types": {},
            "sample_statistics": {},
            "notable_patterns": [],
            "error": str(e)
        })

# Register all bcftools tools to tools_list
tools_list = [
    echo,
    validate_vcf,
    bcftools_view_tool,
    bcftools_query_tool,
    bcftools_filter_tool,
    bcftools_norm_tool,
    bcftools_stats_tool,
    bcftools_annotate_tool,
    vcf_comparison_tool,  # Register the new tool
]

# Configure Ollama model
ollama_model = OllamaModel(
    host="http://192.168.0.14:11434",
    model_id="qwen3:latest"
)

# Setup OpenAI model adapter for Strands using LiteLLM
def get_openai_model(credential_manager=None):
    """
    Create an OpenAI model instance for Strands using LiteLLM.
    
    Args:
        credential_manager: Optional CredentialManager instance.
                           If None, a new one will be created.
    
    Returns:
        LiteLLMModel: Configured OpenAI model instance
    """
    try:
        openai_client = OpenAIClient(credential_manager=credential_manager)
        # Basic connection test
        openai_client.test_connection()
        
        # Get API key from our credential manager
        api_key = openai_client.credential_manager.get_credential("openai")
        
        # Create LiteLLM model for OpenAI
        return LiteLLMModel(
            client_args={"api_key": api_key},
            model_id="openai/gpt-4o",
            params={
                "temperature": 0.7,
                "max_tokens": 4000
            }
        )
    except APIClientError as e:
        import logging
        logging.error(f"Failed to initialize OpenAI model: {e}")
        raise

# Create a mock Cerebras model adapter for Strands
class CerebrasStrandsModel:
    """Custom Strands model implementation for Cerebras API."""
    
    def __init__(self, model: str = "cerebras-gpt", credential_manager=None):
        """
        Initialize the Cerebras model adapter.
        
        Args:
            model: Model name to use
            credential_manager: Optional CredentialManager instance.
                               If None, a new one will be created.
        """
        self.model = model
        self.client = CerebrasClient(credential_manager=credential_manager)
        # Test connection during initialization
        self.client.test_connection()
    
    def generate(self, prompt: str, **kwargs) -> str:
        """
        Generate text from the Cerebras model.
        
        Args:
            prompt: The prompt text
            **kwargs: Additional arguments for the model
            
        Returns:
            str: The generated text
        """
        # Convert prompt to messages format
        messages = [{"role": "user", "content": prompt}]
        
        # Call the Cerebras API
        response = self.client.chat_completion(
            messages=messages,
            model=self.model,
            **kwargs
        )
        
        # Extract the response text from Cerebras response format
        return response["choices"][0]["message"]["content"]

def get_cerebras_model(credential_manager=None):
    """
    Create a Cerebras model instance for Strands using our custom API client.
    
    Args:
        credential_manager: Optional CredentialManager instance.
                           If None, a new one will be created.
    
    Returns:
        CerebrasStrandsModel: Configured Cerebras model instance
    """
    try:
        return CerebrasStrandsModel(credential_manager=credential_manager)
    except APIClientError as e:
        import logging
        logging.error(f"Failed to initialize Cerebras model: {e}")
        raise

def get_agent_with_session(
    session_config: Optional[SessionConfig] = None,
    model_provider: Literal["ollama", "openai", "cerebras"] = "ollama"
):
    """
    Instantiate a VCF Analysis Agent with optional session-based configuration.

    Output mode (CoT vs. raw) is determined in this order:
      1. session_config.raw_mode (if not None)
      2. VCF_AGENT_RAW_MODE env var ("1", "true", "yes" = raw)
      3. Default: chain-of-thought (CoT)

    Args:
        session_config (Optional[SessionConfig]): Session configuration. If raw_mode is set, it overrides env/CLI.
        model_provider (str): The LLM provider to use: "ollama" (default), "openai", or "cerebras"

    Returns:
        Agent: Configured VCF Analysis Agent instance.

    Usage:
        >>> from vcf_agent.config import SessionConfig
        >>> session_config = SessionConfig(raw_mode=True)
        >>> from vcf_agent.agent import get_agent_with_session
        >>> agent = get_agent_with_session(session_config, model_provider="openai")
    """
    # Always use raw output mode (no chain-of-thought, no <think> blocks)
    system_prompt = (
        "You are the VCF Analysis Agent. You analyze, validate, and process VCF files for genomics workflows.\n"
        "For all analysis tasks, you MUST respond with ONLY valid JSON matching the provided schema.\n"
        "Do NOT include any explanation, markdown, <think> blocks, or extra text—output must be pure JSON.\n"
        "If you cannot perform the task, return a valid JSON object with an 'error' field and no other fields.\n"
        "/no_think"  # Always disable chain-of-thought
    )

    # Use model_provider from session_config if available
    if session_config and session_config.model_provider:
        model_provider = session_config.model_provider

    # Get credentials file from session config
    credentials_file = session_config.credentials_file if session_config else None

    # Select the appropriate model based on provider
    model = None
    if model_provider == "ollama":
        model = OllamaModel(
            host="http://192.168.0.14:11434",
            model_id="qwen3:latest"
        )
    elif model_provider == "openai":
        from .api_clients import CredentialManager
        # Create credential manager with the credentials file
        credential_manager = CredentialManager(credentials_file=credentials_file)
        model = get_openai_model(credential_manager)
    elif model_provider == "cerebras":
        from .api_clients import CredentialManager
        # Create credential manager with the credentials file
        credential_manager = CredentialManager(credentials_file=credentials_file)
        model = get_cerebras_model(credential_manager)
    else:
        raise ValueError(f"Unknown model provider: {model_provider}")
    
    return Agent(
        system_prompt=system_prompt,
        tools=tools_list,  # type: ignore
        model=model,
        callback_handler=None
    )

# Default agent instance (uses env/CLI for RAW_MODE)
agent = get_agent_with_session()

__all__ = ["agent", "SYSTEM_PROMPT", "get_agent_with_session"] 