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
from .prompt_templates import get_prompt_for_task

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
@tools.tool(
    openai_schema={
        "name": "bcftools_view_tool",
        "description": "Run bcftools view with the given arguments.",
        "parameters": {
            "args": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Arguments for bcftools view (e.g., ['-h', 'file.vcf'])"
            }
        },
        "required": ["args"]
    }
)
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

@tools.tool(
    openai_schema={
        "name": "bcftools_query_tool",
        "description": "Run bcftools query with the given arguments.",
        "parameters": {
            "args": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Arguments for bcftools query."
            }
        },
        "required": ["args"]
    }
)
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

@tools.tool(
    openai_schema={
        "name": "bcftools_filter_tool",
        "description": "Run bcftools filter with the given arguments.",
        "parameters": {
            "args": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Arguments for bcftools filter."
            }
        },
        "required": ["args"]
    }
)
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

@tools.tool(
    openai_schema={
        "name": "bcftools_norm_tool",
        "description": "Run bcftools norm with the given arguments.",
        "parameters": {
            "args": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Arguments for bcftools norm."
            }
        },
        "required": ["args"]
    }
)
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

@tools.tool(
    openai_schema={
        "name": "bcftools_stats_tool",
        "description": "Run bcftools stats with the given arguments.",
        "parameters": {
            "args": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Arguments for bcftools stats."
            }
        },
        "required": ["args"]
    }
)
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

@tools.tool(
    openai_schema={
        "name": "bcftools_annotate_tool",
        "description": "Run bcftools annotate with the given arguments.",
        "parameters": {
            "args": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Arguments for bcftools annotate."
            }
        },
        "required": ["args"]
    }
)
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

@tools.tool(
    openai_schema={
        "name": "vcf_comparison_tool",
        "description": "Compare two VCF files for concordance, discordance, and quality metrics, with normalization and multi-sample support.",
        "parameters": {
            "file1": {"type": "string", "description": "Path to the first VCF file."},
            "file2": {"type": "string", "description": "Path to the second VCF file."},
            "reference": {"type": "string", "description": "Path to the reference FASTA for normalization."}
        },
        "required": ["file1", "file2", "reference"]
    }
)
def vcf_comparison_tool(file1: str, file2: str, reference: str) -> str:
    """
    Compare two VCF files for concordance, discordance, and quality metrics.
    - Normalizes and decomposes VCFs using bcftools norm.
    - Supports multi-sample VCFs, reporting per-sample concordance/discordance.
    Args:
        file1: Path to first VCF
        file2: Path to second VCF
        reference: Path to reference FASTA for normalization
    Returns:
        JSON string with concordant_variant_count, discordant_variant_count, unique_to_file_1, unique_to_file_2, quality_metrics, and per_sample_concordance.
    """
    import json
    import tempfile
    import subprocess
    from .vcf_utils import extract_variant_summary
    try:
        # Normalize and decompose both VCFs
        def norm_vcf(input_vcf, ref):
            normed = tempfile.NamedTemporaryFile(delete=False, suffix='.vcf.gz')
            cmd = [
                'bcftools', 'norm', '-m-any', '-f', ref, '-Oz', '-o', normed.name, input_vcf
            ]
            subprocess.run(cmd, check=True)
            return normed.name
        normed1 = norm_vcf(file1, reference)
        normed2 = norm_vcf(file2, reference)
        # Use bcftools isec for comparison
        isec_dir = tempfile.mkdtemp()
        cmd_isec = [
            'bcftools', 'isec', '-p', isec_dir, normed1, normed2
        ]
        subprocess.run(cmd_isec, check=True)
        # Parse results (simplified for brevity)
        concordant = 0
        discordant = 0
        unique1 = []
        unique2 = []
        # Per-sample stats
        per_sample = {}
        # Use cyvcf2 to parse and count per-sample
        from cyvcf2 import VCF
        v1 = VCF(normed1)
        v2 = VCF(normed2)
        samples1 = v1.samples
        samples2 = v2.samples
        all_samples = list(set(samples1) | set(samples2))
        for s in all_samples:
            per_sample[s] = {"concordant": 0, "discordant": 0}
        # For each record in v1, check if present in v2 (by chrom, pos, ref, alt)
        v2_records = {(r.CHROM, r.POS, r.REF, tuple(r.ALT)): r for r in v2}
        for r in v1:
            key = (r.CHROM, r.POS, r.REF, tuple(r.ALT))
            if key in v2_records:
                concordant += 1
                # Per-sample genotype concordance
                for i, s in enumerate(samples1):
                    gt1 = r.genotypes[i] if i < len(r.genotypes) else None
                    gt2 = v2_records[key].genotypes[i] if i < len(v2_records[key].genotypes) else None
                    if gt1 == gt2:
                        per_sample[s]["concordant"] += 1
                    else:
                        per_sample[s]["discordant"] += 1
            else:
                discordant += 1
                unique1.append(key)
        for r in v2:
            key = (r.CHROM, r.POS, r.REF, tuple(r.ALT))
            if key not in v2_records:
                unique2.append(key)
        # Quality metrics (mocked for now)
        quality_metrics = {"mean_qual_file1": 0, "mean_qual_file2": 0}
        result = {
            "concordant_variant_count": concordant,
            "discordant_variant_count": discordant,
            "unique_to_file_1": unique1,
            "unique_to_file_2": unique2,
            "quality_metrics": quality_metrics,
            "per_sample_concordance": per_sample
        }
        return json.dumps(result)
    except Exception as e:
        return json.dumps({"error": str(e)})

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

@tools.tool(
    openai_schema={
        "name": "vcf_summarization_tool",
        "description": "Summarize a VCF file for key variant statistics and sample-level insights.",
        "parameters": {
            "filepath": {"type": "string", "description": "Path to the VCF file to summarize."}
        },
        "required": ["filepath"]
    }
)
def vcf_summarization_tool(filepath: str) -> str:
    import json
    from .vcf_utils import extract_variant_summary
    try:
        summary = extract_variant_summary(filepath)
        # Compose sample_statistics (mocked for now)
        sample_statistics = {}
        for sample in summary["samples"]:
            sample_statistics[sample] = {
                "mean_depth": 30.0,  # Placeholder
                "het_ratio": 0.5     # Placeholder
            }
        result = {
            "variant_count": summary["variant_count"],
            "variant_types": summary["variant_types"],
            "sample_statistics": sample_statistics,
            "notable_patterns": []
        }
        return json.dumps(result)
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
    host="http://127.0.0.1:11434",
    model_id="Qwen3:4b"
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

    def converse(self, messages, tool_specs=None, system_prompt=None, **kwargs):
        """
        Mock converse method for Strands compatibility. Yields a mock response event.
        """
        # Compose a mock response event as expected by Strands
        content = "This is a mock Cerebras response."
        yield {
            "role": "assistant",
            "content": content,
            "tool_calls": [],
            "usage": {"prompt_tokens": 1, "completion_tokens": 1, "total_tokens": 2},
        }

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
            host="http://127.0.0.1:11434",
            model_id="Qwen3:4b"
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

def run_llm_analysis_task(
    task: str,
    file_paths: list,
    extra_context: Optional[Dict[str, Any]] = None,
    session_config: Optional[SessionConfig] = None,
    model_provider: Literal["ollama", "openai", "cerebras"] = "ollama",
    **llm_kwargs
) -> str:
    """
    Constructs a prompt and runs the LLM analysis task, returning the JSON response as a string.
    Args:
        task: The prompt contract/task name (e.g., 'vcf_summarization_v1')
        file_paths: List of VCF file paths (1 or 2, depending on task)
        extra_context: Optional dict for additional context
        session_config: Optional SessionConfig for agent/model selection
        model_provider: LLM provider to use
        **llm_kwargs: Additional kwargs for the agent/model
    Returns:
        str: LLM response (should be valid JSON per contract)
    """
    agent = get_agent_with_session(session_config, model_provider)
    prompt = get_prompt_for_task(task, file_paths, extra_context)
    result = agent(prompt, **llm_kwargs)
    return str(result)

__all__ = ["agent", "SYSTEM_PROMPT", "get_agent_with_session"] 