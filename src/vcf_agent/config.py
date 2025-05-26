"""
Configuration module for VCF Analysis Agent.

- Stores global configuration for compliance tools and references
- Can be extended for environment management and advanced settings
"""

CONFIG = {
    # Compliance checker settings
    "compliance_tool": "bcftools",  # Options: "bcftools", "gatk", or manufacturer/tool name
    "gatk_reference": None,           # Path to reference FASTA for GATK (if needed)
    "manufacturer_tool": None,        # Name of manufacturer-specific tool (e.g., "ont_vcf_validator"),
    # Add more configuration options as needed
}
"""
Global configuration dictionary for the VCF Analysis Agent.

Keys:
    compliance_tool (str): Compliance tool to use ('bcftools', 'gatk', or manufacturer/tool name)
    gatk_reference (Optional[str]): Path to reference FASTA for GATK
    manufacturer_tool (Optional[str]): Name of manufacturer-specific tool

Usage:
    >>> from vcf_agent.config import CONFIG
    >>> CONFIG["compliance_tool"]
    'bcftools'
"""

from typing import Optional, Literal

class SessionConfig:
    """
    Session configuration for VCF Analysis Agent.

    Controls session-specific settings, including output mode toggling (chain-of-thought vs. raw output),
    model provider selection, and reference FASTA for normalization.

    Attributes:
        raw_mode (Optional[bool]):
            - True: disables chain-of-thought reasoning (raw output mode)
            - False: enables chain-of-thought (CoT, default)
            - None: uses environment variable or CLI flag (see README)
        model_provider (str):
            - "ollama": Use local Ollama models (default)
            - "openai": Use OpenAI API
            - "cerebras": Use Cerebras API
        credentials_file (Optional[str]):
            - Path to JSON file with API credentials
            - If None, will use environment variables
        reference_fasta (Optional[str]):
            - Path to reference FASTA for normalization (required for some tools)

    Usage:
        >>> from vcf_agent.config import SessionConfig
        >>> session_config = SessionConfig(raw_mode=True, model_provider="openai", reference_fasta="ref.fa")
        >>> from vcf_agent.agent import get_agent_with_session
        >>> agent = get_agent_with_session(session_config)

    See the README for details on all output mode toggling methods and model providers.
    """
    def __init__(
        self, 
        raw_mode: Optional[bool] = None,
        model_provider: Literal["ollama", "openai", "cerebras"] = "ollama",
        credentials_file: Optional[str] = None,
        reference_fasta: Optional[str] = None,
        ollama_model_name: Optional[str] = "qwen:4b"  # Default to qwen:4b
    ):
        self.raw_mode = raw_mode
        self.model_provider = model_provider
        self.credentials_file = credentials_file
        self.reference_fasta = reference_fasta
        self.ollama_model_name = ollama_model_name # Store it

    def __repr__(self) -> str:
        """String representation of the session config."""
        return (
            f"SessionConfig(raw_mode={self.raw_mode}, "
            f"model_provider='{self.model_provider}', "
            f"credentials_file={repr(self.credentials_file)}, "
            f"reference_fasta={repr(self.reference_fasta)}, "
            f"ollama_model_name='{self.ollama_model_name}')"
        ) 