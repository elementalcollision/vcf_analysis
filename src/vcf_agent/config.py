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

from typing import Optional

class SessionConfig:
    """
    Session configuration for VCF Analysis Agent.

    Controls session-specific settings, including output mode toggling (chain-of-thought vs. raw output).

    Attributes:
        raw_mode (Optional[bool]):
            - True: disables chain-of-thought reasoning (raw output mode)
            - False: enables chain-of-thought (CoT, default)
            - None: uses environment variable or CLI flag (see README)

    Usage:
        >>> from vcf_agent.config import SessionConfig
        >>> session_config = SessionConfig(raw_mode=True)
        >>> from vcf_agent.agent import get_agent_with_session
        >>> agent = get_agent_with_session(session_config)

    See the README for details on all output mode toggling methods.
    """
    def __init__(self, raw_mode: Optional[bool] = None):
        self.raw_mode = raw_mode 