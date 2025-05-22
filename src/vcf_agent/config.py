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