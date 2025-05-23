"""
vcf_agent package

- Provides tools and wrappers for VCF/BCF analysis, validation, and compliance
- Integrates bcftools, GATK, and AI models
- Exposes CLI and API interfaces for genomics workflows
"""

__all__ = [
    "agent",
    "SYSTEM_PROMPT",
    "get_agent_with_session",
    "run_llm_analysis_task"
] 