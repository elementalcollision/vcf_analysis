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

class MemoryOptimizationConfig:
    """
    Configuration for memory optimization features in VCF Analysis Agent.
    
    This replaces the artificial "phase" approach with production-ready optimization levels.
    """
    def __init__(
        self,
        optimization_level: Literal["basic", "standard", "aggressive"] = "standard",
        memory_management_enabled: bool = True,
        dimension_reduction_enabled: bool = True,
        target_dimensions: int = 768,
        caching_strategy: Literal["simple", "memory_aware"] = "memory_aware",
        cache_max_size_mb: int = 40,
        cache_max_entries: int = 1000,
        streaming_batch_size: int = 25,
        memory_cleanup_threshold_mb: int = 40
    ):
        """
        Initialize memory optimization configuration.
        
        Args:
            optimization_level: Overall optimization level
                - "basic": PyArrow optimizations only
                - "standard": Memory management + basic dimension reduction
                - "aggressive": Full optimizations including streaming
            memory_management_enabled: Enable memory-aware caching and cleanup
            dimension_reduction_enabled: Enable PCA dimension reduction for embeddings
            target_dimensions: Target dimensions for reduced embeddings (768 for 50% reduction)
            caching_strategy: Caching approach for embeddings
            cache_max_size_mb: Maximum cache size in megabytes
            cache_max_entries: Maximum number of cached entries
            streaming_batch_size: Batch size for streaming embedding generation
            memory_cleanup_threshold_mb: Memory threshold for automatic cleanup
        """
        self.optimization_level = optimization_level
        self.memory_management_enabled = memory_management_enabled
        self.dimension_reduction_enabled = dimension_reduction_enabled
        self.target_dimensions = target_dimensions
        self.caching_strategy = caching_strategy
        self.cache_max_size_mb = cache_max_size_mb
        self.cache_max_entries = cache_max_entries
        self.streaming_batch_size = streaming_batch_size
        self.memory_cleanup_threshold_mb = memory_cleanup_threshold_mb
        
        # Auto-configure based on optimization level
        if optimization_level == "basic":
            self.memory_management_enabled = False
            self.dimension_reduction_enabled = False
            self.caching_strategy = "simple"
        elif optimization_level == "aggressive":
            self.memory_management_enabled = True
            self.dimension_reduction_enabled = True
            self.caching_strategy = "memory_aware"
            self.streaming_batch_size = 50  # Larger batches for aggressive mode
    
    def get_embedding_dimensions(self) -> int:
        """Get the appropriate embedding dimensions based on configuration."""
        if self.dimension_reduction_enabled:
            return self.target_dimensions
        return 1536  # Original OpenAI embedding dimensions
    
    def is_optimized_model_enabled(self) -> bool:
        """Check if optimized vector models should be used."""
        return self.dimension_reduction_enabled

class SessionConfig:
    """
    Session configuration for VCF Analysis Agent.

    Controls session-specific settings, including output mode toggling (chain-of-thought vs. raw output),
    model provider selection, reference FASTA for normalization, and memory optimizations.

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
        memory_optimization (MemoryOptimizationConfig):
            - Configuration for memory optimization features

    Usage:
        >>> from vcf_agent.config import SessionConfig, MemoryOptimizationConfig
        >>> memory_config = MemoryOptimizationConfig(optimization_level="aggressive")
        >>> session_config = SessionConfig(raw_mode=True, model_provider="openai", memory_optimization=memory_config)
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
        ollama_model_name: Optional[str] = "qwen3:4b",  # Default to qwen3:4b
        ollama_base_url: Optional[str] = "http://localhost:11434",
        memory_optimization: Optional[MemoryOptimizationConfig] = None
    ):
        self.raw_mode = raw_mode
        self.model_provider = model_provider
        self.credentials_file = credentials_file
        self.reference_fasta = reference_fasta
        self.ollama_model_name = ollama_model_name # Store it
        self.ollama_base_url = ollama_base_url
        self.memory_optimization = memory_optimization or MemoryOptimizationConfig()

    def __repr__(self) -> str:
        """String representation of the session config."""
        return (
            f"SessionConfig(raw_mode={self.raw_mode}, "
            f"model_provider='{self.model_provider}', "
            f"credentials_file={repr(self.credentials_file)}, "
            f"reference_fasta={repr(self.reference_fasta)}, "
            f"ollama_model_name='{self.ollama_model_name}', "
            f"ollama_base_url={repr(self.ollama_base_url)}, "
            f"memory_optimization={self.memory_optimization.optimization_level})"
        ) 