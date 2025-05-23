"""
VCF Analysis Agent CLI Entrypoint

- Provides a command-line interface for interacting with the agent
- Accepts prompt strings for tool invocation (e.g., validation, echo)
- Supports mock response for testing
- Supports multiple model providers (Ollama, OpenAI, Cerebras)
"""

import argparse
import os
from typing import Literal, cast

def main():
    """
    Command-line interface for the VCF Analysis Agent.

    Accepts a prompt string to invoke agent tools (e.g., validation, echo).
    Supports mock response for testing via the VCF_AGENT_CLI_MOCK_RESPONSE environment variable.
    Supports --raw / --no-think flag to disable chain-of-thought reasoning (raw output mode).
    Supports --model to select the model provider (ollama, openai, cerebras).

    Args:
        None (uses sys.argv for CLI arguments)

    Returns:
        None (prints result to stdout)

    Example:
        $ python -m vcf_agent.cli "validate_vcf: sample_data/HG00098.vcf.gz"
        VALID: sample_data/HG00098.vcf.gz passed all checks.
        $ python -m vcf_agent.cli --raw "bcftools_view_tool: ['-h']"
        $ python -m vcf_agent.cli --model openai "echo: Hello from OpenAI!"
    """
    mock_response = os.environ.get("VCF_AGENT_CLI_MOCK_RESPONSE")
    if mock_response is not None:
        print(mock_response)
        return
        
    parser = argparse.ArgumentParser(description="VCF Analysis Agent CLI")
    parser.add_argument("prompt", type=str, help="Prompt for the agent")
    parser.add_argument(
        "--raw", "--no-think", action="store_true", 
        help="Disable chain-of-thought reasoning (raw output mode)"
    )
    parser.add_argument(
        "--model", type=str, choices=["ollama", "openai", "cerebras"], default="ollama",
        help="Select model provider (default: ollama)"
    )
    parser.add_argument(
        "--credentials", type=str, 
        help="Path to JSON credentials file for API access"
    )
    parser.add_argument(
        "--reference", "--fasta", type=str, default=None,
        help="Path to reference FASTA for normalization (required for some tools)"
    )
    args = parser.parse_args()

    if args.raw:
        os.environ["VCF_AGENT_RAW_MODE"] = "1"

    # Import after setting env vars
    from vcf_agent.agent import get_agent_with_session
    from vcf_agent.config import SessionConfig
    
    # Create session config with CLI options
    model_provider = cast(Literal["ollama", "openai", "cerebras"], args.model)
    session_config = SessionConfig(
        raw_mode=args.raw if args.raw else None,
        model_provider=model_provider,
        credentials_file=args.credentials,
        reference_fasta=args.reference
    )
    
    # Get agent with specified model
    agent = get_agent_with_session(
        session_config=session_config,
        model_provider=model_provider
    )
    
    # Run prompt
    response = agent(args.prompt)
    print(response)

if __name__ == "__main__":
    main() 