"""
VCF Analysis Agent CLI Entrypoint

- Provides a command-line interface for interacting with the agent
- Accepts prompt strings for tool invocation (e.g., validation, echo)
- Supports mock response for testing
"""

import argparse
import os

def main():
    """
    Command-line interface for the VCF Analysis Agent.

    Accepts a prompt string to invoke agent tools (e.g., validation, echo).
    Supports mock response for testing via the VCF_AGENT_CLI_MOCK_RESPONSE environment variable.
    Supports --raw / --no-think flag to disable chain-of-thought reasoning (raw output mode).

    Args:
        None (uses sys.argv for CLI arguments)

    Returns:
        None (prints result to stdout)

    Example:
        $ python -m vcf_agent.cli "validate_vcf: sample_data/HG00098.vcf.gz"
        VALID: sample_data/HG00098.vcf.gz passed all checks.
        $ python -m vcf_agent.cli --raw "bcftools_view_tool: ['-h']"
    """
    mock_response = os.environ.get("VCF_AGENT_CLI_MOCK_RESPONSE")
    if mock_response is not None:
        print(mock_response)
        return
    parser = argparse.ArgumentParser(description="VCF Analysis Agent CLI")
    parser.add_argument("prompt", type=str, help="Prompt for the agent")
    parser.add_argument(
        "--raw", "--no-think", action="store_true", help="Disable chain-of-thought reasoning (raw output mode)"
    )
    args = parser.parse_args()

    if args.raw:
        os.environ["VCF_AGENT_RAW_MODE"] = "1"

    from vcf_agent.agent import agent  # Import after setting env var
    response = agent(args.prompt)
    print(response)

if __name__ == "__main__":
    main() 