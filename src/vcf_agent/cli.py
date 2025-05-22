"""
VCF Analysis Agent CLI Entrypoint

- Provides a command-line interface for interacting with the agent
- Accepts prompt strings for tool invocation (e.g., validation, echo)
- Supports mock response for testing
"""

import argparse
import os
from vcf_agent.agent import agent

def main():
    """
    Command-line interface for the VCF Analysis Agent.

    Accepts a prompt string to invoke agent tools (e.g., validation, echo).
    Supports mock response for testing via the VCF_AGENT_CLI_MOCK_RESPONSE environment variable.

    Args:
        None (uses sys.argv for CLI arguments)

    Returns:
        None (prints result to stdout)

    Example:
        $ python -m vcf_agent.cli "validate_vcf: sample_data/HG00098.vcf.gz"
        VALID: sample_data/HG00098.vcf.gz passed all checks.
    """
    mock_response = os.environ.get("VCF_AGENT_CLI_MOCK_RESPONSE")
    if mock_response is not None:
        print(mock_response)
        return
    parser = argparse.ArgumentParser(description="VCF Analysis Agent CLI")
    parser.add_argument("prompt", type=str, help="Prompt for the agent")
    args = parser.parse_args()

    response = agent(args.prompt)
    print(response)

if __name__ == "__main__":
    main() 