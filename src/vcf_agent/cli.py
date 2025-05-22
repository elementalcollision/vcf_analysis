import argparse
from vcf_agent.agent import agent

def main():
    parser = argparse.ArgumentParser(description="VCF Analysis Agent CLI")
    parser.add_argument("prompt", type=str, help="Prompt for the agent")
    args = parser.parse_args()

    response = agent(args.prompt)
    print(response)

if __name__ == "__main__":
    main() 