#!/usr/bin/env python3
"""Test direct tool calling functionality."""

import sys
import tempfile
import os
from pathlib import Path

# Add the src directory to the Python path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from vcf_agent.agent import get_agent_with_session
from vcf_agent.config import SessionConfig

def test_direct_tools():
    """Test direct tool calling."""
    print("=== Testing Direct Tool Calling ===")
    
    # Create a test VCF file
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	60	PASS	.
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_content)
        test_vcf = f.name
    
    print(f"Created test VCF: {test_vcf}")
    
    try:
        # Create agent
        agent = get_agent_with_session(
            session_config=SessionConfig(raw_mode=False),
            model_provider="ollama"
        )
        
        # Test direct tool calling
        print("\nTesting direct tool call...")
        try:
            result = agent.tool.validate_vcf(filepath=test_vcf)
            print(f"✅ Direct tool call result: {result}")
        except Exception as e:
            print(f"❌ Direct tool call failed: {e}")
        
        # Test echo tool
        print("\nTesting echo tool...")
        try:
            echo_result = agent.tool.echo(text="Hello from direct call!")
            print(f"✅ Echo result: {echo_result}")
        except Exception as e:
            print(f"❌ Echo failed: {e}")
            
        # Test available tools
        print("\nChecking available tools...")
        try:
            if hasattr(agent, 'tool'):
                print("✅ agent.tool exists")
                tool_attrs = [attr for attr in dir(agent.tool) if not attr.startswith('_')]
                print(f"Available tool methods: {tool_attrs}")
            else:
                print("❌ agent.tool does not exist")
        except Exception as e:
            print(f"❌ Error checking tools: {e}")
            
    finally:
        # Cleanup
        os.unlink(test_vcf)
        print(f"Cleaned up test file: {test_vcf}")

if __name__ == "__main__":
    test_direct_tools() 