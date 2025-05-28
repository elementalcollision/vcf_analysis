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

def test_direct_tool_calling():
    """Test direct tool calling."""
    print("=== Testing Direct Tool Calling ===")
    
    # Create a simple test VCF
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
        agent = get_agent_with_session(SessionConfig(raw_mode=False), "ollama")
        print("✅ Agent created successfully")
        
        # Check if agent has tool attribute
        print(f"Agent has 'tool' attribute: {hasattr(agent, 'tool')}")
        
        # Check if validate_vcf is available as direct attribute
        print(f"Agent has 'validate_vcf' attribute: {hasattr(agent, 'validate_vcf')}")
        
        # List all agent attributes
        agent_attrs = [attr for attr in dir(agent) if not attr.startswith('_')]
        print(f"Agent attributes: {agent_attrs}")
        
        # Try direct tool call via attribute
        if hasattr(agent, 'validate_vcf'):
            print("\nTesting direct tool call via agent.validate_vcf...")
            result = agent.validate_vcf(test_vcf)
            print(f"✅ Direct tool call result: {result[:200]}...")
        else:
            print("❌ validate_vcf not available as direct attribute")
            
        # Try direct tool call via agent.tool if it exists
        if hasattr(agent, 'tool'):
            print("\nTesting direct tool call via agent.tool...")
            if hasattr(agent.tool, 'validate_vcf'):
                result = agent.tool.validate_vcf(filepath=test_vcf)
                print(f"✅ agent.tool.validate_vcf result: {str(result)[:200]}...")
            else:
                print("❌ agent.tool.validate_vcf not available")
                tool_attrs = [attr for attr in dir(agent.tool) if not attr.startswith('_')]
                print(f"Available tool methods: {tool_attrs}")
        
        # Test echo tool
        print("\nTesting echo tool...")
        if hasattr(agent, 'echo'):
            echo_result = agent.echo("Hello from direct call!")
            print(f"✅ Echo result: {echo_result}")
        elif hasattr(agent, 'tool') and hasattr(agent.tool, 'echo'):
            echo_result = agent.tool.echo(text="Hello from agent.tool!")
            print(f"✅ agent.tool.echo result: {echo_result}")
        else:
            print("❌ Echo tool not available")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Cleanup
        os.unlink(test_vcf)
        print(f"Cleaned up test file: {test_vcf}")

if __name__ == "__main__":
    test_direct_tool_calling() 