#!/usr/bin/env python3
"""Test automatic tool execution functionality."""

import sys
import tempfile
import os
from pathlib import Path

# Add the src directory to the Python path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from vcf_agent.agent import get_agent_with_session
from vcf_agent.config import SessionConfig

def test_auto_execution():
    """Test automatic tool execution."""
    print("=== Testing Automatic Tool Execution ===")
    
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
        
        # Test 1: Direct tool call (should work)
        print("\n--- Test 1: Direct Tool Call ---")
        result = agent.validate_vcf(test_vcf)
        print(f"✅ Direct call result: {result[:100]}...")
        
        # Test 2: Natural language with explicit tool mention
        print("\n--- Test 2: Explicit Tool Request ---")
        response = agent(f"Use the validate_vcf tool to check this file: {test_vcf}")
        print(f"Response type: {type(response)}")
        print(f"Response: {str(response)[:300]}...")
        
        # Test 3: Natural language without explicit tool mention
        print("\n--- Test 3: Implicit Tool Request ---")
        response2 = agent(f"Please validate this VCF file: {test_vcf}")
        print(f"Response type: {type(response2)}")
        print(f"Response: {str(response2)[:300]}...")
        
        # Test 4: Check if agent has tool execution capabilities
        print("\n--- Test 4: Agent Capabilities ---")
        print(f"Agent tool names: {agent.tool_names}")
        print(f"Agent has tool_caller: {hasattr(agent, 'tool_caller')}")
        print(f"Agent has ToolCaller: {hasattr(agent, 'ToolCaller')}")
        
        # Test 5: Try echo tool for comparison
        print("\n--- Test 5: Echo Tool Test ---")
        echo_response = agent("Use the echo tool to say 'Hello World'")
        print(f"Echo response: {str(echo_response)[:200]}...")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Cleanup
        os.unlink(test_vcf)
        print(f"Cleaned up test file: {test_vcf}")

if __name__ == "__main__":
    test_auto_execution() 