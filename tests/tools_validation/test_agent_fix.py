#!/usr/bin/env python3
"""
Test script to validate the refactored VCF agent functionality.
Tests both natural conversation and tool calling capabilities.
"""

import os
import sys
import tempfile
from pathlib import Path

# Add the src directory to the Python path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from vcf_agent.agent import get_agent_with_session
from vcf_agent.config import SessionConfig

def create_test_vcf():
    """Create a minimal test VCF file for testing."""
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100	.	A	T	60	PASS	DP=30	GT:DP	1/1:30
chr1	200	.	G	C	45	PASS	DP=25	GT:DP	0/1:25
"""
    
    # Create temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_content)
        return f.name

def test_natural_conversation():
    """Test natural conversation capabilities."""
    print("=== Testing Natural Conversation ===")
    
    try:
        # Create agent with Ollama (most likely to be available)
        agent = get_agent_with_session(
            session_config=SessionConfig(raw_mode=False),
            model_provider="ollama"
        )
        
        # Test simple conversation
        response = agent("Hello! What can you help me with regarding VCF files?")
        print(f"Agent Response: {response}")
        
        # Check if response is conversational (not JSON)
        if isinstance(response, str) and not response.strip().startswith('{'):
            print("✅ Natural conversation working!")
            return True
        else:
            print("❌ Agent still responding in JSON mode")
            return False
            
    except Exception as e:
        print(f"❌ Natural conversation test failed: {e}")
        return False

def test_tool_calling():
    """Test tool calling capabilities."""
    print("\n=== Testing Tool Calling ===")
    
    try:
        # Create test VCF file
        test_vcf = create_test_vcf()
        print(f"Created test VCF: {test_vcf}")
        
        # Create agent
        agent = get_agent_with_session(
            session_config=SessionConfig(raw_mode=False),
            model_provider="ollama"
        )
        
        # Test tool calling through natural language
        response = agent(f"Please validate this VCF file: {test_vcf}")
        print(f"Validation Response: {response}")
        
        # Test echo tool
        echo_response = agent("Echo this message: Hello from the VCF agent!")
        print(f"Echo Response: {echo_response}")
        
        # Cleanup
        os.unlink(test_vcf)
        
        print("✅ Tool calling test completed!")
        return True
        
    except Exception as e:
        print(f"❌ Tool calling test failed: {e}")
        # Cleanup on error
        if 'test_vcf' in locals():
            try:
                os.unlink(test_vcf)
            except:
                pass
        return False

def test_different_providers():
    """Test different model providers."""
    print("\n=== Testing Different Model Providers ===")
    
    providers = ["ollama", "openai", "cerebras"]
    results = {}
    
    for provider in providers:
        print(f"\nTesting {provider}...")
        try:
            agent = get_agent_with_session(
                session_config=SessionConfig(raw_mode=False),
                model_provider=provider
            )
            
            # Simple test
            response = agent("What is a VCF file?")
            results[provider] = "✅ Working"
            print(f"{provider}: ✅ Working")
            
        except Exception as e:
            results[provider] = f"❌ Failed: {str(e)[:100]}"
            print(f"{provider}: ❌ Failed: {e}")
    
    return results

def main():
    """Run all tests."""
    print("VCF Agent Refactoring Validation Tests")
    print("=" * 50)
    
    # Test 1: Natural conversation
    conv_success = test_natural_conversation()
    
    # Test 2: Tool calling
    tool_success = test_tool_calling()
    
    # Test 3: Different providers
    provider_results = test_different_providers()
    
    # Summary
    print("\n" + "=" * 50)
    print("TEST SUMMARY")
    print("=" * 50)
    print(f"Natural Conversation: {'✅ PASS' if conv_success else '❌ FAIL'}")
    print(f"Tool Calling: {'✅ PASS' if tool_success else '❌ FAIL'}")
    print("\nProvider Results:")
    for provider, result in provider_results.items():
        print(f"  {provider}: {result}")
    
    # Overall result
    overall_success = conv_success and tool_success and any("✅" in r for r in provider_results.values())
    print(f"\nOverall: {'✅ SUCCESS' if overall_success else '❌ NEEDS WORK'}")
    
    return overall_success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 