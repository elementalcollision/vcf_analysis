#!/usr/bin/env python3
"""
Demo validation script for prompt contracts.
Tests all three prompt contracts with sample data to ensure they work for the client demo.
"""

import os
import sys
import json
import yaml
from pathlib import Path
from jsonschema import validate, ValidationError

# Add src to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'src')))

from vcf_agent.agent import run_llm_analysis_task, get_agent_with_session
from vcf_agent.config import SessionConfig
from vcf_agent.prompt_templates import load_prompt_contract, get_prompt_for_task

def test_prompt_contract_integration():
    """Test that prompt contracts work end-to-end with the agent system."""
    
    print("🧪 Testing Prompt Contracts for Demo Readiness")
    print("=" * 60)
    
    # Test data
    test_vcf = "sample_data/small_valid.vcf"
    test_reference = "sample_data/minimal_reference.fa"
    
    if not os.path.exists(test_vcf):
        print(f"❌ Test VCF file not found: {test_vcf}")
        print("Creating minimal test VCF...")
        create_minimal_test_vcf(test_vcf)
    
    if not os.path.exists(test_reference):
        print(f"❌ Test reference file not found: {test_reference}")
        print("Creating minimal reference file...")
        create_minimal_reference(test_reference)
    
    # Test all three contracts
    contracts_to_test = [
        ("vcf_analysis_summary_v1", [test_vcf]),
        ("vcf_summarization_v1", [test_vcf]),
        ("vcf_comparison_v1", [test_vcf, test_vcf])  # Compare file with itself for testing
    ]
    
    results = {}
    
    for contract_name, file_paths in contracts_to_test:
        print(f"\n📋 Testing Contract: {contract_name}")
        print("-" * 40)
        
        try:
            # Test 1: Load contract
            print("1. Loading contract...")
            contract = load_prompt_contract(contract_name)
            print(f"   ✅ Contract loaded: {contract['id']} v{contract['version']}")
            
            # Test 2: Generate prompt
            print("2. Generating prompt...")
            prompt = get_prompt_for_task(contract_name, file_paths)
            print(f"   ✅ Prompt generated ({len(prompt)} chars)")
            
            # Test 3: Test with Ollama (local model)
            print("3. Testing with Ollama...")
            try:
                config = SessionConfig(raw_mode=False, model_provider="ollama")
                result = run_llm_analysis_task(
                    task=contract_name,
                    file_paths=file_paths,
                    model_provider="ollama"
                )
                
                # Parse and validate JSON
                result_data = json.loads(result)
                
                # Validate against schema
                validate(instance=result_data, schema=contract['json_schema'])
                
                print(f"   ✅ Ollama test passed")
                print(f"   📊 Result keys: {list(result_data.keys())}")
                
                results[f"{contract_name}_ollama"] = {
                    "status": "success",
                    "result": result_data,
                    "schema_valid": True
                }
                
            except Exception as e:
                print(f"   ⚠️  Ollama test failed: {e}")
                results[f"{contract_name}_ollama"] = {
                    "status": "error",
                    "error": str(e)
                }
            
            # Test 4: Test direct agent calling (demo style)
            print("4. Testing direct agent calling...")
            try:
                agent = get_agent_with_session(SessionConfig(raw_mode=False), "ollama")
                
                # Test natural language that should trigger the tool
                if contract_name == "vcf_analysis_summary_v1":
                    response = agent(f"Please analyze the VCF file at {file_paths[0]} and provide a comprehensive summary")
                elif contract_name == "vcf_summarization_v1":
                    response = agent(f"Summarize the variants in {file_paths[0]}")
                elif contract_name == "vcf_comparison_v1":
                    response = agent(f"Compare the VCF files {file_paths[0]} and {file_paths[1]} using reference {test_reference}")
                
                print(f"   ✅ Direct agent test passed")
                print(f"   📝 Response length: {len(str(response))} chars")
                
                results[f"{contract_name}_direct"] = {
                    "status": "success",
                    "response": str(response)[:200] + "..." if len(str(response)) > 200 else str(response)
                }
                
            except Exception as e:
                print(f"   ⚠️  Direct agent test failed: {e}")
                results[f"{contract_name}_direct"] = {
                    "status": "error",
                    "error": str(e)
                }
                
        except Exception as e:
            print(f"   ❌ Contract test failed: {e}")
            results[contract_name] = {
                "status": "error",
                "error": str(e)
            }
    
    # Summary
    print("\n" + "=" * 60)
    print("📊 DEMO READINESS SUMMARY")
    print("=" * 60)
    
    success_count = sum(1 for r in results.values() if r.get("status") == "success")
    total_count = len(results)
    
    print(f"✅ Successful tests: {success_count}/{total_count}")
    
    if success_count == total_count:
        print("🎉 ALL TESTS PASSED - DEMO READY!")
    else:
        print("⚠️  Some tests failed - needs attention before demo")
        
        print("\n❌ Failed tests:")
        for test_name, result in results.items():
            if result.get("status") == "error":
                print(f"   - {test_name}: {result.get('error', 'Unknown error')}")
    
    # Save results for reference
    with open("prompt_contracts_demo_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n📄 Detailed results saved to: prompt_contracts_demo_results.json")
    
    return success_count == total_count

def create_minimal_test_vcf(filepath):
    """Create a minimal valid VCF file for testing."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    vcf_content = """##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100	.	A	G	99	PASS	DP=30	GT:DP	0/1:25
chr1	200	.	T	C	95	PASS	DP=35	GT:DP	1/1:30
chr2	300	.	G	A	88	PASS	DP=28	GT:DP	0/1:22
"""
    
    with open(filepath, 'w') as f:
        f.write(vcf_content)
    
    print(f"✅ Created test VCF: {filepath}")

def create_minimal_reference(filepath):
    """Create a minimal reference genome file for testing."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    ref_content = """>chr1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>chr2
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
"""
    
    with open(filepath, 'w') as f:
        f.write(ref_content)
    
    print(f"✅ Created test reference: {filepath}")

def test_reproducibility():
    """Test that prompt contracts produce reproducible results."""
    print("\n🔄 Testing Reproducibility")
    print("-" * 40)
    
    test_vcf = "sample_data/small_valid.vcf"
    contract_name = "vcf_analysis_summary_v1"
    
    results = []
    
    # Run the same analysis multiple times
    for i in range(3):
        print(f"Run {i+1}/3...")
        try:
            result = run_llm_analysis_task(
                task=contract_name,
                file_paths=[test_vcf],
                model_provider="ollama"
            )
            result_data = json.loads(result)
            results.append(result_data)
            print(f"   ✅ Run {i+1} completed")
        except Exception as e:
            print(f"   ❌ Run {i+1} failed: {e}")
            return False
    
    # Check if results are similar (allowing for some variation in AI responses)
    if len(results) >= 2:
        # Compare key fields that should be deterministic
        first_result = results[0]
        consistent_fields = []
        
        for key in ["variant_count", "variant_types"]:
            if key in first_result:
                values = [r.get(key) for r in results if key in r]
                if len(set(str(v) for v in values)) == 1:
                    consistent_fields.append(key)
        
        print(f"   📊 Consistent fields: {consistent_fields}")
        
        if len(consistent_fields) >= 1:
            print("   ✅ Reproducibility test passed (core metrics consistent)")
            return True
        else:
            print("   ⚠️  Reproducibility test failed (no consistent fields)")
            return False
    
    return False

if __name__ == "__main__":
    print("🚀 VCF Agent Prompt Contracts Demo Validation")
    print("=" * 60)
    
    # Test integration
    integration_success = test_prompt_contract_integration()
    
    # Test reproducibility
    reproducibility_success = test_reproducibility()
    
    print("\n" + "=" * 60)
    print("🎯 FINAL DEMO READINESS STATUS")
    print("=" * 60)
    
    if integration_success and reproducibility_success:
        print("🎉 DEMO READY - All prompt contracts working!")
        print("✅ Integration tests passed")
        print("✅ Reproducibility tests passed")
        sys.exit(0)
    else:
        print("❌ DEMO NOT READY - Issues found:")
        if not integration_success:
            print("   - Integration tests failed")
        if not reproducibility_success:
            print("   - Reproducibility tests failed")
        sys.exit(1) 