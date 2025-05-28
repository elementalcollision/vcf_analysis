#!/usr/bin/env python3
"""
Simple script to run VCF Analysis Agent load tests.

This script runs our comprehensive load testing suite to measure:
- Batch processing performance
- Concurrent analysis capabilities  
- Database performance under load
- Memory usage patterns

Target: 1,000+ variants/second processing capability
"""

import sys
import os
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from tests.performance.test_load_performance import LoadTestConfig, LoadTestRunner

def main():
    """Run load testing suite."""
    print("🚀 VCF Analysis Agent - Load Testing Suite")
    print("=" * 60)
    
    # Configure load tests
    config = LoadTestConfig(
        num_variants=500,  # Start with smaller load
        num_concurrent_users=3,
        batch_size=50,
        target_throughput=500.0,  # 500 variants/sec target
        memory_limit_mb=1024,
        enable_profiling=True
    )
    
    print(f"📋 Test Configuration:")
    print(f"   • Variants: {config.num_variants}")
    print(f"   • Concurrent Users: {config.num_concurrent_users}")
    print(f"   • Batch Size: {config.batch_size}")
    print(f"   • Target Throughput: {config.target_throughput} variants/sec")
    print(f"   • Memory Limit: {config.memory_limit_mb}MB")
    print()
    
    # Initialize runner
    runner = LoadTestRunner(config)
    
    try:
        # Run individual tests
        print("🧪 Running individual load tests...")
        
        # 1. Batch Processing Test
        print("\n1️⃣ Batch Processing Test")
        batch_metrics = runner.run_batch_processing_test()
        
        # 2. Database Load Test  
        print("\n2️⃣ Database Load Test")
        db_metrics = runner.run_database_load_test()
        
        # 3. Memory Stress Test
        print("\n3️⃣ Memory Stress Test")
        memory_metrics = runner.run_memory_stress_test()
        
        # Print summary
        print("\n" + "=" * 60)
        print("📊 LOAD TESTING SUMMARY")
        print("=" * 60)
        
        runner._print_summary()
        
        # Save results
        runner.save_results("performance_reports/load_test_results.json")
        
        print("\n✅ Load testing completed successfully!")
        
        # Check if targets were met
        batch_target_met = batch_metrics.additional_metrics.get("throughput_ratio", 0) >= 1.0
        db_target_met = db_metrics.additional_metrics.get("search_time_target_met", False)
        memory_ok = not memory_metrics.additional_metrics.get("memory_limit_exceeded", True)
        
        if batch_target_met and db_target_met and memory_ok:
            print("🎯 All performance targets MET! System is production-ready.")
            return 0
        else:
            print("⚠️ Some performance targets NOT MET. Review results for optimization opportunities.")
            return 1
            
    except Exception as e:
        print(f"\n❌ Load testing failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main()) 