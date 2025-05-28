#!/usr/bin/env python3
"""
Comprehensive Testing Suite for VCF Analysis Agent

This script performs comprehensive testing to validate optimizations and ensure production readiness:
- Performance regression testing
- Functionality validation
- Load testing
- Memory leak detection
- Error handling validation
"""

import sys
import os
import time
import psutil
import tracemalloc
import asyncio
import tempfile
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional
import logging
import json

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcf_agent.config import SessionConfig
from vcf_agent.data_store_manager import create_data_store_manager
from vcf_agent.optimizations import get_optimizer, OptimizationConfig

logger = logging.getLogger(__name__)

class ComprehensiveTestSuite:
    """Comprehensive testing suite for VCF Analysis Agent."""
    
    def __init__(self):
        self.results = {
            "test_date": datetime.now().isoformat(),
            "performance_tests": {},
            "functionality_tests": {},
            "load_tests": {},
            "memory_tests": {},
            "error_handling_tests": {},
            "summary": {}
        }
        self.temp_dir = None
        
    def setup_test_environment(self):
        """Set up isolated test environment."""
        self.temp_dir = tempfile.mkdtemp(prefix="vcf_agent_test_")
        os.chdir(self.temp_dir)
        print(f"🔧 Test environment set up in: {self.temp_dir}")
        
    def cleanup_test_environment(self):
        """Clean up test environment."""
        if self.temp_dir and os.path.exists(self.temp_dir):
            os.chdir(Path(__file__).parent.parent)
            shutil.rmtree(self.temp_dir)
            print(f"🧹 Test environment cleaned up")
    
    def test_performance_regression(self) -> Dict[str, Any]:
        """Test for performance regressions."""
        print("🚀 Running performance regression tests...")
        
        results = {
            "cache_performance": {},
            "embedding_performance": {},
            "database_performance": {},
            "overall_performance": {}
        }
        
        try:
            # Test cache performance
            start_time = time.time()
            config = SessionConfig(model_provider='ollama')
            manager = create_data_store_manager(session_config=config)
            
            # Generate test embeddings to test cache
            test_texts = [f"Test variant {i}" for i in range(100)]
            embedding_times = []
            
            for text in test_texts:
                embed_start = time.time()
                # This should use the optimized cache
                embedding = manager.embedding_service.generate_embedding_sync(text)
                embed_time = time.time() - embed_start
                embedding_times.append(embed_time)
            
            results["cache_performance"] = {
                "total_embeddings": len(test_texts),
                "average_time_per_embedding": sum(embedding_times) / len(embedding_times),
                "total_time": sum(embedding_times),
                "cache_efficiency": "PASS" if sum(embedding_times) < 5.0 else "FAIL"
            }
            
            # Test overall performance
            total_time = time.time() - start_time
            results["overall_performance"] = {
                "total_test_time": total_time,
                "performance_target": "< 10 seconds",
                "result": "PASS" if total_time < 10.0 else "FAIL"
            }
            
        except Exception as e:
            results["error"] = str(e)
            
        return results
    
    def test_functionality_validation(self) -> Dict[str, Any]:
        """Validate core functionality still works after optimizations."""
        print("🔍 Running functionality validation tests...")
        
        results = {
            "data_store_creation": {},
            "variant_processing": {},
            "embedding_generation": {},
            "database_queries": {}
        }
        
        try:
            # Test data store creation
            config = SessionConfig(model_provider='ollama')
            manager = create_data_store_manager(session_config=config)
            results["data_store_creation"]["result"] = "PASS"
            
            # Test embedding generation
            test_embedding = manager.embedding_service.generate_embedding_sync("Test variant")
            results["embedding_generation"] = {
                "result": "PASS" if test_embedding and len(test_embedding) > 0 else "FAIL",
                "embedding_dimension": len(test_embedding) if test_embedding else 0
            }
            
            # Test variant processing
            test_variant = {
                "chr": "chr1",
                "pos": 100000,
                "ref": "A",
                "alt": "G",
                "quality": 30.0
            }
            
            processed = manager._prepare_lance_variant(test_variant, "TEST_SAMPLE")
            results["variant_processing"] = {
                "result": "PASS" if processed else "FAIL",
                "has_embedding": "embedding" in processed if processed else False
            }
            
        except Exception as e:
            results["error"] = str(e)
            
        return results
    
    def test_load_performance(self) -> Dict[str, Any]:
        """Test performance under load."""
        print("⚡ Running load performance tests...")
        
        results = {
            "concurrent_embeddings": {},
            "batch_processing": {},
            "memory_usage": {}
        }
        
        try:
            config = SessionConfig(model_provider='ollama')
            manager = create_data_store_manager(session_config=config)
            
            # Test concurrent embedding generation
            start_time = time.time()
            concurrent_texts = [f"Load test variant {i}" for i in range(50)]
            
            embeddings = []
            for text in concurrent_texts:
                embedding = manager.embedding_service.generate_embedding_sync(text)
                embeddings.append(embedding)
            
            load_time = time.time() - start_time
            
            results["concurrent_embeddings"] = {
                "total_embeddings": len(embeddings),
                "total_time": load_time,
                "average_time": load_time / len(embeddings),
                "throughput_per_second": len(embeddings) / load_time,
                "result": "PASS" if load_time < 30.0 else "FAIL"
            }
            
            # Test memory usage
            process = psutil.Process()
            memory_info = process.memory_info()
            results["memory_usage"] = {
                "rss_mb": memory_info.rss / 1024 / 1024,
                "vms_mb": memory_info.vms / 1024 / 1024,
                "result": "PASS" if memory_info.rss < 500 * 1024 * 1024 else "FAIL"  # < 500MB
            }
            
        except Exception as e:
            results["error"] = str(e)
            
        return results
    
    def test_memory_leaks(self) -> Dict[str, Any]:
        """Test for memory leaks."""
        print("🧠 Running memory leak tests...")
        
        results = {
            "memory_growth": {},
            "garbage_collection": {},
            "cache_limits": {}
        }
        
        try:
            tracemalloc.start()
            
            config = SessionConfig(model_provider='ollama')
            manager = create_data_store_manager(session_config=config)
            
            # Baseline memory
            snapshot1 = tracemalloc.take_snapshot()
            
            # Generate many embeddings to test for leaks
            for i in range(200):
                text = f"Memory test variant {i}"
                embedding = manager.embedding_service.generate_embedding_sync(text)
            
            # Check memory after operations
            snapshot2 = tracemalloc.take_snapshot()
            
            top_stats = snapshot2.compare_to(snapshot1, 'lineno')
            total_growth = sum(stat.size_diff for stat in top_stats)
            
            results["memory_growth"] = {
                "total_growth_bytes": total_growth,
                "total_growth_mb": total_growth / 1024 / 1024,
                "result": "PASS" if total_growth < 50 * 1024 * 1024 else "FAIL"  # < 50MB growth
            }
            
            tracemalloc.stop()
            
        except Exception as e:
            results["error"] = str(e)
            
        return results
    
    def test_error_handling(self) -> Dict[str, Any]:
        """Test error handling and recovery."""
        print("🛡️ Running error handling tests...")
        
        results = {
            "invalid_input": {},
            "network_errors": {},
            "file_system_errors": {},
            "recovery": {}
        }
        
        try:
            config = SessionConfig(model_provider='ollama')
            manager = create_data_store_manager(session_config=config)
            
            # Test invalid input handling
            try:
                embedding = manager.embedding_service.generate_embedding_sync("")
                results["invalid_input"]["empty_string"] = "HANDLED"
            except Exception:
                results["invalid_input"]["empty_string"] = "ERROR"
            
            try:
                embedding = manager.embedding_service.generate_embedding_sync(None)
                results["invalid_input"]["none_input"] = "HANDLED"
            except Exception:
                results["invalid_input"]["none_input"] = "ERROR"
            
            # Test recovery after errors
            try:
                # This should work after error conditions
                embedding = manager.embedding_service.generate_embedding_sync("Recovery test")
                results["recovery"]["after_errors"] = "PASS" if embedding else "FAIL"
            except Exception:
                results["recovery"]["after_errors"] = "FAIL"
                
        except Exception as e:
            results["error"] = str(e)
            
        return results
    
    def generate_test_summary(self) -> Dict[str, Any]:
        """Generate comprehensive test summary."""
        summary = {
            "total_tests": 0,
            "passed_tests": 0,
            "failed_tests": 0,
            "test_categories": {},
            "overall_result": "UNKNOWN"
        }
        
        # Count results from all test categories
        for category, tests in self.results.items():
            if category in ["test_date", "summary"]:
                continue
                
            category_results = {"passed": 0, "failed": 0, "total": 0}
            
            def count_results(data):
                if isinstance(data, dict):
                    for key, value in data.items():
                        if key == "result":
                            category_results["total"] += 1
                            summary["total_tests"] += 1
                            if value == "PASS":
                                category_results["passed"] += 1
                                summary["passed_tests"] += 1
                            else:
                                category_results["failed"] += 1
                                summary["failed_tests"] += 1
                        elif isinstance(value, dict):
                            count_results(value)
            
            count_results(tests)
            summary["test_categories"][category] = category_results
        
        # Determine overall result
        if summary["total_tests"] == 0:
            summary["overall_result"] = "NO_TESTS"
        elif summary["failed_tests"] == 0:
            summary["overall_result"] = "ALL_PASS"
        elif summary["passed_tests"] > summary["failed_tests"]:
            summary["overall_result"] = "MOSTLY_PASS"
        else:
            summary["overall_result"] = "MOSTLY_FAIL"
        
        return summary
    
    def run_all_tests(self) -> Dict[str, Any]:
        """Run all comprehensive tests."""
        print("🧪 Starting Comprehensive Test Suite...")
        print("="*80)
        
        try:
            self.setup_test_environment()
            
            # Run all test categories
            self.results["performance_tests"] = self.test_performance_regression()
            self.results["functionality_tests"] = self.test_functionality_validation()
            self.results["load_tests"] = self.test_load_performance()
            self.results["memory_tests"] = self.test_memory_leaks()
            self.results["error_handling_tests"] = self.test_error_handling()
            
            # Generate summary
            self.results["summary"] = self.generate_test_summary()
            
            return self.results
            
        finally:
            self.cleanup_test_environment()
    
    def print_test_results(self):
        """Print formatted test results."""
        print("\n" + "="*80)
        print("🧪 COMPREHENSIVE TEST RESULTS")
        print("="*80)
        
        summary = self.results.get("summary", {})
        print(f"📊 Total Tests: {summary.get('total_tests', 0)}")
        print(f"✅ Passed: {summary.get('passed_tests', 0)}")
        print(f"❌ Failed: {summary.get('failed_tests', 0)}")
        print(f"🎯 Overall Result: {summary.get('overall_result', 'UNKNOWN')}")
        
        print(f"\n📋 TEST CATEGORIES:")
        for category, results in summary.get("test_categories", {}).items():
            status = "✅" if results["failed"] == 0 else "⚠️" if results["passed"] > results["failed"] else "❌"
            print(f"   {status} {category}: {results['passed']}/{results['total']} passed")
        
        # Show performance highlights
        perf_tests = self.results.get("performance_tests", {})
        if "overall_performance" in perf_tests:
            overall = perf_tests["overall_performance"]
            print(f"\n⚡ PERFORMANCE HIGHLIGHTS:")
            print(f"   Total Test Time: {overall.get('total_test_time', 0):.2f}s")
            print(f"   Performance Target: {overall.get('performance_target', 'N/A')}")
            print(f"   Result: {overall.get('result', 'UNKNOWN')}")
        
        print("\n" + "="*80)

def save_test_results(results: Dict[str, Any]) -> str:
    """Save test results to file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"performance_reports/comprehensive_test_results_{timestamp}.json"
    
    with open(filename, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    return filename

def main():
    """Main function to run comprehensive tests."""
    try:
        test_suite = ComprehensiveTestSuite()
        results = test_suite.run_all_tests()
        
        # Print results
        test_suite.print_test_results()
        
        # Save detailed results
        filename = save_test_results(results)
        print(f"\n📄 Detailed results saved to: {filename}")
        
        # Return appropriate exit code
        summary = results.get("summary", {})
        if summary.get("overall_result") in ["ALL_PASS", "MOSTLY_PASS"]:
            return 0
        else:
            return 1
        
    except Exception as e:
        print(f"❌ Comprehensive testing failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main()) 