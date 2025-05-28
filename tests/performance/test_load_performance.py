"""
Comprehensive Load Testing and Performance Profiling for VCF Analysis Agent

This module implements load testing scenarios to validate system performance
under realistic workloads and identify bottlenecks.

Target Performance Goals:
- 1,000+ variants/second processing capability
- <100ms vector search response time
- <500ms graph query response time
- Support for concurrent users (10+ simultaneous sessions)
"""

import asyncio
import time
import threading
import statistics
import psutil
import pytest
import json
import tempfile
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Tuple, Optional
from dataclasses import dataclass, asdict
from pathlib import Path
import pandas as pd
import sys
import random

# Add src to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

# VCF Agent imports
from vcf_agent.agent import get_agent_with_session
from vcf_agent.config import SessionConfig
from vcf_agent.data_store_manager import create_data_store_manager
from vcf_agent.optimizations import PerformanceOptimizer, OptimizationConfig
from vcf_agent.vcf_ingestion import VCFIngestionPipeline, IngestionConfig
from vcf_agent.lancedb_integration import get_db, get_or_create_table, add_variants


@dataclass
class PerformanceMetrics:
    """Performance metrics for load testing."""
    test_name: str
    start_time: float
    end_time: float
    duration: float
    throughput: float  # operations per second
    memory_usage_mb: float
    cpu_usage_percent: float
    success_count: int
    error_count: int
    response_times: List[float]
    avg_response_time: float
    p95_response_time: float
    p99_response_time: float
    additional_metrics: Dict[str, Any]


@dataclass
class LoadTestConfig:
    """Configuration for load testing scenarios."""
    num_variants: int = 1000
    num_concurrent_users: int = 10
    test_duration_seconds: int = 60
    batch_size: int = 100
    target_throughput: float = 1000.0  # variants/second
    memory_limit_mb: int = 2048
    enable_profiling: bool = True


class PerformanceMonitor:
    """Monitor system performance during load tests."""
    
    def __init__(self):
        self.process = psutil.Process()
        self.monitoring = False
        self.metrics = []
        self._monitor_thread = None
    
    def start_monitoring(self, interval: float = 0.1):
        """Start performance monitoring."""
        self.monitoring = True
        self.metrics = []
        self._monitor_thread = threading.Thread(target=self._monitor_loop, args=(interval,))
        self._monitor_thread.start()
    
    def stop_monitoring(self) -> Dict[str, Any]:
        """Stop monitoring and return aggregated metrics."""
        self.monitoring = False
        if self._monitor_thread:
            self._monitor_thread.join()
        
        if not self.metrics:
            return {"error": "No metrics collected"}
        
        cpu_values = [m['cpu_percent'] for m in self.metrics]
        memory_values = [m['memory_mb'] for m in self.metrics]
        
        return {
            "cpu_percent": {
                "avg": statistics.mean(cpu_values),
                "max": max(cpu_values),
                "min": min(cpu_values),
                "p95": statistics.quantiles(cpu_values, n=20)[18] if len(cpu_values) > 20 else max(cpu_values)
            },
            "memory_mb": {
                "avg": statistics.mean(memory_values),
                "max": max(memory_values),
                "min": min(memory_values),
                "p95": statistics.quantiles(memory_values, n=20)[18] if len(memory_values) > 20 else max(memory_values)
            },
            "sample_count": len(self.metrics)
        }
    
    def _monitor_loop(self, interval: float):
        """Monitor loop running in separate thread."""
        while self.monitoring:
            try:
                cpu_percent = self.process.cpu_percent()
                memory_info = self.process.memory_info()
                memory_mb = memory_info.rss / 1024 / 1024
                
                self.metrics.append({
                    "timestamp": time.time(),
                    "cpu_percent": cpu_percent,
                    "memory_mb": memory_mb
                })
                
                time.sleep(interval)
            except Exception as e:
                print(f"Monitoring error: {e}")
                break


class VCFDataGenerator:
    """Generate synthetic VCF data for load testing."""
    
    @staticmethod
    def generate_variant_data(num_variants: int, sample_prefix: str = "test") -> List[Dict[str, Any]]:
        """Generate synthetic variant data for load testing."""
        # Add timestamp and random component to ensure uniqueness across test runs
        timestamp = int(time.time() * 1000)  # milliseconds
        random_suffix = random.randint(1000, 9999)
        
        variants = []
        chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
        bases = ['A', 'T', 'G', 'C']
        
        for i in range(num_variants):
            # Create unique variant ID with timestamp and random suffix
            variant_id = f"{sample_prefix}_variant_{timestamp}_{random_suffix}_{i}"
            
            variant = {
                'id': variant_id,
                'chr': random.choice(chromosomes),
                'pos': random.randint(1000, 100000000),
                'ref': random.choice(bases),
                'alt': random.choice(bases),
                'variant_type': random.choice(['SNV', 'INDEL', 'CNV']),
                'quality_score': random.uniform(10.0, 99.0),
                'filter_status': random.choice(['PASS', 'FAIL', 'UNKNOWN']),
                'allele_frequency': random.uniform(0.001, 0.5)
            }
            variants.append(variant)
        
        return variants
    
    @staticmethod
    def create_test_vcf_file(num_variants: int, output_path: str) -> str:
        """Create a test VCF file with specified number of variants."""
        header = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
"""
        
        chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
        
        with open(output_path, 'w') as f:
            f.write(header)
            
            for i in range(num_variants):
                chrom = chromosomes[i % len(chromosomes)]
                pos = 1000000 + (i * 100)
                variant_id = f"rs{i}"
                ref = "A"
                alt = "G"
                qual = 30.0 + (i % 50)
                filter_val = "PASS"
                info = f"DP={50 + (i % 100)};AF={0.1 + (i % 10) * 0.1}"
                format_val = "GT:DP"
                sample = f"{'0/1' if i % 2 == 0 else '1/1'}:{50 + (i % 100)}"
                
                line = f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\t{format_val}\t{sample}\n"
                f.write(line)
        
        return output_path


class LoadTestRunner:
    """Main load testing runner with comprehensive scenarios."""
    
    def __init__(self, config: LoadTestConfig):
        self.config = config
        self.monitor = PerformanceMonitor()
        self.results: List[PerformanceMetrics] = []
        
        # Initialize performance optimizer
        opt_config = OptimizationConfig(
            enable_embedding_cache=True,
            enable_query_batching=True,
            enable_async_processing=True,
            max_concurrent_tasks=config.num_concurrent_users
        )
        self.optimizer = PerformanceOptimizer(opt_config)
    
    def run_batch_processing_test(self) -> PerformanceMetrics:
        """Test batch VCF processing performance."""
        print(f"🚀 Running batch processing test with {self.config.num_variants} variants...")
        
        # Generate test data
        variants_data = VCFDataGenerator.generate_variant_data(self.config.num_variants)
        
        # Create temporary VCF file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        VCFDataGenerator.create_test_vcf_file(self.config.num_variants, vcf_path)
        
        try:
            # Start monitoring
            self.monitor.start_monitoring()
            start_time = time.time()
            
            # Initialize data store manager
            manager = create_data_store_manager()
            
            # Process variants in batches
            batch_size = self.config.batch_size
            response_times = []
            success_count = 0
            error_count = 0
            
            for i in range(0, len(variants_data), batch_size):
                batch_start = time.time()
                batch = variants_data[i:i + batch_size]
                
                try:
                    # Add batch to data stores
                    result = manager.add_sample_with_variants(
                        sample_data={"id": f"LOAD_TEST_SAMPLE_{i}", "name": f"Load Test Sample {i}"},
                        variants_data=batch
                    )
                    
                    if result.get("success"):
                        success_count += len(batch)
                    else:
                        error_count += len(batch)
                    
                    batch_time = time.time() - batch_start
                    response_times.append(batch_time)
                    
                except Exception as e:
                    print(f"Batch processing error: {e}")
                    error_count += len(batch)
                    response_times.append(time.time() - batch_start)
            
            end_time = time.time()
            duration = end_time - start_time
            
            # Stop monitoring
            system_metrics = self.monitor.stop_monitoring()
            
            # Calculate metrics
            throughput = success_count / duration if duration > 0 else 0
            avg_response_time = statistics.mean(response_times) if response_times else 0
            p95_response_time = statistics.quantiles(response_times, n=20)[18] if len(response_times) > 20 else max(response_times) if response_times else 0
            p99_response_time = statistics.quantiles(response_times, n=100)[98] if len(response_times) > 100 else max(response_times) if response_times else 0
            
            metrics = PerformanceMetrics(
                test_name="batch_processing",
                start_time=start_time,
                end_time=end_time,
                duration=duration,
                throughput=throughput,
                memory_usage_mb=system_metrics.get("memory_mb", {}).get("max", 0),
                cpu_usage_percent=system_metrics.get("cpu_percent", {}).get("avg", 0),
                success_count=success_count,
                error_count=error_count,
                response_times=response_times,
                avg_response_time=avg_response_time,
                p95_response_time=p95_response_time,
                p99_response_time=p99_response_time,
                additional_metrics={
                    "variants_per_second": throughput,
                    "target_throughput": self.config.target_throughput,
                    "throughput_ratio": throughput / self.config.target_throughput,
                    "system_metrics": system_metrics
                }
            )
            
            self.results.append(metrics)
            return metrics
            
        finally:
            # Cleanup
            try:
                os.unlink(vcf_path)
            except:
                pass
    
    def run_concurrent_analysis_test(self) -> PerformanceMetrics:
        """Test concurrent AI analysis performance."""
        print(f"🤖 Running concurrent analysis test with {self.config.num_concurrent_users} users...")
        
        # Create test VCF file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        VCFDataGenerator.create_test_vcf_file(100, vcf_path)  # Smaller file for AI analysis
        
        try:
            # Start monitoring
            self.monitor.start_monitoring()
            start_time = time.time()
            
            response_times = []
            success_count = 0
            error_count = 0
            
            def run_analysis_session(session_id: int) -> Tuple[bool, float]:
                """Run a single analysis session."""
                session_start = time.time()
                try:
                    # Create agent for this session
                    config = SessionConfig(raw_mode=False)
                    agent = get_agent_with_session(config, "ollama")
                    
                    # Run analysis
                    result = agent.vcf_analysis_summary_tool(vcf_path)
                    
                    session_time = time.time() - session_start
                    return True, session_time
                    
                except Exception as e:
                    print(f"Analysis session {session_id} error: {e}")
                    session_time = time.time() - session_start
                    return False, session_time
            
            # Run concurrent sessions
            with ThreadPoolExecutor(max_workers=self.config.num_concurrent_users) as executor:
                futures = [
                    executor.submit(run_analysis_session, i) 
                    for i in range(self.config.num_concurrent_users)
                ]
                
                for future in as_completed(futures):
                    success, response_time = future.result()
                    response_times.append(response_time)
                    if success:
                        success_count += 1
                    else:
                        error_count += 1
            
            end_time = time.time()
            duration = end_time - start_time
            
            # Stop monitoring
            system_metrics = self.monitor.stop_monitoring()
            
            # Calculate metrics
            throughput = success_count / duration if duration > 0 else 0
            avg_response_time = statistics.mean(response_times) if response_times else 0
            p95_response_time = statistics.quantiles(response_times, n=20)[18] if len(response_times) > 20 else max(response_times) if response_times else 0
            p99_response_time = statistics.quantiles(response_times, n=100)[98] if len(response_times) > 100 else max(response_times) if response_times else 0
            
            metrics = PerformanceMetrics(
                test_name="concurrent_analysis",
                start_time=start_time,
                end_time=end_time,
                duration=duration,
                throughput=throughput,
                memory_usage_mb=system_metrics.get("memory_mb", {}).get("max", 0),
                cpu_usage_percent=system_metrics.get("cpu_percent", {}).get("avg", 0),
                success_count=success_count,
                error_count=error_count,
                response_times=response_times,
                avg_response_time=avg_response_time,
                p95_response_time=p95_response_time,
                p99_response_time=p99_response_time,
                additional_metrics={
                    "concurrent_users": self.config.num_concurrent_users,
                    "analyses_per_second": throughput,
                    "system_metrics": system_metrics
                }
            )
            
            self.results.append(metrics)
            return metrics
            
        finally:
            # Cleanup
            try:
                os.unlink(vcf_path)
            except:
                pass
    
    def run_database_load_test(self) -> PerformanceMetrics:
        """Test database performance under load."""
        print("🗄️ Running database load test...")
        
        # Generate test data
        variants_data = VCFDataGenerator.generate_variant_data(self.config.num_variants)
        
        # Start monitoring
        self.monitor.start_monitoring()
        start_time = time.time()
        
        response_times = []
        success_count = 0
        error_count = 0
        
        try:
            # Initialize data store manager
            manager = create_data_store_manager()
            
            # Test vector search performance
            search_queries = [
                "pathogenic variant",
                "BRCA1 mutation",
                "missense variant",
                "high quality variant",
                "chromosome 1 variant"
            ]
            
            # Add some data first
            sample_data = {"id": "DB_LOAD_TEST", "name": "Database Load Test"}
            manager.add_sample_with_variants(sample_data, variants_data[:100])
            
            # Run concurrent searches
            def run_search(query: str) -> Tuple[bool, float]:
                search_start = time.time()
                try:
                    results = manager.search_variants(query, limit=10)
                    search_time = time.time() - search_start
                    return True, search_time
                except Exception as e:
                    print(f"Search error: {e}")
                    search_time = time.time() - search_start
                    return False, search_time
            
            # Run multiple rounds of searches
            with ThreadPoolExecutor(max_workers=5) as executor:
                for round_num in range(10):  # 10 rounds of searches
                    futures = [
                        executor.submit(run_search, query) 
                        for query in search_queries
                    ]
                    
                    for future in as_completed(futures):
                        success, response_time = future.result()
                        response_times.append(response_time)
                        if success:
                            success_count += 1
                        else:
                            error_count += 1
            
            end_time = time.time()
            duration = end_time - start_time
            
            # Stop monitoring
            system_metrics = self.monitor.stop_monitoring()
            
            # Calculate metrics
            throughput = success_count / duration if duration > 0 else 0
            avg_response_time = statistics.mean(response_times) if response_times else 0
            p95_response_time = statistics.quantiles(response_times, n=20)[18] if len(response_times) > 20 else max(response_times) if response_times else 0
            p99_response_time = statistics.quantiles(response_times, n=100)[98] if len(response_times) > 100 else max(response_times) if response_times else 0
            
            metrics = PerformanceMetrics(
                test_name="database_load",
                start_time=start_time,
                end_time=end_time,
                duration=duration,
                throughput=throughput,
                memory_usage_mb=system_metrics.get("memory_mb", {}).get("max", 0),
                cpu_usage_percent=system_metrics.get("cpu_percent", {}).get("avg", 0),
                success_count=success_count,
                error_count=error_count,
                response_times=response_times,
                avg_response_time=avg_response_time,
                p95_response_time=p95_response_time,
                p99_response_time=p99_response_time,
                additional_metrics={
                    "searches_per_second": throughput,
                    "target_search_time_ms": 100,  # Target <100ms
                    "avg_search_time_ms": avg_response_time * 1000,
                    "search_time_target_met": avg_response_time < 0.1,
                    "system_metrics": system_metrics
                }
            )
            
            self.results.append(metrics)
            return metrics
            
        except Exception as e:
            print(f"Database load test error: {e}")
            end_time = time.time()
            duration = end_time - start_time
            system_metrics = self.monitor.stop_monitoring()
            
            metrics = PerformanceMetrics(
                test_name="database_load",
                start_time=start_time,
                end_time=end_time,
                duration=duration,
                throughput=0,
                memory_usage_mb=system_metrics.get("memory_mb", {}).get("max", 0),
                cpu_usage_percent=system_metrics.get("cpu_percent", {}).get("avg", 0),
                success_count=0,
                error_count=1,
                response_times=[],
                avg_response_time=0,
                p95_response_time=0,
                p99_response_time=0,
                additional_metrics={"error": str(e)}
            )
            
            self.results.append(metrics)
            return metrics
    
    def run_memory_stress_test(self) -> PerformanceMetrics:
        """Test memory usage under sustained load."""
        print("💾 Running memory stress test...")
        
        # Start monitoring
        self.monitor.start_monitoring()
        start_time = time.time()
        
        response_times = []
        success_count = 0
        error_count = 0
        
        try:
            # Initialize data store manager
            manager = create_data_store_manager()
            
            # Process multiple batches to stress memory
            total_variants = 0
            for batch_num in range(10):  # 10 batches
                batch_start = time.time()
                
                # Generate large batch
                variants_data = VCFDataGenerator.generate_variant_data(self.config.num_variants // 10)
                
                try:
                    sample_data = {"id": f"MEMORY_TEST_{batch_num}", "name": f"Memory Test Batch {batch_num}"}
                    result = manager.add_sample_with_variants(sample_data, variants_data)
                    
                    if result.get("success"):
                        success_count += len(variants_data)
                        total_variants += len(variants_data)
                    else:
                        error_count += len(variants_data)
                    
                    batch_time = time.time() - batch_start
                    response_times.append(batch_time)
                    
                    # Check memory usage
                    current_memory = psutil.Process().memory_info().rss / 1024 / 1024
                    if current_memory > self.config.memory_limit_mb:
                        print(f"⚠️ Memory limit exceeded: {current_memory:.1f}MB > {self.config.memory_limit_mb}MB")
                    
                except Exception as e:
                    print(f"Memory stress batch {batch_num} error: {e}")
                    error_count += len(variants_data)
                    response_times.append(time.time() - batch_start)
            
            end_time = time.time()
            duration = end_time - start_time
            
            # Stop monitoring
            system_metrics = self.monitor.stop_monitoring()
            
            # Calculate metrics
            throughput = success_count / duration if duration > 0 else 0
            avg_response_time = statistics.mean(response_times) if response_times else 0
            p95_response_time = statistics.quantiles(response_times, n=20)[18] if len(response_times) > 20 else max(response_times) if response_times else 0
            p99_response_time = statistics.quantiles(response_times, n=100)[98] if len(response_times) > 100 else max(response_times) if response_times else 0
            
            max_memory = system_metrics.get("memory_mb", {}).get("max", 0)
            memory_limit_exceeded = max_memory > self.config.memory_limit_mb
            
            metrics = PerformanceMetrics(
                test_name="memory_stress",
                start_time=start_time,
                end_time=end_time,
                duration=duration,
                throughput=throughput,
                memory_usage_mb=max_memory,
                cpu_usage_percent=system_metrics.get("cpu_percent", {}).get("avg", 0),
                success_count=success_count,
                error_count=error_count,
                response_times=response_times,
                avg_response_time=avg_response_time,
                p95_response_time=p95_response_time,
                p99_response_time=p99_response_time,
                additional_metrics={
                    "total_variants_processed": total_variants,
                    "max_memory_mb": max_memory,
                    "memory_limit_mb": self.config.memory_limit_mb,
                    "memory_limit_exceeded": memory_limit_exceeded,
                    "system_metrics": system_metrics
                }
            )
            
            self.results.append(metrics)
            return metrics
            
        except Exception as e:
            print(f"Memory stress test error: {e}")
            end_time = time.time()
            duration = end_time - start_time
            system_metrics = self.monitor.stop_monitoring()
            
            metrics = PerformanceMetrics(
                test_name="memory_stress",
                start_time=start_time,
                end_time=end_time,
                duration=duration,
                throughput=0,
                memory_usage_mb=system_metrics.get("memory_mb", {}).get("max", 0),
                cpu_usage_percent=system_metrics.get("cpu_percent", {}).get("avg", 0),
                success_count=0,
                error_count=1,
                response_times=[],
                avg_response_time=0,
                p95_response_time=0,
                p99_response_time=0,
                additional_metrics={"error": str(e)}
            )
            
            self.results.append(metrics)
            return metrics
    
    def run_all_tests(self) -> List[PerformanceMetrics]:
        """Run all load testing scenarios."""
        print("🚀 Starting comprehensive load testing suite...")
        
        test_methods = [
            self.run_batch_processing_test,
            self.run_concurrent_analysis_test,
            self.run_database_load_test,
            self.run_memory_stress_test
        ]
        
        for test_method in test_methods:
            try:
                print(f"\n{'='*60}")
                metrics = test_method()
                self._print_test_results(metrics)
            except Exception as e:
                print(f"Test {test_method.__name__} failed: {e}")
        
        print(f"\n{'='*60}")
        print("📊 LOAD TESTING COMPLETE")
        self._print_summary()
        
        return self.results
    
    def _print_test_results(self, metrics: PerformanceMetrics):
        """Print formatted test results."""
        print(f"\n📊 {metrics.test_name.upper()} RESULTS:")
        print(f"   Duration: {metrics.duration:.2f}s")
        print(f"   Throughput: {metrics.throughput:.2f} ops/sec")
        print(f"   Success Rate: {metrics.success_count}/{metrics.success_count + metrics.error_count} ({100 * metrics.success_count / (metrics.success_count + metrics.error_count) if (metrics.success_count + metrics.error_count) > 0 else 0:.1f}%)")
        print(f"   Avg Response Time: {metrics.avg_response_time:.3f}s")
        print(f"   P95 Response Time: {metrics.p95_response_time:.3f}s")
        print(f"   Max Memory: {metrics.memory_usage_mb:.1f}MB")
        print(f"   Avg CPU: {metrics.cpu_usage_percent:.1f}%")
        
        # Test-specific metrics
        if metrics.test_name == "batch_processing":
            target_met = metrics.additional_metrics.get("throughput_ratio", 0) >= 1.0
            print(f"   Target Throughput: {'✅ MET' if target_met else '❌ NOT MET'} ({metrics.additional_metrics.get('variants_per_second', 0):.1f}/{self.config.target_throughput} variants/sec)")
        
        elif metrics.test_name == "database_load":
            target_met = metrics.additional_metrics.get("search_time_target_met", False)
            print(f"   Search Time Target: {'✅ MET' if target_met else '❌ NOT MET'} ({metrics.additional_metrics.get('avg_search_time_ms', 0):.1f}ms)")
        
        elif metrics.test_name == "memory_stress":
            memory_ok = not metrics.additional_metrics.get("memory_limit_exceeded", True)
            print(f"   Memory Limit: {'✅ OK' if memory_ok else '❌ EXCEEDED'} ({metrics.memory_usage_mb:.1f}/{self.config.memory_limit_mb}MB)")
    
    def _print_summary(self):
        """Print overall test summary."""
        if not self.results:
            print("❌ No test results available")
            return
        
        total_success = sum(m.success_count for m in self.results)
        total_errors = sum(m.error_count for m in self.results)
        overall_success_rate = total_success / (total_success + total_errors) if (total_success + total_errors) > 0 else 0
        
        print(f"\n🎯 OVERALL PERFORMANCE SUMMARY:")
        print(f"   Tests Run: {len(self.results)}")
        print(f"   Overall Success Rate: {overall_success_rate:.1%}")
        print(f"   Total Operations: {total_success + total_errors}")
        
        # Performance targets assessment
        batch_test = next((m for m in self.results if m.test_name == "batch_processing"), None)
        if batch_test:
            throughput_target_met = batch_test.additional_metrics.get("throughput_ratio", 0) >= 1.0
            print(f"   Throughput Target: {'✅ MET' if throughput_target_met else '❌ NOT MET'}")
        
        db_test = next((m for m in self.results if m.test_name == "database_load"), None)
        if db_test:
            search_target_met = db_test.additional_metrics.get("search_time_target_met", False)
            print(f"   Search Time Target: {'✅ MET' if search_target_met else '❌ NOT MET'}")
        
        memory_test = next((m for m in self.results if m.test_name == "memory_stress"), None)
        if memory_test:
            memory_ok = not memory_test.additional_metrics.get("memory_limit_exceeded", True)
            print(f"   Memory Usage: {'✅ OK' if memory_ok else '❌ EXCEEDED'}")
        
        print(f"\n💾 Results saved to: performance_reports/load_test_results.json")
    
    def save_results(self, output_path: str = "performance_reports/load_test_results.json"):
        """Save test results to JSON file."""
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        results_data = {
            "test_config": asdict(self.config),
            "test_timestamp": time.time(),
            "results": [asdict(metrics) for metrics in self.results]
        }
        
        with open(output_path, 'w') as f:
            json.dump(results_data, f, indent=2)
        
        print(f"📊 Results saved to {output_path}")


# Pytest test functions
@pytest.mark.performance
def test_batch_processing_performance():
    """Test batch processing performance."""
    config = LoadTestConfig(num_variants=1000, batch_size=100)
    runner = LoadTestRunner(config)
    
    metrics = runner.run_batch_processing_test()
    
    # Assertions
    assert metrics.success_count > 0, "No successful operations"
    assert metrics.throughput > 0, "Zero throughput"
    assert metrics.error_count == 0, f"Errors occurred: {metrics.error_count}"
    
    # Performance targets
    assert metrics.throughput >= 500, f"Throughput too low: {metrics.throughput} < 500 variants/sec"
    assert metrics.avg_response_time < 5.0, f"Response time too high: {metrics.avg_response_time}s"


@pytest.mark.performance
def test_concurrent_analysis_performance():
    """Test concurrent analysis performance."""
    config = LoadTestConfig(num_concurrent_users=5)
    runner = LoadTestRunner(config)
    
    metrics = runner.run_concurrent_analysis_test()
    
    # Assertions
    assert metrics.success_count > 0, "No successful operations"
    assert metrics.error_count == 0, f"Errors occurred: {metrics.error_count}"
    assert metrics.avg_response_time < 30.0, f"Response time too high: {metrics.avg_response_time}s"


@pytest.mark.performance
def test_database_load_performance():
    """Test database performance under load."""
    config = LoadTestConfig(num_variants=500)
    runner = LoadTestRunner(config)
    
    metrics = runner.run_database_load_test()
    
    # Assertions
    assert metrics.success_count > 0, "No successful operations"
    assert metrics.avg_response_time < 0.5, f"Search time too high: {metrics.avg_response_time}s"


@pytest.mark.performance
def test_memory_stress():
    """Test memory usage under stress."""
    config = LoadTestConfig(num_variants=2000, memory_limit_mb=1024)
    runner = LoadTestRunner(config)
    
    metrics = runner.run_memory_stress_test()
    
    # Assertions
    assert metrics.success_count > 0, "No successful operations"
    assert not metrics.additional_metrics.get("memory_limit_exceeded", True), "Memory limit exceeded"


@pytest.mark.performance
def test_comprehensive_load_testing():
    """Run comprehensive load testing suite."""
    config = LoadTestConfig(
        num_variants=1000,
        num_concurrent_users=5,
        batch_size=100,
        target_throughput=1000.0,
        memory_limit_mb=2048
    )
    
    runner = LoadTestRunner(config)
    results = runner.run_all_tests()
    
    # Save results
    runner.save_results()
    
    # Overall assertions
    assert len(results) > 0, "No test results"
    
    total_success = sum(m.success_count for m in results)
    total_errors = sum(m.error_count for m in results)
    success_rate = total_success / (total_success + total_errors) if (total_success + total_errors) > 0 else 0
    
    assert success_rate >= 0.95, f"Overall success rate too low: {success_rate:.1%}"


if __name__ == "__main__":
    # Run load tests directly
    config = LoadTestConfig(
        num_variants=1000,
        num_concurrent_users=10,
        batch_size=100,
        target_throughput=1000.0,
        memory_limit_mb=2048
    )
    
    runner = LoadTestRunner(config)
    results = runner.run_all_tests()
    runner.save_results() 