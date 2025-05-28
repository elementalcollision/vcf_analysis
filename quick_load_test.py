#!/usr/bin/env python3
"""
Quick Load Test for VCF Analysis Agent

A lightweight test to validate load testing framework and measure basic performance.
"""

import time
import statistics
import psutil
import tempfile
import os
from pathlib import Path

def generate_test_variants(num_variants: int):
    """Generate simple test variant data."""
    variants = []
    for i in range(num_variants):
        variant = {
            "id": f"test_variant_{i}",
            "chr": str((i % 22) + 1),
            "pos": 1000000 + (i * 100),
            "ref": "A",
            "alt": "G",
            "sample_id": f"SAMPLE_{i % 5}"
        }
        variants.append(variant)
    return variants

def measure_performance(func, *args, **kwargs):
    """Measure function performance."""
    process = psutil.Process()
    
    # Baseline memory
    mem_before = process.memory_info().rss / 1024 / 1024
    
    # Run function with timing
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    
    # Final memory
    mem_after = process.memory_info().rss / 1024 / 1024
    
    duration = end_time - start_time
    memory_used = mem_after - mem_before
    
    return {
        "result": result,
        "duration": duration,
        "memory_used_mb": memory_used,
        "memory_before_mb": mem_before,
        "memory_after_mb": mem_after
    }

def test_data_generation():
    """Test variant data generation performance."""
    print("🧪 Testing data generation...")
    
    def generate_data():
        return generate_test_variants(1000)
    
    metrics = measure_performance(generate_data)
    variants = metrics["result"]
    
    print(f"   ✅ Generated {len(variants)} variants")
    print(f"   ⏱️  Duration: {metrics['duration']:.3f}s")
    print(f"   💾 Memory used: {metrics['memory_used_mb']:.1f}MB")
    print(f"   🚀 Throughput: {len(variants)/metrics['duration']:.0f} variants/sec")
    
    return metrics

def test_batch_processing():
    """Test batch processing simulation."""
    print("\n🧪 Testing batch processing...")
    
    def process_batches():
        variants = generate_test_variants(500)
        batch_size = 50
        processed = 0
        
        for i in range(0, len(variants), batch_size):
            batch = variants[i:i+batch_size]
            # Simulate processing time
            time.sleep(0.001)  # 1ms per batch
            processed += len(batch)
        
        return processed
    
    metrics = measure_performance(process_batches)
    processed = metrics["result"]
    
    print(f"   ✅ Processed {processed} variants in batches")
    print(f"   ⏱️  Duration: {metrics['duration']:.3f}s")
    print(f"   💾 Memory used: {metrics['memory_used_mb']:.1f}MB")
    print(f"   🚀 Throughput: {processed/metrics['duration']:.0f} variants/sec")
    
    return metrics

def test_concurrent_simulation():
    """Test concurrent processing simulation."""
    print("\n🧪 Testing concurrent processing simulation...")
    
    import threading
    import queue
    
    def worker(work_queue, results_queue):
        """Worker thread function."""
        while True:
            try:
                item = work_queue.get(timeout=1)
                if item is None:
                    break
                # Simulate work
                time.sleep(0.001)
                results_queue.put(f"processed_{item}")
                work_queue.task_done()
            except queue.Empty:
                break
    
    def concurrent_processing():
        work_queue = queue.Queue()
        results_queue = queue.Queue()
        
        # Add work items
        for i in range(100):
            work_queue.put(f"variant_{i}")
        
        # Start workers
        workers = []
        for i in range(3):  # 3 concurrent workers
            t = threading.Thread(target=worker, args=(work_queue, results_queue))
            t.start()
            workers.append(t)
        
        # Wait for completion
        work_queue.join()
        
        # Stop workers
        for _ in workers:
            work_queue.put(None)
        for t in workers:
            t.join()
        
        # Collect results
        results = []
        while not results_queue.empty():
            results.append(results_queue.get())
        
        return len(results)
    
    metrics = measure_performance(concurrent_processing)
    processed = metrics["result"]
    
    print(f"   ✅ Processed {processed} items concurrently")
    print(f"   ⏱️  Duration: {metrics['duration']:.3f}s")
    print(f"   💾 Memory used: {metrics['memory_used_mb']:.1f}MB")
    print(f"   🚀 Throughput: {processed/metrics['duration']:.0f} items/sec")
    
    return metrics

def main():
    """Run quick load tests."""
    print("🚀 VCF Analysis Agent - Quick Load Test")
    print("=" * 50)
    
    # Run tests
    test1 = test_data_generation()
    test2 = test_batch_processing()
    test3 = test_concurrent_simulation()
    
    # Summary
    print("\n" + "=" * 50)
    print("📊 QUICK LOAD TEST SUMMARY")
    print("=" * 50)
    
    total_duration = test1["duration"] + test2["duration"] + test3["duration"]
    max_memory = max(test1["memory_after_mb"], test2["memory_after_mb"], test3["memory_after_mb"])
    
    print(f"🕒 Total test duration: {total_duration:.3f}s")
    print(f"💾 Peak memory usage: {max_memory:.1f}MB")
    print(f"🎯 Data generation: {1000/test1['duration']:.0f} variants/sec")
    print(f"🎯 Batch processing: {500/test2['duration']:.0f} variants/sec")
    print(f"🎯 Concurrent processing: {100/test3['duration']:.0f} items/sec")
    
    # Performance assessment
    data_gen_good = (1000/test1['duration']) > 500  # >500 variants/sec
    batch_good = (500/test2['duration']) > 200      # >200 variants/sec
    memory_good = max_memory < 500                   # <500MB
    
    if data_gen_good and batch_good and memory_good:
        print("\n✅ All performance targets MET! Framework is working well.")
        return 0
    else:
        print("\n⚠️ Some performance targets not optimal. Framework needs tuning.")
        return 1

if __name__ == "__main__":
    exit(main()) 