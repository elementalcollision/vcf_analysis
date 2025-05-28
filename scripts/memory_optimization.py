#!/usr/bin/env python3
"""
Memory Optimization Script for VCF Analysis Agent

This script identifies and implements memory optimizations:
- Memory leak detection
- Object size analysis
- Garbage collection optimization
- Memory pool implementation
"""

import sys
import gc
import tracemalloc
import psutil
import time
import weakref
from pathlib import Path
from typing import Dict, List, Any, Optional
from collections import defaultdict
import logging

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcf_agent.config import SessionConfig
from vcf_agent.data_store_manager import create_data_store_manager

logger = logging.getLogger(__name__)

class MemoryProfiler:
    """Memory profiling and optimization utilities."""
    
    def __init__(self):
        self.snapshots = []
        self.object_tracker = weakref.WeakSet()
        self.start_memory = None
        
    def start_profiling(self):
        """Start memory profiling."""
        tracemalloc.start()
        self.start_memory = psutil.Process().memory_info().rss
        gc.collect()  # Clean start
        
    def take_snapshot(self, label: str = ""):
        """Take a memory snapshot."""
        snapshot = tracemalloc.take_snapshot()
        current_memory = psutil.Process().memory_info().rss
        
        self.snapshots.append({
            'label': label,
            'snapshot': snapshot,
            'memory_rss': current_memory,
            'timestamp': time.time()
        })
        
        return snapshot
    
    def analyze_memory_growth(self) -> Dict[str, Any]:
        """Analyze memory growth between snapshots."""
        if len(self.snapshots) < 2:
            return {"error": "Need at least 2 snapshots"}
        
        latest = self.snapshots[-1]
        previous = self.snapshots[-2]
        
        # Compare snapshots
        top_stats = latest['snapshot'].compare_to(previous['snapshot'], 'lineno')
        
        # Memory growth
        memory_growth = latest['memory_rss'] - previous['memory_rss']
        
        analysis = {
            'memory_growth_bytes': memory_growth,
            'memory_growth_mb': memory_growth / 1024 / 1024,
            'time_elapsed': latest['timestamp'] - previous['timestamp'],
            'top_allocations': []
        }
        
        # Top memory allocations
        for stat in top_stats[:10]:
            analysis['top_allocations'].append({
                'file': stat.traceback.format()[0] if stat.traceback else 'unknown',
                'size_diff_mb': stat.size_diff / 1024 / 1024,
                'count_diff': stat.count_diff,
                'size_mb': stat.size / 1024 / 1024
            })
        
        return analysis
    
    def detect_memory_leaks(self) -> List[Dict[str, Any]]:
        """Detect potential memory leaks."""
        if len(self.snapshots) < 3:
            return []
        
        leaks = []
        
        # Compare first and last snapshots
        first = self.snapshots[0]
        last = self.snapshots[-1]
        
        top_stats = last['snapshot'].compare_to(first['snapshot'], 'lineno')
        
        # Look for consistently growing allocations
        for stat in top_stats:
            if stat.size_diff > 1024 * 1024:  # > 1MB growth
                leaks.append({
                    'file': stat.traceback.format()[0] if stat.traceback else 'unknown',
                    'size_growth_mb': stat.size_diff / 1024 / 1024,
                    'count_growth': stat.count_diff,
                    'current_size_mb': stat.size / 1024 / 1024,
                    'severity': 'high' if stat.size_diff > 10 * 1024 * 1024 else 'medium'
                })
        
        return leaks
    
    def optimize_garbage_collection(self):
        """Optimize garbage collection settings."""
        # Get current GC stats
        gc_stats = gc.get_stats()
        
        # Force full garbage collection
        collected = gc.collect()
        
        # Adjust GC thresholds for better performance
        # Default is (700, 10, 10), we'll make it more aggressive
        gc.set_threshold(500, 8, 8)
        
        return {
            'objects_collected': collected,
            'gc_stats_before': gc_stats,
            'new_thresholds': gc.get_threshold()
        }

class ObjectSizeAnalyzer:
    """Analyze object sizes and memory usage patterns."""
    
    @staticmethod
    def get_object_size(obj) -> int:
        """Get deep size of an object."""
        size = sys.getsizeof(obj)
        
        if isinstance(obj, dict):
            size += sum(ObjectSizeAnalyzer.get_object_size(k) + ObjectSizeAnalyzer.get_object_size(v) 
                       for k, v in obj.items())
        elif isinstance(obj, (list, tuple, set)):
            size += sum(ObjectSizeAnalyzer.get_object_size(item) for item in obj)
        
        return size
    
    @staticmethod
    def analyze_large_objects() -> List[Dict[str, Any]]:
        """Find large objects in memory."""
        large_objects = []
        
        # Get all objects
        for obj in gc.get_objects():
            try:
                size = sys.getsizeof(obj)
                if size > 1024 * 1024:  # > 1MB
                    large_objects.append({
                        'type': type(obj).__name__,
                        'size_mb': size / 1024 / 1024,
                        'id': id(obj),
                        'repr': str(obj)[:100] + '...' if len(str(obj)) > 100 else str(obj)
                    })
            except (TypeError, ReferenceError):
                continue
        
        return sorted(large_objects, key=lambda x: x['size_mb'], reverse=True)

class MemoryOptimizer:
    """Main memory optimization coordinator."""
    
    def __init__(self):
        self.profiler = MemoryProfiler()
        self.analyzer = ObjectSizeAnalyzer()
        
    def run_memory_optimization_test(self) -> Dict[str, Any]:
        """Run comprehensive memory optimization test."""
        print("🧠 Starting memory optimization analysis...")
        
        results = {
            'initial_memory': None,
            'snapshots': [],
            'optimizations': [],
            'recommendations': []
        }
        
        # Start profiling
        self.profiler.start_profiling()
        initial_memory = psutil.Process().memory_info().rss / 1024 / 1024
        results['initial_memory'] = initial_memory
        
        print(f"Initial memory usage: {initial_memory:.1f} MB")
        
        # Take baseline snapshot
        self.profiler.take_snapshot("baseline")
        
        # Test data store manager initialization
        print("Testing data store manager initialization...")
        config = SessionConfig(model_provider="ollama")
        manager = create_data_store_manager(session_config=config)
        
        # Snapshot after initialization
        self.profiler.take_snapshot("after_init")
        
        # Test data ingestion
        print("Testing data ingestion...")
        sample_data = {
            "id": "MEMORY_TEST_001",
            "name": "Memory Test Sample",
            "type": "germline"
        }
        
        variants_data = [
            {
                "id": f"chr1-{300000 + i}-A-G",
                "chr": "1",
                "pos": 300000 + i,
                "ref": "A",
                "alt": "G",
                "gene_symbol": f"GENE_{i % 5}",
                "clinical_significance": "Benign",
                "genotype": "0/1",
                "quality_score": 99.0
            }
            for i in range(200)  # Test with 200 variants
        ]
        
        manager.add_sample_with_variants(sample_data, variants_data)
        
        # Snapshot after ingestion
        self.profiler.take_snapshot("after_ingestion")
        
        # Test multiple searches
        print("Testing multiple searches...")
        for i in range(10):
            manager.search_variants(f"test query {i}", limit=10)
        
        # Snapshot after searches
        self.profiler.take_snapshot("after_searches")
        
        # Analyze memory growth
        growth_analysis = self.profiler.analyze_memory_growth()
        results['snapshots'].append(growth_analysis)
        
        # Detect memory leaks
        leaks = self.profiler.detect_memory_leaks()
        if leaks:
            results['memory_leaks'] = leaks
        
        # Optimize garbage collection
        gc_optimization = self.profiler.optimize_garbage_collection()
        results['optimizations'].append({
            'type': 'garbage_collection',
            'details': gc_optimization
        })
        
        # Analyze large objects
        large_objects = self.analyzer.analyze_large_objects()
        if large_objects:
            results['large_objects'] = large_objects[:10]  # Top 10
        
        # Final memory check
        final_memory = psutil.Process().memory_info().rss / 1024 / 1024
        results['final_memory'] = final_memory
        results['total_memory_increase'] = final_memory - initial_memory
        
        # Generate recommendations
        recommendations = self._generate_memory_recommendations(results)
        results['recommendations'] = recommendations
        
        return results
    
    def _generate_memory_recommendations(self, results: Dict[str, Any]) -> List[str]:
        """Generate memory optimization recommendations."""
        recommendations = []
        
        # Memory growth recommendations
        if results.get('total_memory_increase', 0) > 100:  # > 100MB
            recommendations.append(
                "High memory usage detected. Consider implementing object pooling and memory limits."
            )
        
        # Memory leak recommendations
        if results.get('memory_leaks'):
            recommendations.append(
                "Potential memory leaks detected. Review object lifecycle management and implement weak references where appropriate."
            )
        
        # Large object recommendations
        if results.get('large_objects'):
            recommendations.append(
                "Large objects found in memory. Consider implementing lazy loading and data streaming for large datasets."
            )
        
        # Garbage collection recommendations
        gc_stats = results.get('optimizations', [{}])[0].get('details', {})
        if gc_stats.get('objects_collected', 0) > 1000:
            recommendations.append(
                "High garbage collection activity. Consider optimizing object creation patterns and implementing object reuse."
            )
        
        return recommendations
    
    def implement_memory_optimizations(self) -> Dict[str, Any]:
        """Implement specific memory optimizations."""
        optimizations = {}
        
        # 1. Optimize garbage collection
        gc.set_threshold(500, 8, 8)  # More aggressive GC
        optimizations['gc_optimization'] = "Implemented aggressive garbage collection thresholds"
        
        # 2. Enable debug mode for memory tracking
        if hasattr(gc, 'set_debug'):
            gc.set_debug(gc.DEBUG_STATS)
            optimizations['gc_debug'] = "Enabled garbage collection debugging"
        
        # 3. Force garbage collection
        collected = gc.collect()
        optimizations['gc_collection'] = f"Collected {collected} objects"
        
        return optimizations

def main():
    """Main function to run memory optimization."""
    optimizer = MemoryOptimizer()
    
    try:
        # Run optimization test
        results = optimizer.run_memory_optimization_test()
        
        # Print summary
        print("\n" + "="*60)
        print("🧠 MEMORY OPTIMIZATION SUMMARY")
        print("="*60)
        
        print(f"Initial Memory: {results['initial_memory']:.1f} MB")
        print(f"Final Memory: {results['final_memory']:.1f} MB")
        print(f"Memory Increase: {results['total_memory_increase']:.1f} MB")
        
        if results.get('memory_leaks'):
            print(f"\n⚠️  Memory Leaks Detected: {len(results['memory_leaks'])}")
            for leak in results['memory_leaks'][:3]:
                print(f"   - {leak['file']}: +{leak['size_growth_mb']:.1f} MB")
        
        if results.get('large_objects'):
            print(f"\n📦 Large Objects: {len(results['large_objects'])}")
            for obj in results['large_objects'][:3]:
                print(f"   - {obj['type']}: {obj['size_mb']:.1f} MB")
        
        if results.get('recommendations'):
            print(f"\n💡 Recommendations:")
            for i, rec in enumerate(results['recommendations'], 1):
                print(f"   {i}. {rec}")
        
        # Implement optimizations
        print(f"\n🔧 Implementing optimizations...")
        optimizations = optimizer.implement_memory_optimizations()
        for opt_type, description in optimizations.items():
            print(f"   ✅ {description}")
        
        print("\n" + "="*60)
        
        return 0
        
    except Exception as e:
        print(f"❌ Memory optimization failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main()) 