#!/usr/bin/env python3
"""
Performance Analysis Script for VCF Analysis Agent

This script performs comprehensive performance analysis including:
- Memory usage profiling
- Database query performance
- AI model response times
- Code hotspot identification
- Resource utilization analysis
"""

import sys
import os
import time
import psutil
import cProfile
import pstats
import io
import json
import tracemalloc
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcf_agent.config import SessionConfig
from vcf_agent.data_store_manager import create_data_store_manager, DataStoreConfig
from vcf_agent.lancedb_integration import VariantEmbeddingService

class PerformanceAnalyzer:
    """Comprehensive performance analyzer for VCF Agent."""
    
    def __init__(self, output_dir: str = "performance_reports"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = {}
        self.start_time = time.time()
        
    def analyze_memory_usage(self) -> Dict[str, Any]:
        """Analyze memory usage patterns."""
        print("🔍 Analyzing memory usage...")
        
        # Start memory tracing
        tracemalloc.start()
        
        # Get initial memory state
        process = psutil.Process()
        initial_memory = process.memory_info()
        
        # Test data store manager initialization
        config = SessionConfig(model_provider="ollama")
        data_config = DataStoreConfig(performance_monitoring=True)
        
        # Memory snapshot before initialization
        snapshot1 = tracemalloc.take_snapshot()
        
        # Initialize manager
        manager = create_data_store_manager(
            lancedb_path="./test_lancedb",
            kuzu_path="./test_kuzu_db",
            session_config=config,
            data_store_config=data_config
        )
        
        # Memory snapshot after initialization
        snapshot2 = tracemalloc.take_snapshot()
        
        # Calculate memory differences
        top_stats = snapshot2.compare_to(snapshot1, 'lineno')
        
        final_memory = process.memory_info()
        
        memory_analysis = {
            "initial_memory_mb": initial_memory.rss / 1024 / 1024,
            "final_memory_mb": final_memory.rss / 1024 / 1024,
            "memory_increase_mb": (final_memory.rss - initial_memory.rss) / 1024 / 1024,
            "top_memory_allocations": [
                {
                    "file": stat.traceback.format()[0],
                    "size_mb": stat.size / 1024 / 1024,
                    "count": stat.count
                }
                for stat in top_stats[:10]
            ]
        }
        
        tracemalloc.stop()
        return memory_analysis
    
    def analyze_database_performance(self) -> Dict[str, Any]:
        """Analyze database operation performance."""
        print("🔍 Analyzing database performance...")
        
        config = SessionConfig(model_provider="ollama")
        data_config = DataStoreConfig(
            performance_monitoring=True,
            max_workers=4,
            sync_batch_size=100
        )
        
        manager = create_data_store_manager(
            lancedb_path="./test_lancedb",
            kuzu_path="./test_kuzu_db",
            session_config=config,
            data_store_config=data_config
        )
        
        # Test data
        sample_data = {
            "id": "PERF_TEST_001",
            "name": "Performance Test Sample",
            "type": "germline"
        }
        
        variants_data = [
            {
                "id": f"chr1-{100000 + i}-A-G",
                "chr": "1",
                "pos": 100000 + i,
                "ref": "A",
                "alt": "G",
                "gene_symbol": f"GENE_{i % 10}",
                "clinical_significance": "Benign" if i % 2 == 0 else "Pathogenic",
                "genotype": "0/1",
                "quality_score": 99.0
            }
            for i in range(100)  # Test with 100 variants
        ]
        
        # Measure ingestion performance
        start_time = time.time()
        result = manager.add_sample_with_variants(
            sample_data=sample_data,
            variants_data=variants_data
        )
        ingestion_time = time.time() - start_time
        
        # Measure search performance
        start_time = time.time()
        search_results = manager.search_variants(
            query="pathogenic variants",
            search_type="hybrid",
            limit=10
        )
        search_time = time.time() - start_time
        
        # Measure analysis performance
        start_time = time.time()
        analysis = manager.get_sample_analysis("PERF_TEST_001")
        analysis_time = time.time() - start_time
        
        # Get performance statistics
        stats = manager._get_performance_statistics()
        
        return {
            "ingestion": {
                "time_seconds": ingestion_time,
                "variants_per_second": len(variants_data) / ingestion_time,
                "total_variants": len(variants_data)
            },
            "search": {
                "time_seconds": search_time,
                "results_count": len(search_results.get("results", [])),
                "queries_per_second": 1 / search_time if search_time > 0 else 0
            },
            "analysis": {
                "time_seconds": analysis_time,
                "variants_analyzed": len(analysis.get("variants", []))
            },
            "performance_stats": stats
        }
    
    def analyze_ai_performance(self) -> Dict[str, Any]:
        """Analyze AI model performance."""
        print("🔍 Analyzing AI model performance...")
        
        # Test embedding generation performance
        embedding_service = VariantEmbeddingService(
            provider="ollama",
            model="qwen3:4b"
        )
        
        test_descriptions = [
            "Pathogenic variant in BRCA1 gene",
            "Benign variant in TP53 gene",
            "Variant of uncertain significance in EGFR",
            "Likely pathogenic variant in CFTR gene",
            "Benign variant in APOE gene"
        ]
        
        embedding_times = []
        for desc in test_descriptions:
            start_time = time.time()
            try:
                embedding = embedding_service.generate_embedding_sync(desc)
                embedding_time = time.time() - start_time
                embedding_times.append(embedding_time)
            except Exception as e:
                print(f"Warning: Embedding generation failed: {e}")
                embedding_times.append(None)
        
        valid_times = [t for t in embedding_times if t is not None]
        
        return {
            "embedding_generation": {
                "total_tests": len(test_descriptions),
                "successful_tests": len(valid_times),
                "avg_time_seconds": sum(valid_times) / len(valid_times) if valid_times else 0,
                "min_time_seconds": min(valid_times) if valid_times else 0,
                "max_time_seconds": max(valid_times) if valid_times else 0,
                "embeddings_per_second": len(valid_times) / sum(valid_times) if valid_times and sum(valid_times) > 0 else 0
            }
        }
    
    def analyze_code_hotspots(self) -> Dict[str, Any]:
        """Analyze code performance hotspots using cProfile."""
        print("🔍 Analyzing code hotspots...")
        
        # Create profiler
        profiler = cProfile.Profile()
        
        # Profile a complete workflow
        profiler.enable()
        
        try:
            config = SessionConfig(model_provider="ollama")
            manager = create_data_store_manager(
                lancedb_path="./test_lancedb",
                kuzu_path="./test_kuzu_db",
                session_config=config
            )
            
            # Simulate typical workflow
            sample_data = {"id": "HOTSPOT_TEST", "name": "Test", "type": "germline"}
            variants_data = [
                {
                    "id": f"chr1-{200000 + i}-A-G",
                    "chr": "1",
                    "pos": 200000 + i,
                    "ref": "A",
                    "alt": "G",
                    "gene_symbol": "TEST_GENE",
                    "clinical_significance": "Benign"
                }
                for i in range(50)
            ]
            
            manager.add_sample_with_variants(sample_data, variants_data)
            manager.search_variants("test variants", limit=5)
            manager.get_sample_analysis("HOTSPOT_TEST")
            
        except Exception as e:
            print(f"Warning: Profiling workflow failed: {e}")
        
        profiler.disable()
        
        # Analyze profiling results
        s = io.StringIO()
        ps = pstats.Stats(profiler, stream=s)
        ps.sort_stats('cumulative')
        ps.print_stats(20)  # Top 20 functions
        
        profile_output = s.getvalue()
        
        # Extract top functions
        lines = profile_output.split('\n')
        hotspots = []
        for line in lines[5:25]:  # Skip header lines
            if line.strip() and not line.startswith('Ordered by'):
                parts = line.split()
                if len(parts) >= 6:
                    hotspots.append({
                        "ncalls": parts[0],
                        "tottime": parts[1],
                        "cumtime": parts[3],
                        "function": ' '.join(parts[5:])
                    })
        
        return {
            "profile_output": profile_output,
            "top_hotspots": hotspots[:10]
        }
    
    def analyze_resource_utilization(self) -> Dict[str, Any]:
        """Analyze system resource utilization."""
        print("🔍 Analyzing resource utilization...")
        
        process = psutil.Process()
        
        # CPU usage
        cpu_percent = process.cpu_percent(interval=1)
        
        # Memory usage
        memory_info = process.memory_info()
        memory_percent = process.memory_percent()
        
        # Disk I/O
        io_counters = process.io_counters()
        
        # System-wide stats
        system_cpu = psutil.cpu_percent(interval=1)
        system_memory = psutil.virtual_memory()
        system_disk = psutil.disk_usage('/')
        
        return {
            "process": {
                "cpu_percent": cpu_percent,
                "memory_mb": memory_info.rss / 1024 / 1024,
                "memory_percent": memory_percent,
                "read_bytes": io_counters.read_bytes,
                "write_bytes": io_counters.write_bytes
            },
            "system": {
                "cpu_percent": system_cpu,
                "memory_total_gb": system_memory.total / 1024 / 1024 / 1024,
                "memory_available_gb": system_memory.available / 1024 / 1024 / 1024,
                "memory_percent": system_memory.percent,
                "disk_total_gb": system_disk.total / 1024 / 1024 / 1024,
                "disk_free_gb": system_disk.free / 1024 / 1024 / 1024,
                "disk_percent": (system_disk.used / system_disk.total) * 100
            }
        }
    
    def generate_optimization_recommendations(self) -> List[str]:
        """Generate optimization recommendations based on analysis."""
        recommendations = []
        
        # Memory optimization
        if self.results.get("memory", {}).get("memory_increase_mb", 0) > 100:
            recommendations.append(
                "Consider implementing memory pooling for database connections to reduce initialization overhead"
            )
        
        # Database performance
        db_perf = self.results.get("database", {})
        if db_perf.get("ingestion", {}).get("variants_per_second", 0) < 1000:
            recommendations.append(
                "Optimize batch ingestion by increasing batch size or implementing parallel processing"
            )
        
        if db_perf.get("search", {}).get("time_seconds", 0) > 1.0:
            recommendations.append(
                "Consider adding database indexes or optimizing query patterns for faster search"
            )
        
        # AI performance
        ai_perf = self.results.get("ai", {}).get("embedding_generation", {})
        if ai_perf.get("avg_time_seconds", 0) > 2.0:
            recommendations.append(
                "Implement embedding caching to reduce AI model call overhead"
            )
        
        # Resource utilization
        resource = self.results.get("resources", {})
        if resource.get("system", {}).get("memory_percent", 0) > 80:
            recommendations.append(
                "System memory usage is high - consider implementing memory limits or garbage collection optimization"
            )
        
        if resource.get("system", {}).get("cpu_percent", 0) > 80:
            recommendations.append(
                "High CPU usage detected - consider implementing async processing or reducing computational complexity"
            )
        
        return recommendations
    
    def run_full_analysis(self) -> Dict[str, Any]:
        """Run complete performance analysis."""
        print("🚀 Starting comprehensive performance analysis...")
        
        try:
            self.results["memory"] = self.analyze_memory_usage()
        except Exception as e:
            print(f"Memory analysis failed: {e}")
            self.results["memory"] = {"error": str(e)}
        
        try:
            self.results["database"] = self.analyze_database_performance()
        except Exception as e:
            print(f"Database analysis failed: {e}")
            self.results["database"] = {"error": str(e)}
        
        try:
            self.results["ai"] = self.analyze_ai_performance()
        except Exception as e:
            print(f"AI analysis failed: {e}")
            self.results["ai"] = {"error": str(e)}
        
        try:
            self.results["hotspots"] = self.analyze_code_hotspots()
        except Exception as e:
            print(f"Hotspot analysis failed: {e}")
            self.results["hotspots"] = {"error": str(e)}
        
        try:
            self.results["resources"] = self.analyze_resource_utilization()
        except Exception as e:
            print(f"Resource analysis failed: {e}")
            self.results["resources"] = {"error": str(e)}
        
        # Generate recommendations
        self.results["recommendations"] = self.generate_optimization_recommendations()
        
        # Add metadata
        self.results["metadata"] = {
            "analysis_date": datetime.now().isoformat(),
            "total_duration_seconds": time.time() - self.start_time,
            "python_version": sys.version,
            "platform": os.uname().sysname if hasattr(os, 'uname') else 'unknown'
        }
        
        return self.results
    
    def save_results(self, filename: Optional[str] = None) -> str:
        """Save analysis results to file."""
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"performance_analysis_{timestamp}.json"
        
        filepath = self.output_dir / filename
        
        with open(filepath, 'w') as f:
            json.dump(self.results, f, indent=2, default=str)
        
        return str(filepath)
    
    def print_summary(self):
        """Print analysis summary."""
        print("\n" + "="*60)
        print("📊 PERFORMANCE ANALYSIS SUMMARY")
        print("="*60)
        
        # Memory summary
        memory = self.results.get("memory", {})
        if "error" not in memory:
            print(f"💾 Memory Usage: {memory.get('memory_increase_mb', 0):.1f} MB increase")
        
        # Database summary
        database = self.results.get("database", {})
        if "error" not in database:
            ingestion = database.get("ingestion", {})
            search = database.get("search", {})
            print(f"🗄️  Database Performance:")
            print(f"   - Ingestion: {ingestion.get('variants_per_second', 0):.0f} variants/sec")
            print(f"   - Search: {search.get('time_seconds', 0):.3f} seconds")
        
        # AI summary
        ai = self.results.get("ai", {})
        if "error" not in ai:
            embedding = ai.get("embedding_generation", {})
            print(f"🤖 AI Performance:")
            print(f"   - Embeddings: {embedding.get('avg_time_seconds', 0):.3f} sec avg")
        
        # Resource summary
        resources = self.results.get("resources", {})
        if "error" not in resources:
            system = resources.get("system", {})
            print(f"⚡ Resource Usage:")
            print(f"   - CPU: {system.get('cpu_percent', 0):.1f}%")
            print(f"   - Memory: {system.get('memory_percent', 0):.1f}%")
        
        # Recommendations
        recommendations = self.results.get("recommendations", [])
        if recommendations:
            print(f"\n💡 Optimization Recommendations:")
            for i, rec in enumerate(recommendations, 1):
                print(f"   {i}. {rec}")
        
        print("\n" + "="*60)

def main():
    """Main function to run performance analysis."""
    analyzer = PerformanceAnalyzer()
    
    try:
        results = analyzer.run_full_analysis()
        filepath = analyzer.save_results()
        
        analyzer.print_summary()
        print(f"\n📄 Full results saved to: {filepath}")
        
        # Cleanup test databases
        import shutil
        for path in ["./test_lancedb", "./test_kuzu_db"]:
            if os.path.exists(path):
                shutil.rmtree(path)
        
        return 0
        
    except Exception as e:
        print(f"❌ Performance analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main()) 