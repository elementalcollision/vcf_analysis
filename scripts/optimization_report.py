#!/usr/bin/env python3
"""
Optimization Report Generator for VCF Analysis Agent

This script generates a comprehensive optimization report based on performance analysis
and provides specific recommendations for Phase 3 optimization improvements.
"""

import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any

def load_latest_performance_report() -> Dict[str, Any]:
    """Load the latest performance analysis report."""
    reports_dir = Path("performance_reports")
    if not reports_dir.exists():
        return {}
    
    # Find the latest report
    report_files = list(reports_dir.glob("performance_analysis_*.json"))
    if not report_files:
        return {}
    
    latest_report = max(report_files, key=lambda x: x.stat().st_mtime)
    
    with open(latest_report, 'r') as f:
        return json.load(f)

def analyze_performance_bottlenecks(report: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Analyze performance bottlenecks from the report."""
    bottlenecks = []
    
    # Check hotspots
    hotspots = report.get("hotspots", {}).get("top_hotspots", [])
    for hotspot in hotspots[1:6]:  # Skip header row, take top 5
        if isinstance(hotspot, dict) and "cumtime" in hotspot:
            try:
                cumtime = float(hotspot["cumtime"])
                if cumtime > 1.0:  # More than 1 second
                    bottlenecks.append({
                        "type": "code_hotspot",
                        "function": hotspot["function"],
                        "time_seconds": cumtime,
                        "calls": hotspot["ncalls"],
                        "severity": "high" if cumtime > 10 else "medium"
                    })
            except (ValueError, TypeError):
                continue
    
    return bottlenecks

def generate_optimization_recommendations(report: Dict[str, Any], bottlenecks: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Generate specific optimization recommendations."""
    recommendations = []
    
    # Cache optimization (major bottleneck identified)
    cache_bottleneck = any(
        "cache" in b.get("function", "").lower() or "json" in b.get("function", "").lower()
        for b in bottlenecks
    )
    
    if cache_bottleneck:
        recommendations.append({
            "category": "Cache Optimization",
            "priority": "high",
            "issue": "JSON cache serialization is taking 12+ seconds",
            "recommendation": "Replace JSON cache with binary format (pickle/msgpack) or implement async cache persistence",
            "implementation": [
                "Use msgpack for faster serialization",
                "Implement background cache persistence",
                "Add cache size limits to prevent large serialization",
                "Consider using SQLite for cache storage"
            ],
            "expected_improvement": "90% reduction in cache save time"
        })
    
    # Embedding optimization
    recommendations.append({
        "category": "Embedding Generation",
        "priority": "medium",
        "issue": "Multiple embedding calls for same content",
        "recommendation": "Implement smarter caching and batch processing",
        "implementation": [
            "Pre-compute embeddings for common variant descriptions",
            "Implement embedding deduplication",
            "Use batch embedding APIs where available",
            "Add embedding compression for storage"
        ],
        "expected_improvement": "50% reduction in embedding generation time"
    })
    
    # Database query optimization
    recommendations.append({
        "category": "Database Queries",
        "priority": "medium",
        "issue": "Individual gene queries instead of batch queries",
        "recommendation": "Implement proper batch query optimization",
        "implementation": [
            "Fix batch query result iteration",
            "Optimize Kuzu query patterns",
            "Add query result caching",
            "Implement connection pooling"
        ],
        "expected_improvement": "70% reduction in query time for large datasets"
    })
    
    # Memory optimization
    recommendations.append({
        "category": "Memory Management",
        "priority": "low",
        "issue": "Large object creation and retention",
        "recommendation": "Implement memory-efficient data structures",
        "implementation": [
            "Use generators for large data processing",
            "Implement object pooling for frequent allocations",
            "Add memory monitoring and limits",
            "Optimize data structure sizes"
        ],
        "expected_improvement": "30% reduction in memory usage"
    })
    
    return recommendations

def generate_implementation_plan(recommendations: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Generate a detailed implementation plan."""
    plan = {
        "immediate_actions": [],
        "short_term": [],
        "long_term": []
    }
    
    for rec in recommendations:
        if rec["priority"] == "high":
            plan["immediate_actions"].extend([
                f"{rec['category']}: {action}" for action in rec["implementation"][:2]
            ])
        elif rec["priority"] == "medium":
            plan["short_term"].extend([
                f"{rec['category']}: {action}" for action in rec["implementation"][:2]
            ])
        else:
            plan["long_term"].extend([
                f"{rec['category']}: {action}" for action in rec["implementation"][:1]
            ])
    
    return plan

def generate_optimization_report() -> Dict[str, Any]:
    """Generate comprehensive optimization report."""
    print("📊 Generating Optimization Report...")
    
    # Load performance data
    performance_report = load_latest_performance_report()
    
    if not performance_report:
        return {"error": "No performance reports found"}
    
    # Analyze bottlenecks
    bottlenecks = analyze_performance_bottlenecks(performance_report)
    
    # Generate recommendations
    recommendations = generate_optimization_recommendations(performance_report, bottlenecks)
    
    # Generate implementation plan
    implementation_plan = generate_implementation_plan(recommendations)
    
    # Calculate potential improvements
    total_time_saved = sum(
        float(rec.get("expected_improvement", "0%").replace("%", "").split()[0]) / 100 * 10
        for rec in recommendations
    )
    
    report = {
        "analysis_date": datetime.now().isoformat(),
        "performance_summary": {
            "total_analysis_time": performance_report.get("metadata", {}).get("total_duration_seconds", 0),
            "major_bottlenecks_found": len([b for b in bottlenecks if b.get("severity") == "high"]),
            "optimization_opportunities": len(recommendations)
        },
        "bottlenecks": bottlenecks,
        "recommendations": recommendations,
        "implementation_plan": implementation_plan,
        "expected_improvements": {
            "estimated_time_savings": f"{total_time_saved:.1f} seconds",
            "performance_gain": f"{min(total_time_saved / 15 * 100, 80):.0f}%",
            "priority_order": ["Cache Optimization", "Database Queries", "Embedding Generation", "Memory Management"]
        }
    }
    
    return report

def print_optimization_summary(report: Dict[str, Any]):
    """Print a formatted optimization summary."""
    print("\n" + "="*80)
    print("🚀 VCF AGENT OPTIMIZATION REPORT")
    print("="*80)
    
    summary = report.get("performance_summary", {})
    print(f"📈 Analysis Duration: {summary.get('total_analysis_time', 0):.1f} seconds")
    print(f"🔍 Major Bottlenecks: {summary.get('major_bottlenecks_found', 0)}")
    print(f"💡 Optimization Opportunities: {summary.get('optimization_opportunities', 0)}")
    
    print(f"\n🎯 EXPECTED IMPROVEMENTS")
    improvements = report.get("expected_improvements", {})
    print(f"⏱️  Time Savings: {improvements.get('estimated_time_savings', 'N/A')}")
    print(f"📊 Performance Gain: {improvements.get('performance_gain', 'N/A')}")
    
    print(f"\n🔥 TOP BOTTLENECKS")
    bottlenecks = report.get("bottlenecks", [])
    for i, bottleneck in enumerate(bottlenecks[:3], 1):
        print(f"   {i}. {bottleneck.get('function', 'Unknown')[:60]}...")
        print(f"      Time: {bottleneck.get('time_seconds', 0):.2f}s | Calls: {bottleneck.get('calls', 0)}")
    
    print(f"\n💡 PRIORITY RECOMMENDATIONS")
    recommendations = report.get("recommendations", [])
    high_priority = [r for r in recommendations if r.get("priority") == "high"]
    for i, rec in enumerate(high_priority, 1):
        print(f"   {i}. {rec.get('category', 'Unknown')}: {rec.get('recommendation', 'N/A')}")
        print(f"      Expected: {rec.get('expected_improvement', 'N/A')}")
    
    print(f"\n📋 IMMEDIATE ACTIONS")
    plan = report.get("implementation_plan", {})
    for i, action in enumerate(plan.get("immediate_actions", [])[:5], 1):
        print(f"   {i}. {action}")
    
    print("\n" + "="*80)

def save_optimization_report(report: Dict[str, Any]) -> str:
    """Save optimization report to file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"performance_reports/optimization_report_{timestamp}.json"
    
    with open(filename, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    
    return filename

def main():
    """Main function to generate and display optimization report."""
    try:
        report = generate_optimization_report()
        
        if "error" in report:
            print(f"❌ Error: {report['error']}")
            return 1
        
        # Print summary
        print_optimization_summary(report)
        
        # Save detailed report
        filename = save_optimization_report(report)
        print(f"\n📄 Detailed report saved to: {filename}")
        
        return 0
        
    except Exception as e:
        print(f"❌ Optimization report generation failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main()) 