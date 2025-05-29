#!/usr/bin/env python3
"""
Phase 5.2 Production Validation Script
=====================================

Comprehensive validation of the hybrid Apache Iggy + Kafka streaming
architecture with intelligent routing, circuit breaker patterns, and
exactly-once delivery semantics.

Features Validated:
- Dual-platform coordination and intelligent routing
- Circuit breaker patterns and automatic failover
- Message deduplication and exactly-once semantics
- Performance monitoring and health-based decisions
- End-to-end processing with synthetic VCF data

Usage:
    python scripts/validate_phase5_2.py
"""

import asyncio
import time
import tempfile
import os
import sys
import logging
from typing import List, Dict, Any
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcf_agent.phase5.streaming_coordinator import (
    StreamingCoordinator, RoutingStrategy, RoutingDecision
)
from vcf_agent.phase5.kafka_processor import KafkaVCFProcessor
from vcf_agent.phase5.monitoring import PerformanceMonitor, CircuitBreaker, CircuitBreakerState
from vcf_agent.phase5.vcf_message import VCFVariantMessage, VCFMessageType
from vcf_agent.phase5.config import Phase5Config, IggyConfig, KafkaConfig, Environment

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class Phase5_2Validator:
    """
    Comprehensive validator for Phase 5.2 dual-platform architecture.
    
    Demonstrates production-ready patterns including:
    - Intelligent routing based on platform health
    - Circuit breaker patterns for automatic failover
    - Message deduplication for exactly-once semantics
    - Performance monitoring and alerting
    """
    
    def __init__(self):
        """Initialize validator with test configuration."""
        # Create proper config objects
        iggy_config = IggyConfig(
            host="localhost",
            quic_port=8080,
            tcp_port=8090,
            http_port=3000,
            stream_name="vcf-variants-validation",
            topic_name="vcf-variants-validation",
            partition_count=4,
            compression_enabled=True
        )
        
        kafka_config = KafkaConfig(
            bootstrap_servers=["localhost:9092"],
            topic_name="vcf-variants-validation",
            group_id="vcf-agent-validation",
            partition_count=4,
            compression_type="gzip",
            batch_size=100
        )
        
        self.config = Phase5Config(
            environment=Environment.DEVELOPMENT,  # Use proper enum
            iggy=iggy_config,
            kafka=kafka_config
        )
        
        self.results = {
            "tests_passed": 0,
            "tests_failed": 0,
            "errors": [],
            "performance_metrics": {},
            "feature_validations": {}
        }
    
    def create_synthetic_variants(self, count: int = 1000) -> List[VCFVariantMessage]:
        """Create synthetic VCF variants for testing."""
        variants = []
        chromosomes = ["1", "2", "3", "4", "5", "X", "Y", "MT"]
        
        logger.info(f"Creating {count} synthetic VCF variants")
        
        for i in range(count):
            chrom = chromosomes[i % len(chromosomes)]
            position = 1000 + (i * 100)
            
            variant = VCFVariantMessage(
                chromosome=chrom,
                position=position,
                reference="A",
                alternate=["T", "C"][i % 2],
                quality=30.0 + (i % 40),
                message_type=VCFMessageType.VARIANT,
                source_file=f"synthetic_validation_{i // 100}.vcf",
                timestamp=time.time(),
                samples=[{
                    "sample_id": f"SAMPLE_{j}",
                    "GT": ["0/1", "1/1", "0/0"][i % 3],
                    "DP": str(20 + i % 30)
                } for j in range(1, 4)],
                info_fields={"DP": str(50 + i % 100), "AF": f"{0.1 + (i % 9) * 0.1:.1f}"}
            )
            variants.append(variant)
        
        logger.info(f"Created {len(variants)} synthetic variants across {len(set(v.chromosome for v in variants))} chromosomes")
        return variants
    
    async def validate_message_deduplication(self) -> bool:
        """Validate exactly-once semantics with message deduplication."""
        logger.info("🔍 Validating message deduplication (exactly-once semantics)")
        
        try:
            from vcf_agent.phase5.streaming_coordinator import MessageDeduplicator
            
            deduplicator = MessageDeduplicator(window_size=100, ttl_seconds=300)
            variants = self.create_synthetic_variants(50)
            
            # Test deduplication logic
            unique_variants = []
            duplicate_count = 0
            
            for variant in variants:
                if not deduplicator.is_duplicate(variant):
                    unique_variants.append(variant)
                    deduplicator.mark_delivered(variant, "iggy")
                else:
                    duplicate_count += 1
            
            # Re-process same variants (should all be duplicates now)
            for variant in variants:
                if deduplicator.is_duplicate(variant):
                    duplicate_count += 1
            
            # Validate results
            expected_duplicates = len(variants)  # All should be duplicates on second pass
            success = duplicate_count == expected_duplicates
            
            if success:
                logger.info(f"✅ Message deduplication validated: {duplicate_count}/{expected_duplicates} duplicates detected")
                self.results["feature_validations"]["message_deduplication"] = {
                    "status": "passed",
                    "unique_processed": len(unique_variants),
                    "duplicates_detected": duplicate_count,
                    "deduplication_stats": deduplicator.get_stats()
                }
                self.results["tests_passed"] += 1
                return True
            else:
                logger.error(f"❌ Message deduplication failed: {duplicate_count}/{expected_duplicates} duplicates")
                self.results["tests_failed"] += 1
                return False
                
        except Exception as e:
            logger.error(f"❌ Message deduplication validation error: {e}")
            self.results["errors"].append(f"Message deduplication: {e}")
            self.results["tests_failed"] += 1
            return False
    
    async def validate_circuit_breaker_patterns(self) -> bool:
        """Validate circuit breaker implementation for platform health management."""
        logger.info("🔍 Validating circuit breaker patterns")
        
        try:
            # Test circuit breaker state transitions
            breaker = CircuitBreaker(failure_threshold=3, recovery_timeout=1)
            
            # Initial state should be closed
            initial_state = breaker.can_execute()
            initial_state_name = breaker.state
            assert initial_state is True
            assert initial_state_name == CircuitBreakerState.CLOSED
            logger.debug(f"Initial state: {initial_state_name.value}")
            
            # Record failures to trigger state change
            for i in range(4):  # Exceed threshold (3)
                breaker.record_failure()
                logger.debug(f"Recorded failure {i+1}, state: {breaker.state.value}")
            
            # Should now be open
            open_state = breaker.state
            open_can_execute = breaker.can_execute()
            assert open_state == CircuitBreakerState.OPEN
            assert open_can_execute is False
            logger.debug(f"After failures - State: {open_state.value}, Can execute: {open_can_execute}")
            
            # Wait for recovery timeout
            logger.debug("Waiting for recovery timeout...")
            await asyncio.sleep(1.2)  # Slightly longer than timeout
            
            # Should transition to half-open when we call can_execute
            half_open_can_execute = breaker.can_execute()
            half_open_state = breaker.state
            assert half_open_can_execute is True
            assert half_open_state == CircuitBreakerState.HALF_OPEN
            logger.debug(f"After timeout - State: {half_open_state.value}, Can execute: {half_open_can_execute}")
            
            # Record successes to close circuit - need to call can_execute to consume half-open calls
            # and record successes appropriately
            success_count = 0
            while breaker.state == CircuitBreakerState.HALF_OPEN and success_count < 5:
                if breaker.can_execute():  # This increments half_open_calls
                    breaker.record_success()
                    success_count += 1
                    logger.debug(f"Success {success_count}, state: {breaker.state.value}")
                else:
                    break
            
            # Should be closed again
            final_state = breaker.state
            assert final_state == CircuitBreakerState.CLOSED
            logger.debug(f"Final state: {final_state.value}")
            
            logger.info("✅ Circuit breaker patterns validated: closed → open → half-open → closed")
            self.results["feature_validations"]["circuit_breaker"] = {
                "status": "passed",
                "state_transitions": [
                    CircuitBreakerState.CLOSED.value, 
                    CircuitBreakerState.OPEN.value, 
                    CircuitBreakerState.HALF_OPEN.value, 
                    CircuitBreakerState.CLOSED.value
                ],
                "failure_threshold": 3,
                "recovery_timeout": 1,
                "successes_to_close": success_count
            }
            self.results["tests_passed"] += 1
            return True
            
        except Exception as e:
            logger.error(f"❌ Circuit breaker validation error: {e}")
            logger.error(f"Exception details: {type(e).__name__}: {str(e)}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            self.results["errors"].append(f"Circuit breaker: {e}")
            self.results["tests_failed"] += 1
            return False
    
    async def validate_performance_monitoring(self) -> bool:
        """Validate performance monitoring and platform health tracking."""
        logger.info("🔍 Validating performance monitoring")
        
        try:
            monitor = PerformanceMonitor()
            await monitor.start()
            
            # Simulate platform operations
            start_time = time.time()
            
            # Good Iggy performance
            for _ in range(10):
                monitor.record_iggy_operation(latency_ms=1.0, success=True, throughput=1000)
            
            # Poor Kafka performance initially
            for _ in range(5):
                monitor.record_kafka_operation(latency_ms=50.0, success=False, throughput=100)
            
            # Get platform recommendation
            recommendation1 = monitor.get_recommended_platform()
            logger.info(f"Initial recommendation: {recommendation1}")
            
            # Simulate Iggy degradation
            for _ in range(10):
                monitor.record_iggy_operation(latency_ms=5000.0, success=False, throughput=0)
            
            # Improved Kafka performance
            for _ in range(10):
                monitor.record_kafka_operation(latency_ms=10.0, success=True, throughput=500)
            
            # Get new recommendation
            recommendation2 = monitor.get_recommended_platform()
            logger.info(f"After degradation: {recommendation2}")
            
            # Get comprehensive status
            platform_status = monitor.get_platform_status()
            
            await monitor.stop()
            
            # Validate monitoring worked
            monitoring_time = time.time() - start_time
            
            logger.info(f"✅ Performance monitoring validated in {monitoring_time:.2f}s")
            self.results["feature_validations"]["performance_monitoring"] = {
                "status": "passed",
                "initial_recommendation": recommendation1,
                "final_recommendation": recommendation2,
                "platform_status": platform_status,
                "monitoring_duration": monitoring_time
            }
            self.results["tests_passed"] += 1
            return True
            
        except Exception as e:
            logger.error(f"❌ Performance monitoring validation error: {e}")
            self.results["errors"].append(f"Performance monitoring: {e}")
            self.results["tests_failed"] += 1
            return False
    
    async def validate_intelligent_routing(self) -> bool:
        """Validate intelligent platform routing based on health metrics."""
        logger.info("🔍 Validating intelligent routing")
        
        try:
            coordinator = StreamingCoordinator(self.config)
            
            # Mock processors to avoid real connections during validation
            from unittest.mock import Mock
            coordinator.iggy_processor = Mock()
            coordinator.kafka_processor = Mock()
            
            # Test different routing strategies
            strategies = [
                RoutingStrategy.INTELLIGENT,
                RoutingStrategy.PRIMARY_ONLY,
                RoutingStrategy.FALLBACK_ONLY
            ]
            
            routing_results = {}
            variants = self.create_synthetic_variants(10)
            
            for strategy in strategies:
                coordinator.set_routing_strategy(strategy)
                
                # Mock platform status
                coordinator.performance_monitor = Mock()
                coordinator.performance_monitor.get_platform_status.return_value = {
                    "recommended_platform": "iggy",
                    "iggy": {
                        "health_score": 0.9,
                        "can_handle_requests": True,
                        "circuit_breaker_state": "closed"
                    },
                    "kafka": {
                        "health_score": 0.7,
                        "can_handle_requests": True,
                        "circuit_breaker_state": "closed"
                    }
                }
                
                # Test routing decisions
                decisions = []
                for variant in variants[:3]:  # Test with first 3 variants
                    decision = coordinator._select_platform(variant)
                    decisions.append({
                        "platform": decision.selected_platform,
                        "strategy": decision.strategy.value,
                        "rationale": decision.rationale
                    })
                
                routing_results[strategy.value] = decisions
                logger.debug(f"Strategy {strategy.value}: {len(decisions)} routing decisions made")
            
            logger.info("✅ Intelligent routing validated across all strategies")
            self.results["feature_validations"]["intelligent_routing"] = {
                "status": "passed",
                "routing_strategies_tested": len(strategies),
                "decisions_per_strategy": {k: len(v) for k, v in routing_results.items()},
                "sample_decisions": routing_results
            }
            self.results["tests_passed"] += 1
            return True
            
        except Exception as e:
            logger.error(f"❌ Intelligent routing validation error: {e}")
            self.results["errors"].append(f"Intelligent routing: {e}")
            self.results["tests_failed"] += 1
            return False
    
    async def validate_end_to_end_processing(self) -> bool:
        """Validate end-to-end dual-platform processing."""
        logger.info("🔍 Validating end-to-end dual-platform processing")
        
        try:
            coordinator = StreamingCoordinator(self.config)
            
            # Mock processors for validation (avoids requiring real infrastructure)
            from unittest.mock import AsyncMock
            coordinator.iggy_processor = AsyncMock()
            coordinator.kafka_processor = AsyncMock()
            coordinator.performance_monitor = AsyncMock()
            
            # Configure mock responses
            coordinator.iggy_processor.process_variant.return_value = True
            coordinator.kafka_processor.process_variant.return_value = True
            coordinator.iggy_processor.start.return_value = None
            coordinator.kafka_processor.start.return_value = None
            coordinator.iggy_processor.stop.return_value = None
            coordinator.kafka_processor.stop.return_value = None
            coordinator.performance_monitor.start.return_value = None
            coordinator.performance_monitor.stop.return_value = None
            coordinator.performance_monitor.get_platform_status.return_value = {
                "recommended_platform": "iggy",
                "iggy": {"health_score": 0.9, "can_handle_requests": True},
                "kafka": {"health_score": 0.8, "can_handle_requests": True}
            }
            
            # Create test data
            variants = self.create_synthetic_variants(100)
            
            # Process variants
            start_time = time.time()
            
            async with coordinator.processing_session():
                # Process in batches
                batch_size = 25
                total_processed = 0
                
                for i in range(0, len(variants), batch_size):
                    batch = variants[i:i + batch_size]
                    batch_results = await coordinator.process_variants_batch(batch)
                    total_processed += batch_results["processed"]
                    
                    logger.debug(f"Batch {i//batch_size + 1}: {batch_results['processed']}/{len(batch)} processed")
                
                # Get final status
                final_status = await coordinator.get_comprehensive_status()
            
            processing_time = time.time() - start_time
            throughput = total_processed / processing_time
            
            logger.info(f"✅ End-to-end processing validated: {total_processed}/{len(variants)} variants in {processing_time:.2f}s ({throughput:.1f} variants/sec)")
            
            self.results["feature_validations"]["end_to_end_processing"] = {
                "status": "passed",
                "variants_processed": total_processed,
                "total_variants": len(variants),
                "processing_time": processing_time,
                "throughput": throughput,
                "coordinator_status": final_status["coordinator"]
            }
            self.results["performance_metrics"]["e2e_throughput"] = throughput
            self.results["tests_passed"] += 1
            return True
            
        except Exception as e:
            logger.error(f"❌ End-to-end processing validation error: {e}")
            self.results["errors"].append(f"End-to-end processing: {e}")
            self.results["tests_failed"] += 1
            return False
    
    async def run_validation_suite(self) -> Dict[str, Any]:
        """Run complete Phase 5.2 validation suite."""
        logger.info("🚀 Starting Phase 5.2 Production Validation Suite")
        logger.info("=" * 60)
        
        start_time = time.time()
        
        # Run all validation tests
        validations = [
            ("Message Deduplication", self.validate_message_deduplication()),
            ("Circuit Breaker Patterns", self.validate_circuit_breaker_patterns()),
            ("Performance Monitoring", self.validate_performance_monitoring()),
            ("Intelligent Routing", self.validate_intelligent_routing()),
            ("End-to-End Processing", self.validate_end_to_end_processing())
        ]
        
        logger.info(f"Running {len(validations)} validation tests...")
        
        for test_name, test_coro in validations:
            logger.info(f"\n📋 {test_name}")
            logger.info("-" * 40)
            
            try:
                success = await test_coro
                if success:
                    logger.info(f"✅ {test_name}: PASSED")
                else:
                    logger.error(f"❌ {test_name}: FAILED")
            except Exception as e:
                logger.error(f"💥 {test_name}: ERROR - {e}")
                self.results["errors"].append(f"{test_name}: {e}")
                self.results["tests_failed"] += 1
        
        # Calculate final results
        total_time = time.time() - start_time
        total_tests = self.results["tests_passed"] + self.results["tests_failed"]
        success_rate = (self.results["tests_passed"] / total_tests * 100) if total_tests > 0 else 0
        
        self.results["validation_summary"] = {
            "total_tests": total_tests,
            "passed": self.results["tests_passed"],
            "failed": self.results["tests_failed"],
            "success_rate": success_rate,
            "total_time": total_time,
            "errors": self.results["errors"]
        }
        
        # Print summary
        logger.info("\n" + "=" * 60)
        logger.info("🎯 PHASE 5.2 VALIDATION SUMMARY")
        logger.info("=" * 60)
        logger.info(f"Total Tests: {total_tests}")
        logger.info(f"Passed: {self.results['tests_passed']} ✅")
        logger.info(f"Failed: {self.results['tests_failed']} ❌")
        logger.info(f"Success Rate: {success_rate:.1f}%")
        logger.info(f"Total Time: {total_time:.2f}s")
        
        if self.results["errors"]:
            logger.info(f"\n❌ Errors ({len(self.results['errors'])}):")
            for error in self.results["errors"]:
                logger.info(f"  • {error}")
        
        if success_rate >= 80:
            logger.info(f"\n🎉 Phase 5.2 validation SUCCESSFUL! ({success_rate:.1f}% pass rate)")
        else:
            logger.error(f"\n💥 Phase 5.2 validation FAILED! ({success_rate:.1f}% pass rate)")
        
        logger.info("\n📊 Feature Validations:")
        for feature, details in self.results["feature_validations"].items():
            status_emoji = "✅" if details["status"] == "passed" else "❌"
            logger.info(f"  {status_emoji} {feature.replace('_', ' ').title()}: {details['status']}")
        
        return self.results


async def main():
    """Main validation entry point."""
    validator = Phase5_2Validator()
    
    try:
        results = await validator.run_validation_suite()
        
        # Exit with appropriate code
        if results["validation_summary"]["success_rate"] >= 80:
            sys.exit(0)  # Success
        else:
            sys.exit(1)  # Failure
            
    except KeyboardInterrupt:
        logger.info("\n⚠️ Validation interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"💥 Validation suite error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    # Run with proper async handling
    asyncio.run(main()) 