#!/usr/bin/env python3
"""
Quick test script to start the Prometheus metrics server and simulate AI agent interactions.
Run this to verify Prometheus can scrape /metrics.
"""
from src.vcf_agent.logging_metrics import start_metrics_server, log_ai_interaction
import time

if __name__ == "__main__":
    # Start metrics server on all interfaces
    start_metrics_server(8000)
    print("Metrics server running on :8000/metrics. Simulating AI interactions...")
    user = "test_user"
    model = "test-model"
    for i in range(5):
        prompt = f"Test prompt {i}"
        response = f"Test response {i}"
        tokens = 10 + i
        latency = 0.1 * (i + 1)
        log_ai_interaction(user, model, prompt, response, latency, tokens, success=True)
        time.sleep(1)
    print("Done. Visit /metrics to verify Prometheus scraping.")
    # Keep the server running for manual testing
    while True:
        time.sleep(10) 