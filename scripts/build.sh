#!/bin/bash
# Build the Docker image for VCF Analysis Agent (OrbStack/Docker)
set -e

docker build -t vcf_agent_dev -f Dockerfile . 