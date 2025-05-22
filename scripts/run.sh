#!/bin/bash
# Run the VCF Analysis Agent container (OrbStack/Docker)
set -e

docker run --rm -it \
  --name vcf_agent_dev \
  -v "$(pwd)":/app \
  -p 8000:8000 \
  -e ORBSTACK_SERVER=192.168.0.225 \
  vcf_agent_dev bash 