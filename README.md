# VCF Analysis Agent

## Overview
The VCF Analysis Agent is an automated system for processing, analyzing, and extracting insights from Variant Call Format (VCF) and BCF files. It leverages the Strands Agent SDK, bcftools, and AI models (OpenAI, Cerebras) for comparative genomic analysis, with a focus on reproducibility, observability, and extensibility.

## Key Features
- Automated VCF/BCF parsing, annotation, filtering, and comparison
- Modular bcftools integration
- AI-powered comparative analysis (OpenAI, Cerebras)
- Containerized and orchestrated with Kestra, OrbStack, and Docker
- Comprehensive metrics and observability (Prometheus, Jaeger, Grafana)
- Extensible, testable, and reproducible workflows

## Quickstart
1. **Clone the repository**
2. **Set up Python environment** (see `.context/tasks/planned/TASK-001-02.md`)
3. **Build and run container** (see `.context/tasks/planned/TASK-001-03.md`)
4. **Configure and run agent** (see `.context/tasks/planned/TASK-001-05.md`)
5. **Run tests** (`pytest`)

See the `.context/plan/planning_document.md` and PRD for full requirements and architecture.

## Documentation
- [Product Requirements Document (PRD)](PRD%20-%20%20VCF%20Analysis%20Agent.md)
- [Aegis Framework Structure](.context/AI_INSTRUCTIONS.md)
- [Task Breakdown](.context/tasks/planned/)

## License
[MIT](LICENSE) 