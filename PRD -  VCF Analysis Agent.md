## **Product Requirements Document: VCF Analysis Agent**

Version: 1.1  
Date: May 21, 2025  
Status: Draft  
**1\. Overview & Introduction**

This document outlines the product requirements for the VCF Analysis Agent (also referred to as VCF Peeling Agent). This agent is an automated system designed to process, analyze, "peel apart," and extract insights from Variant Call Format (VCF) and BCF files. It will be built using the Strands Agent SDK to orchestrate tasks, leveraging bcftools for core genomic data manipulation, and integrating with AI models from OpenAI and Cerebras for comparative analysis and insight generation.

The project emphasizes a modern, robust development stack:

* **Python:** Version 3.11+ as the primary programming language.  
* **Package Management:** uv (e.g., v1.0+) for fast, deterministic dependency resolution.  
* **Containerization:** OrbStack (e.g., v2.4+) for local macOS development and testing, with OCI/Docker images for CI/CD and deployment.  
* **Orchestration:** Kestra (e.g., v0.15+) for defining and managing agent tasks and complex workflows.  
* **Monitoring & Observability:** Prometheus (e.g., v3.0+) for metrics collection and Jaeger (e.g., v1.50+) via OpenTelemetry (OTEL) for distributed tracing.

**2\. Goals & Objectives**

**2.1 Primary Goals**

* Develop an intelligent agent capable of parsing, annotating, filtering, and comparing VCF/BCF files with deterministic reproducibility.  
* Provide a flexible and extensible framework for integrating various bcftools operations as tools within the Strands Agent.  
* Enable side-by-side reasoning and comparative analysis of VCF data or derived insights using AI models from OpenAI (e.g., GPT-4o/GPT-4-Turbo) and Cerebras (e.g., Llama-4-Scout-17B).  
* Ensure robust, testable, and highly observable agent performance within a containerized, orchestrated environment.  
* Expose comprehensive metrics and traces suitable for defining and monitoring Service Level Objectives (SLOs), e.g., 95th-percentile latency ≤2 seconds for analysis of ≤1 MB VCF slices.

**2.2 Key Objectives**

* **Automation:** Significantly reduce manual effort in VCF data extraction, manipulation, and preliminary analysis.  
* **Flexibility:** Allow users to define custom VCF processing workflows or leverage AI for dynamic task execution and decision-making.  
* **Insight Generation:** Facilitate deeper understanding of variant data through AI-powered summarization, comparison, and potentially hypothesis generation, including confidence scoring for LLM outputs.  
* **Scalability & Efficiency:** Design the agent for efficient single-node operation, with Kestra handling broader workflow scaling for larger VCF files and higher processing loads.  
* **Reproducibility:** Guarantee consistent results through version-controlled code, containerized environments, locked dependencies, and deterministic workflow orchestration.

**2.3 Non-Goals (for Version 1.0)**

* Full clinical-grade variant interpretation or diagnostic decision-making.  
* Large-scale, on-premise cluster deployment (this may be considered in a later phase).  
* Advanced AI-driven VCF error correction or de novo variant calling.  
* A graphical user interface (GUI); interaction will be primarily command-line, API, or Kestra-driven.  
* Complex multi-agent swarm behaviors beyond the primary VCF analysis agent.  
* Direct processing of alignment files (e.g., BAM/CRAM); pre-processing to VCF/BCF is assumed.

**3\. Target Audience**

* **Bioinformaticians:** Requiring automated, configurable tools for VCF data manipulation, filtering, and initial analysis.  
* **Genomic Researchers:** Seeking to extract specific variant information, explore variant characteristics, and leverage AI for deeper insights.  
* **Computational Biologists:** Developing or utilizing automated pipelines for genomic data processing and analysis.

**4\. Technical Architecture & Specifications**

**4.1 Core Components & Technologies**

* **Agent Runtime:** Strands Agent SDK (e.g., v0.1.1+) utilizing Python 3.11+ (with type hints and asynchronous support).  
* **Variant Operations:** bcftools (e.g., v1.21+) CLI, integrated via Python shims/wrappers (potentially using pysam for robust interaction or direct subprocess calls).  
  * Key bcftools commands: view, query, filter, annotate, stats, norm (including multi-allelic site splitting), consensus.  
* **LLM Backends:**  
  * OpenAI API (e.g., GPT-4o, GPT-4-Turbo) via OpenAI Agents SDK or direct API calls.  
  * Cerebras Inference API (e.g., Llama-4-Scout-17B) via appropriate client/SDK.  
* **Package Management:** uv (e.g., v1.0+).  
* **Containerization:**  
  * OrbStack (e.g., v2.4+) on macOS for local development, ensuring fast iteration and x86/ARM emulation.  
  * OCI images / Docker for CI (e.g., GitHub Actions on Ubuntu) and deployment (Python 3.11-slim base image preferred).  
* **Orchestration / CI:** Kestra (e.g., v0.15+) task-runners for workflow management.  
* **Monitoring & Observability:**  
  * Prometheus (e.g., v3.0+) for metrics collection (prometheus\_client library).  
  * Jaeger (e.g., v1.50+) for distributed tracing via OpenTelemetry (OTEL SDK, OTLP exporter).  
  * Grafana for dashboarding Prometheus metrics.  
* **Data Handling:**  
  * **Input:** Local or network-accessible VCF/BCF files (uncompressed and bgzipped, with CSI index auto-detection).  
  * **Output:** New or modified VCF/BCF files, text-based reports (e.g., variant lists, statistics), annotated per-variant JSON, structured logs containing results from AI analysis. Outputs written to a configurable directory or passed to downstream Kestra tasks.  
  * **Intermediate Data:** Managed within the container's ephemeral storage or a Kestra-designated working directory.  
* **Analytical Storage (for specific features/milestones):**  
  * **Vector Database:** LanceDB (Arrow-native) for storing variant embeddings (e.g., 384-dimensions) and enabling similarity searches. Python client SDK integration.  
  * **Graph Database:** Neo4j for storing gene-variant relationship graphs (e.g., :Variant, :Sample nodes with indexed properties like chr, pos, ref, alt). Cypher query API integration.

**4.2 Agent Logic & Design**

* **Strands Agent SDK:** Serves as the central framework for defining the agent's behavior, managing tools (wrapped bcftools commands, custom Python functions), and orchestrating LLM interactions.  
* **Prompt Engineering:** Develop clear, effective, and deterministic prompt contracts to guide VCF analysis tasks and ensure consistent AI model responses and comparisons.  
* **Tool Definition:** Encapsulate bcftools commands and custom Python data processing functions as Strands tools, allowing for dynamic invocation and parameterization.  
* **Comparative AI Analysis:** Design specific tasks where the agent can send similar queries or processed data to both OpenAI and Cerebras models, then present a structured comparison of their responses, potentially including confidence scores or justifications for differences.

**5\. Functional Requirements**

* **F1: VCF/BCF File Ingestion and Validation:**  
  * Accept one or more VCF/BCF file paths as input.  
  * Automatically detect and handle bgzip compression and associated CSI indexes.  
  * Perform basic validation (e.g., file existence, format integrity check using bcftools stats or validate). Report errors clearly for malformed or missing files.  
* **F2: Configurable bcftools-based VCF Operations ("Peeling"):**  
  * Execute a predefined and extensible set of bcftools commands via Strands tools, with parameters passed dynamically.  
  * Essential operations include:  
    * **Filtering:** Based on INFO, FORMAT fields, quality scores, genomic regions, PASS status, etc. (bcftools view \-i '...').  
    * **Extraction:** Specific fields or annotations into custom formats (bcftools query \-f '...').  
    * **Subsetting:** By samples or genomic regions.  
    * **Statistics:** Generation of summary statistics (bcftools stats).  
    * **Normalization:** Variant normalization, including splitting multi-allelic sites (bcftools norm).  
    * **Annotation:** Applying annotations from external files (e.g., BED, VCF) (bcftools annotate).  
    * **Consensus:** Generation of consensus sequences (bcftools consensus).  
  * Ensure bcftools output (stdout, stderr, generated files) is correctly captured and made available for subsequent steps or reporting.  
* **F3: Variant Analysis and Comparison using AI:**  
  * Format and send selected VCF data snippets, variant summaries, or bcftools outputs to Cerebras and/or OpenAI models for tasks like interpretation, summarization, or comparison.  
  * Example AI tasks:  
    * "Summarize the functional implications of these top N variants: \[variant data\]."  
    * "Compare the types of variants filtered by method A vs. method B, based on these summary statistics: \[stats A\], \[stats B\]."  
    * "Given this VCF excerpt for gene X, what are potential areas of interest or concern for further investigation?"  
    * Generate LLM comparison commentary with a confidence score.  
  * Successfully query both Cerebras and OpenAI APIs, handling authentication and formatting data appropriately for LLM prompts.  
  * Capture, log, and structure AI-generated responses for reporting and downstream use.  
* **F4: Reporting and Output Generation:**  
  * Produce annotated per-variant JSON.  
  * Generate processed VCF/BCF files as specified by the workflow.  
  * Create text-based summaries or reports as configured.  
  * Maintain structured logs (e.g., JSON) of all agent actions, bcftools commands executed, parameters used, and AI API calls/responses.  
* **F5: Data Persistence for Advanced Analytics (as per milestones):**  
  * Store variant embeddings (e.g., 384-dimensions) generated from VCF data into a LanceDB table.  
  * Ingest relevant data into Neo4j, creating nodes (e.g., :Variant, :Sample) with indexed properties (chr, pos, ref, alt) and establishing relationships (e.g., sample\_has\_variant).

**6\. Non-Functional Requirements**

* **Performance:**  
  * End-to-end pipeline processing: Target ≤60 seconds for a 10 MB VCF on representative hardware (e.g., M2 MacBook Pro).  
  * Agent-level SLO: 95th-percentile latency for core analysis tasks on ≤1 MB VCF slices should be ≤2 seconds.  
* **Resource Limits:** Target ≤4 GB RAM and 2 vCPUs per agent execution container under typical load.  
* **Portability:** Ensure the agent runs consistently in OrbStack on macOS (for local development) and in GitHub Actions CI environments (Ubuntu-based Docker).  
* **Reproducibility:** Guarantee deterministic outputs for identical inputs and configurations.  
* **Testability:** Code must be designed for high unit and integration test coverage (target ≥90 for critical modules).  
* **Security:** Adhere to best practices for API key management and data handling. Include a SECURITY.md file.

**7\. Milestones & Timeline (Estimated Total: \~12-14 Weeks)**

This timeline includes buffer for research, unexpected issues, iterative feedback, and approval windows.

* **Phase 0: Foundation & Scaffolding** (Est. 1-2 Weeks)  
  * Deliverables: Project repository setup with uv for environment and dependency management (pinned dependencies, lockfile). Base Python 3.11 OrbStack container definition. Initial Kestra workflow for CI/CD pipeline. Basic Strands agent scaffolding with placeholder tool definitions and dual LLM routing logic.  
  * Client Gate: Kick-off Review / Architecture Review.  
* **Phase 1: Core VCF Processing Engine** (Est. 2-3 Weeks)  
  * Deliverables: Python wrappers/tools for key bcftools commands (view, query, filter, norm, stats, annotate). Robust VCF/BCF file I/O and validation. Comprehensive pytest unit tests for bcftools wrappers (target ≥90 coverage), including golden-file diffs and SAMspec compliance checks.  
  * Client Gate: Technical Demo of bcftools integration.  
* **Phase 2: Strands Agent & AI Integration** (Est. 3-4 Weeks)  
  * Deliverables: Integration of bcftools tools into the Strands Agent. Development of prompt contracts. Connection to OpenAI (GPT-4o/Turbo) and Cerebras (Llama-4-Scout-17B) APIs. Initial AI analysis tasks (e.g., summarization, comparison logic) implemented. Logging of AI interactions. Unit tests mocking LLM calls and validating prompt formatting. Cost/performance metrics system for LLMs initiated.  
  * Client Gate: Phase-1 UAT / Internal Demo of AI capabilities.  
* **Phase 3: Orchestration, Containerization & Observability Setup** (Est. 2-3 Weeks)  
  * Deliverables: Dockerization of the agent. Development of Kestra flows for agent execution, parameterization, and output handling. Testing agent within Kestra environment (local via OrbStack, then target CI). Setup of Prometheus for metrics collection (agent performance, bcftools times, AI latencies, token usage, errors) and Jaeger for tracing agent execution flow (OTEL instrumentation for tool calls and agent decision paths). Initial Grafana/Jaeger dashboards.  
  * Client Gate: Ops Sign-off / Security Audit (preliminary).  
* **Phase 4: Data Stores & End-to-End Testing** (Est. 2-3 Weeks)  
  * Deliverables: Integration with LanceDB for variant embedding storage and Neo4j for graph ingest (sample → variant edges with Cypher loaders). Kestra E2E workflow tests covering file input, agent processing, AI calls, and data store interactions. Load testing (e.g., to 1k variants/sec or representative workload). pytest-memray for memory leak detection.  
  * Client Gate: Data Engineering Review / Phase-2 UAT / Final Client Demo.  
* **Phase 5: Hardening, Documentation & Release Prep** (Est. 1-2 Weeks)  
  * Deliverables: SECURITY.md, load test report, run-books. Developer quickstart guide (OrbStack setup, API key config, LanceDB/Neo4j schema examples). User documentation. Code refinement, optimization, and final code freeze. Potential Helm chart for Kubernetes deployment.  
  * Client Gate: Go-live Approval / Production Approval.

**8\. Testing & CI/CD Strategy**

* **Frameworks/Tools:** pytest (for unit and integration tests), pytest-cov (for coverage), hypothesis (for property-based testing of VCF parsing/manipulation), vcr.py or similar (for mocking external API calls to LLMs), pytest-memray (for memory leak detection).  
* **Environment:** All tests to run inside OrbStack locally before commits/pushes and again in the CI pipeline (GitHub Actions) using Docker containers that mirror the deployment setup.  
* **Test Types & Focus:**  
  * **bcftools Wrappers:** Golden-file diffs against bcftools stdout/output files. Verify correct command generation, output parsing, and error handling for diverse VCF inputs (including edge cases like multi-allelic sites, contig naming, header variations). Ensure SAMspec compliance.  
  * **Strands Agent Logic:** Test prompt handling, tool selection logic, and overall flow control within the agent. Mock AI and bcftools tool responses for deterministic testing.  
  * **AI Model Interaction:** Verify API connectivity, correct request formatting (prompts, parameters), and response parsing for both Cerebras and OpenAI. Mock API endpoints for reliable and cost-effective testing. Validate LLM response structure (e.g., using Jaccard similarity for text, schema validation for structured output).  
  * **Data Stores Integration:** Test LanceDB similarity searches return expected identifiers; Neo4j Cypher queries for path counts and data integrity.  
  * **Kestra Workflow Integration:** Test the invocation of the agent from Kestra tasks, parameter passing, artifact handling, and output retrieval.  
  * **Container Builds:** Use docker scan for vulnerability scanning and hadolint for Dockerfile best practices.  
  * **Data Validation:** Ensure input VCFs are handled correctly and all output formats meet specifications.  
* **Test Data:** A curated set of small, representative VCF/BCF files (including some with known errors, complex variants, or edge cases) will be maintained and version-controlled.

**9\. Monitoring & Observability Strategy**

* **Instrumentation:**  
  * **Metrics (Prometheus):** Utilize prometheus\_client in Python to expose custom counters, gauges, and histograms.  
    * Agent: Execution time (total, per stage/tool), error rates, VCF records processed, variant processing throughput.  
    * bcftools: Command execution duration, CPU/memory utilization (if feasible to extract).  
    * AI APIs: Call latency (per provider), token usage (input/output), API error rates.  
    * Kestra: Task success/failure rates, queue lengths (if applicable).  
    * Container: Resource usage (CPU, memory) via Kestra/container metrics.  
  * **Tracing (Jaeger via OpenTelemetry):**  
    * Implement distributed tracing using the OpenTelemetry SDK. Auto-instrument OpenAI calls where possible (e.g., OTEL generative-AI instrumentation for token counts). Manually instrument other tool calls and critical agent decision paths.  
    * Configure OTLP exporter to a Jaeger all-in-one image (e.g., ports 4317/4318) or a central collector.  
    * Create spans for: overall agent task execution, each bcftools operation, each AI API call, significant internal processing steps, and data store interactions.  
    * Goal: Visualize request flow, identify performance bottlenecks, and facilitate debugging across distributed components. Consider trace sampling (e.g., 1-10% in production) to manage overhead.  
* **Dashboards:**  
  * Grafana: Create dashboards sourcing data from Prometheus to visualize key performance indicators, error rates, and resource utilization.  
  * Jaeger UI: Utilize for trace exploration, filtering by service (e.g., service=vcf-agent), operation, or tags.  
* **Logging:**  
  * Implement structured logging (e.g., JSON format) for all agent activities, executed bcftools commands (with parameters), AI interactions (prompts, sanitized responses), and errors.  
  * Ensure logs are easily parsable and can be integrated with Kestra's logging system and potentially a centralized logging platform (e.g., ELK stack, Grafana Loki).  
  * Include correlation IDs (e.g., trace IDs, request IDs) in logs to link them with distributed traces and facilitate debugging.

**10\. Success Metrics**

**10.1 Quantitative**

* **Processing Speed & Throughput:** Adherence to NFRs (e.g., ≤60s for 10MB VCF, ≤2s for ≤1MB slice analysis). Measure variants processed per second/minute for benchmark datasets.  
* **Accuracy & Correctness:** \>99.9 correctness for bcftools operations (verified against manual checks or known outputs). VCF parsing accuracy (100% SAMspec compliance for supported features). High Jaccard similarity or other semantic measures for AI summarization tasks against golden datasets.  
* **AI Task Completion Rate:** \>98 successful completion of AI analysis tasks (excluding API provider outages).  
* **Resource Utilization:** CPU and memory footprint of the agent within defined limits (≤4GB RAM, 2 vCPU) under typical and peak loads.  
* **Error Rate:** Application error rate \<0.1 during agent operation under normal conditions.  
* **Test Coverage:** Achieve and maintain ≥90 unit test coverage for critical agent logic and bcftools wrappers.

**10.2 Qualitative**

* **Ease of Use:** Positive feedback from target users (bioinformaticians, researchers) on the ease of defining tasks, configuring workflows, and interpreting results.  
* **Reliability & Stability:** Consistent and stable operation of the agent and Kestra workflows over extended periods and varying loads.  
* **Usefulness of AI Insights:** Subjective assessment from users on the value, relevance, and actionability of the AI-generated analyses, summaries, and comparisons.  
* **Clarity of Monitoring Data:** Effectiveness of Prometheus/Jaeger dashboards and structured logs in understanding agent behavior, diagnosing issues, and identifying performance bottlenecks.

**11\. Risks & Mitigation Strategies**

| Risk | Likelihood | Impact | Mitigation Strategy |
| :---- | :---- | :---- | :---- |
| Upstream LLM API Changes/Rate-Limits/Performance | Medium | Medium-High | Implement a robust API interaction layer with versioning awareness. Implement retry mechanisms with exponential backoff, circuit breakers. Throttle requests; batch prompts where feasible. Monitor API documentation. Have fallback options or version pinning for SDKs. Consider local Ollama fallback for non-critical tasks if feasible. |
| VCF Complexity & Edge Cases (e.g., multi-allelic) | Medium | Medium | Utilize bcftools norm extensively for normalization. Implement a SAMtools/bcftools validate layer. Build an extensive regression test suite with diverse and problematic VCF examples. Use robust parsing libraries like pysam. |
| bcftools Integration Complexity | Medium | Medium | Favor pysam for its robust bindings if it meets performance needs; otherwise, carefully manage subprocess calls. Start with a core set of bcftools commands and expand iteratively. Thoroughly test with diverse VCFs and parameters. |
| Performance Bottlenecks with Large VCFs | Medium | High | Optimize bcftools usage (e.g., appropriate indexing, streaming where possible, regional queries). Profile agent performance continuously. Design Kestra workflows to manage batching or chunking of large files if necessary. |
| Difficulty in Effective Prompt Engineering | Medium | Medium | Employ an iterative prompt development process with rigorous testing and versioning. Allow users to customize or select from pre-defined prompt templates if appropriate. Start with specific, narrowly defined AI tasks and expand. |
| Container Environment Divergence (x86/ARM) | Low-Medium | Medium | Leverage OrbStack's multi-arch testing capabilities during local development. Standardize on multi-arch Docker base images. Perform CI tests on relevant architectures if feasible. |
| OpenTelemetry (OTEL) Instrumentation Overhead | Low | Low-Medium | Start with essential instrumentation. Monitor performance impact. Configure trace sampling (e.g., 1-10% of traces in production) if significant overhead is observed. |
| Kestra Integration Challenges | Low | Medium | Utilize Kestra documentation, examples, and community support. Start with simple "hello world" style container execution tasks and incrementally build complexity. Ensure proper handling of inputs, outputs, and state. |
| Dependency Management (uv, Python versions) | Low | Low | Strictly use uv for creating and managing virtual environments and for locking dependencies (uv pip compile for requirements.txt). Clearly document setup procedures. Use Docker for final build to ensure runtime consistency. |
| API Key Security & Management | Medium | High | Utilize secure secret management solutions (e.g., HashiCorp Vault, cloud provider KMS, GitHub Actions secrets). Avoid hardcoding keys. Implement least-privilege access for API keys. |

**12\. Assumptions & Dependencies**

* Access to OpenAI and Cerebras API keys with sufficient quotas will be provisioned.  
* A Kestra environment (local or shared) will be available for development, testing, and deployment.  
* OrbStack is the standard for local macOS development; Docker is standard for CI and other environments.  
* Prometheus and Jaeger instances (or a managed observability platform supporting these) will be available or set up.  
* Input VCF/BCF files are generally syntactically valid according to their respective specifications, though the agent should handle common errors and deviations gracefully.  
* bcftools (specified version) is installable and functions correctly within the containerized Linux environment.  
* Python 3.11+ is the established baseline for the project.  
* uv is the chosen tool for Python package and project management.

**13\. Future Considerations & Potential Extensions (Post V1.0)**

* **Enhanced Annotation:** Integration with Ensembl VEP (REST API or local instance), dbSNP, ClinVar, etc.  
* **Advanced Graph Analytics:** Utilize Neo4j Graph Data Science (GDS) library for complex queries (e.g., node similarity, pathfinding, community detection).  
* **Gene Expression Integration:** Consider integration with tools like BioVDB for correlating variants with gene expression patterns.  
* **Faster Python VCF Parsing:** Explore cyvcf2 for performance-critical random access VCF parsing if pysam proves insufficient.  
* **Interactive Visualization:** Develop components or integrations for JupyterLab (e.g., pysam pile-up widgets) or other tools for ad-hoc QC and exploration.  
* **User-Defined AI Tasks:** Allow users to provide their own prompts or select from a broader library of AI analysis tasks.  
* **Support for Additional File Formats:** Potentially extend to GVF, GVCF, or direct BAM/CRAM processing if demand arises.

**14\. Client Approval & Review Process**

* **Phase 0 Review (End of Phase 0):** Review of project setup, architecture, and initial scaffolding.  
* **Phase 1 Review (End of Phase 1):** Demo of core bcftools wrapper functionality and VCF processing capabilities. Client feedback on core functionality.  
* **Phase 2 Review (End of Phase 2):** Demo of Strands Agent with integrated AI (OpenAI & Cerebras) for initial analysis tasks. Review of AI task design and comparative output.  
* **Phase 3 Review (End of Phase 3):** Demo of the agent running within Kestra, containerized, with initial observability (Prometheus/Jaeger) dashboards. Review of workflow design.  
* **Phase 4 Review & UAT Kick-off (End of Phase 4):** Demo of end-to-end workflows including data store integration. Client team begins User Acceptance Testing (UAT) with their own VCF files and use cases against pre-agreed test scenarios.  
* **Final Sign-off (End of Phase 5):** Formal approval from the client upon successful UAT completion, delivery of all documentation, and final demonstration of all features and monitoring capabilities.

This PRD is intended as a living document and will be updated via a defined change control process as the project progresses and new requirements or insights emerge.

**15\. Next Steps**

1. Review and formally approve this PRD, including the proposed scope, milestones, and timeline.  
2. Establish project boards (e.g., GitHub Projects, Jira) aligned with the defined phases/milestones.  
3. Officially commence Phase 0: Foundation & Scaffolding.