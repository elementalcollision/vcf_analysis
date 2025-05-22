# Prompt Contracts for VCF Analysis Agent

## Overview
This directory contains versioned, auditable prompt contracts for AI-driven VCF analysis. Each contract is designed for determinism, reproducibility, and regulatory compliance.

## Best Practices Incorporated
- **Determinism**: Prompts and test cases are designed to yield identical outputs for the same input and seed.
- **Versioning**: Each contract includes an ID, version, and changelog for traceability.
- **Testability**: Contracts specify test cases and expected outputs for automated validation.
- **Compliance**: Contracts are annotated with relevant standards (e.g., VCFv4.3, FDA, GCP).
- **Documentation**: All fields are explained for transparency and auditability.

---

## Prompt Contract: VCF Summary Analysis
- **Filename:** `vcf_analysis_summary_v1.yaml`
- **Context:** Summarizes a VCF file for variant quality and compliance.
- **Test Case:** `sample_data/HG00098.vcf.gz` → `tests/golden/vcf_summary_HG00098.json`
- **Schema:** Variant count, types, quality metrics, non-compliant records.

## Prompt Contract: VCF Comparison
- **Filename:** `vcf_comparison_v1.yaml`
- **Context:** Compares two VCF files for concordance, discordance, and quality metrics.
- **Test Case:** `sample_data/HG00098.vcf.gz` vs. `sample_data/HG00099.vcf.gz` → `tests/golden/vcf_comparison_HG00098_HG00099.json`
- **Schema:** Concordant/discordant variant counts, unique variants per file, quality metrics.

## Prompt Contract: VCF Summarization
- **Filename:** `vcf_summarization_v1.yaml`
- **Context:** Summarizes a VCF file for key variant statistics and sample-level insights.
- **Test Case:** `sample_data/HG00098.vcf.gz` → `tests/golden/vcf_summarization_HG00098.json`
- **Schema:** Variant count, types, sample statistics (mean depth, het ratio), notable patterns.

---

### Field Explanations (applies to all contracts)
- **id**: Unique identifier for the prompt contract.
- **version**: Semantic version for tracking changes.
- **compliance**: List of standards or regulations addressed.
- **changelog**: History of changes for auditability.
- **role**: The AI's persona or expertise for context.
- **context**: Background and scope of the analysis.
- **instructions**: Step-by-step tasks for the AI.
- **constraints**: Hard requirements (e.g., output format, data usage).
- **evaluation**: Criteria for output validation and determinism.
- **test_cases**: Inputs, seeds, and expected outputs for automated testing.
- **json_schema**: Output schema for validation and downstream integration.

## Review Guidance
- Ensure all fields are present and clearly explained.
- Confirm that test cases and output schemas are sufficient for reproducibility.
- Check that compliance and changelog fields are up to date.

---

## Usage Examples

### VCF Summary Analysis
```python
from vcf_agent.agent import agent
prompt = "vcf_analysis_summary: sample_data/HG00098.vcf.gz"
response = agent(prompt)
print(response)  # Should be valid JSON matching the contract schema
```

### VCF Comparison
```python
from vcf_agent.agent import agent
prompt = "vcf_comparison: sample_data/HG00098.vcf.gz sample_data/HG00099.vcf.gz"
response = agent(prompt)
print(response)  # Should be valid JSON matching the contract schema
```

### VCF Summarization
```python
from vcf_agent.agent import agent
prompt = "vcf_summarization: sample_data/HG00098.vcf.gz"
response = agent(prompt)
print(response)  # Should be valid JSON matching the contract schema
```

---

## Reviewer Checklist
- [ ] All required fields (id, version, compliance, changelog, role, context, instructions, constraints, evaluation, test_cases, json_schema) are present and clearly explained.
- [ ] Instructions explicitly require ONLY valid JSON output, with no extra text or formatting.
- [ ] Test cases are provided and reference golden outputs for reproducibility.
- [ ] Output schema (json_schema) is complete and matches expected downstream usage.
- [ ] Changelog and versioning are up to date and meaningful.
- [ ] Compliance tags are relevant to the analysis context (e.g., VCFv4.3, FDA, GCP).
- [ ] Usage examples are provided for each contract.
- [ ] Automated tests exist and validate both schema and golden output.
- [ ] Documentation is clear for both internal and external reviewers.

---
For questions or suggestions, please open an issue or submit a pull request. 