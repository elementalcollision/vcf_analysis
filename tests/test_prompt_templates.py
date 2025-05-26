import os
import sys
import pytest
import yaml

# Ensure src/ is in the path for import
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from vcf_agent import prompt_templates

PROMPT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../prompts'))
SAMPLE_DIR = os.path.join(os.path.dirname(__file__), '..', 'sample_test_data')

MINIMAL_VCF = os.path.join(SAMPLE_DIR, 'edgecase_minimal.vcf')
MULTIALLELIC_VCF = os.path.join(SAMPLE_DIR, 'edgecase_multiallelic.vcf')


def test_load_prompt_contract_summarization():
    contract = prompt_templates.load_prompt_contract('vcf_summarization_v1')
    assert isinstance(contract, dict)
    assert contract['id'].startswith('vcf-summarization')
    assert 'instructions' in contract


def test_load_prompt_contract_comparison():
    contract = prompt_templates.load_prompt_contract('vcf_comparison_v1')
    assert isinstance(contract, dict)
    assert contract['id'].startswith('vcf-comparison')
    assert 'instructions' in contract


def test_render_prompt_summarization():
    contract = prompt_templates.load_prompt_contract('vcf_summarization_v1')
    prompt = prompt_templates.render_prompt(contract, [MINIMAL_VCF])
    assert 'VCF file:' in prompt
    assert 'Instructions:' in prompt
    assert contract['instructions'].split('\n')[0] in prompt


def test_render_prompt_comparison():
    contract = prompt_templates.load_prompt_contract('vcf_comparison_v1')
    prompt = prompt_templates.render_prompt(contract, [MINIMAL_VCF, MULTIALLELIC_VCF])
    assert 'VCF file 1:' in prompt
    assert 'VCF file 2:' in prompt
    assert 'Instructions:' in prompt
    assert contract['instructions'].split('\n')[0] in prompt


def test_get_prompt_for_task():
    prompt = prompt_templates.get_prompt_for_task('vcf_summarization_v1', [MINIMAL_VCF])
    assert 'VCF file:' in prompt
    prompt2 = prompt_templates.get_prompt_for_task('vcf_comparison_v1', [MINIMAL_VCF, MULTIALLELIC_VCF])
    assert 'VCF file 1:' in prompt2
    assert 'VCF file 2:' in prompt2


def test_load_prompt_contract_missing():
    with pytest.raises(FileNotFoundError):
        prompt_templates.load_prompt_contract('not_a_real_contract')


def test_load_prompt_contract_empty(tmp_path):
    # Create an empty YAML file
    empty_path = tmp_path / 'empty_contract.yaml'
    empty_path.write_text('')
    # Patch PROMPT_DIR to tmp_path for this test
    orig_dir = prompt_templates.PROMPT_DIR
    prompt_templates.PROMPT_DIR = str(tmp_path)
    try:
        with pytest.raises(ValueError):
            prompt_templates.load_prompt_contract('empty_contract')
    finally:
        prompt_templates.PROMPT_DIR = orig_dir 