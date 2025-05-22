import os
import yaml
import json
import pytest
import re
from jsonschema import validate, ValidationError
from vcf_agent.agent import agent

CONTRACT_PATH = os.path.join(os.path.dirname(__file__), '../../prompts/vcf_analysis_summary_v1.yaml')

def load_contract():
    with open(CONTRACT_PATH, 'r') as f:
        return yaml.load(f, Loader=yaml.SafeLoader)

def load_golden(path):
    with open(path, 'r') as f:
        return json.load(f)

def extract_json_from_text(text):
    # Remove <think> blocks
    text = re.sub(r'<think>.*?</think>', '', text, flags=re.DOTALL)
    # Find the first {...} block
    json_pattern = re.compile(r'\{.*?\}', re.DOTALL)
    match = json_pattern.search(text)
    if match:
        json_str = match.group(0)
        # Try to repair common issues (e.g., missing closing brace)
        if json_str.count('{') > json_str.count('}'):
            json_str += '}'
        try:
            return json.loads(json_str)
        except json.JSONDecodeError:
            pass
    raise ValueError('No valid JSON object found in agent output.')

def agent_response(input_path, seed):
    # Construct the prompt according to the contract
    # Note: The current agent/model may not support seeding; document this
    contract = load_contract()
    prompt = f"vcf_analysis_summary: {input_path}"
    # If the agent supports a seed parameter, it should be passed here
    # For now, we ignore the seed (future: add to agent/model call if supported)
    response = agent(prompt)
    # Post-process to extract the first valid JSON object
    try:
        return extract_json_from_text(str(response))
    except Exception as e:
        pytest.fail(f'Agent output is not valid JSON: {e}\nOutput: {response}')

def test_vcf_analysis_summary_contract():
    contract = load_contract()
    test_case = contract['test_cases'][0]
    input_path = test_case['input']
    seed = test_case['seed']
    expected_output_path = test_case['expected_output']

    # Step 1: Get agent/model output
    output = agent_response(input_path, seed)

    # Step 2: Validate output against contract's JSON schema
    try:
        validate(instance=output, schema=contract['json_schema'])
    except ValidationError as e:
        pytest.fail(f'Output does not match JSON schema: {e}')

    # Step 3: Compare output to golden file
    expected = load_golden(expected_output_path)
    assert output == expected, 'Output does not match golden file'

    # Note: Determinism with seed is not enforced unless the agent/model supports it 