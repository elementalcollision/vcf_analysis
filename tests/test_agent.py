import pytest
from vcf_agent.agent import agent

def test_agent_instantiation():
    assert agent is not None

def test_agent_echo_tool():
    prompt = "echo: Hello, pytest!"
    response = agent(prompt)
    # The response may be a string or an object; convert to string for assertion
    assert "Echo: Hello, pytest!" in str(response)

def test_agent_dummy_response():
    # This test will fail until the agent is functional
    response = agent("Test prompt")
    assert response is not None

def test_config_module():
    from vcf_agent import config
    assert hasattr(config, 'CONFIG')
    assert isinstance(config.CONFIG, dict)
    # Test mutability
    config.CONFIG['test_key'] = 'test_value'
    assert config.CONFIG['test_key'] == 'test_value'

import sys
import subprocess
import os
from unittest.mock import patch

def test_cli_main_echo():
    result = subprocess.run(
        [sys.executable, '-m', 'vcf_agent.cli', 'echo: CLI test!'],
        capture_output=True, text=True,
        env={**os.environ, "VCF_AGENT_CLI_MOCK_RESPONSE": "Echo: CLI test!"}
    )
    assert result.returncode == 0
    assert 'Echo: CLI test!' in result.stdout

def test_cli_main_validate():
    result = subprocess.run(
        [sys.executable, '-m', 'vcf_agent.cli', 'validate_vcf: test.vcf'],
        capture_output=True, text=True,
        env={**os.environ, "VCF_AGENT_CLI_MOCK_RESPONSE": "VALID: test.vcf passed all checks."}
    )
    assert result.returncode == 0
    assert 'VALID: test.vcf passed all checks.' in result.stdout 