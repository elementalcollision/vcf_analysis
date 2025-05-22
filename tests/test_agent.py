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