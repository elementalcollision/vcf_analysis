import os
import sys
import pytest
import json
from unittest.mock import patch, MagicMock

# Ensure src/ is in the path for import
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from vcf_agent.agent import run_llm_analysis_task
from vcf_agent.config import SessionConfig

SAMPLE_DIR = os.path.join(os.path.dirname(__file__), '..', 'sample_test_data')
MINIMAL_VCF = os.path.join(SAMPLE_DIR, 'minimal.vcf.gz')

@pytest.mark.parametrize("provider", ["ollama", "openai", "cerebras"])
@patch('vcf_agent.agent.get_agent_with_session')
def test_run_llm_analysis_task_mock(mock_get_agent, monkeypatch, provider):
    # Configure the mock agent to return a specific response when called
    mock_agent_instance = MagicMock()
    mock_response_content = {"mock": True, "provider": provider, "task": "vcf_summarization_v1"}
    # The agent is a callable, so we set its return_value for when it's called with the prompt
    mock_agent_instance.return_value = json.dumps(mock_response_content) 
    mock_get_agent.return_value = mock_agent_instance

    session_config = SessionConfig(model_provider=provider)
    
    response_str = run_llm_analysis_task(
        task="vcf_summarization_v1",
        file_paths=[MINIMAL_VCF],
        session_config=session_config,
        model_provider=provider
    )
    
    # run_llm_analysis_task is expected to return a JSON string, parse it
    response = json.loads(response_str)

    assert isinstance(response, dict)
    if "error" in response:
        # If an error occurred, it means our mock wasn't effective as expected or there's another issue
        pytest.fail(f"Test failed with error in response: {response['error']}")
    else:
        assert response.get("mock") is True
        assert response.get("provider") == provider
        assert response.get("task") == "vcf_summarization_v1"

def test_run_llm_analysis_task_error(monkeypatch):
    # Use a non-existent contract to trigger error
    session_config = SessionConfig(model_provider="openai")
    with pytest.raises(FileNotFoundError):
        run_llm_analysis_task(
            task="not_a_real_contract",
            file_paths=[MINIMAL_VCF],
            session_config=session_config,
            model_provider="openai"
        ) 