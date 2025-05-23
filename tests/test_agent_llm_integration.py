import os
import sys
import pytest

# Ensure src/ is in the path for import
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from vcf_agent.agent import run_llm_analysis_task
from vcf_agent.config import SessionConfig

SAMPLE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../sample_data'))
MINIMAL_VCF = os.path.join(SAMPLE_DIR, 'edgecase_minimal.vcf')

@pytest.mark.parametrize("provider", ["ollama", "openai", "cerebras"])
def test_run_llm_analysis_task_mock(monkeypatch, provider):
    # Patch the agent to return a mock response
    monkeypatch.setenv("VCF_AGENT_CLI_MOCK_RESPONSE", '{"mock": true, "provider": "%s"}' % provider)
    session_config = SessionConfig(model_provider=provider)
    response = run_llm_analysis_task(
        task="vcf_summarization_v1",
        file_paths=[MINIMAL_VCF],
        session_config=session_config,
        model_provider=provider
    )
    assert isinstance(response, str)
    assert 'mock' in response or 'error' in response
    assert provider in response
    # Clean up env var
    monkeypatch.delenv("VCF_AGENT_CLI_MOCK_RESPONSE", raising=False)

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