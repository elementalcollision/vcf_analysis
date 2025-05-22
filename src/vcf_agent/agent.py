from strands import Agent, tools
from strands.tools import tool
from strands.models.ollama import OllamaModel
from typing import Any

SYSTEM_PROMPT = "You are the VCF Analysis Agent. You analyze, validate, and process VCF files for genomics workflows."

# Placeholder for tool imports and configuration
# from strands_tools import calculator, file_read, shell

# Placeholder echo tool using @tools.tool decorator
@tools.tool
def echo(text: str) -> str:
    """
    Echoes the input text back to the user.

    Args:
        text: Text to echo back.
    
    Returns:
        The echoed text.
    """
    return f"Echo: {text}"

# Register the decorated tool directly
tools_list = [echo]

# Configure Ollama model
ollama_model = OllamaModel(
    host="http://192.168.0.225:11434",
    model_id="qwen3:latest"
)

agent = Agent(
    system_prompt=SYSTEM_PROMPT,
    tools=tools_list,  # type: ignore # Suppress linter warning as runtime is OK
    model=ollama_model,
    callback_handler=None
)

__all__ = ["agent", "SYSTEM_PROMPT"] 