from ils4gas.env import get_llm_config
import os

def gen_google_adk_template():
    """
    Generate a Google ADK template for ils4gas.
    
    This function sets up the environment and initializes the agent with the necessary tools.
    """
    llm_model, llm_api_key, llm_base_url = get_llm_config()
    port = os.environ.get("ILS4GAS_PORT", "50001")
    transport = os.environ.get("ILS4GAS_TRANSPORT", "sse")
    if transport == "sse":
        mcp_toolset = "SseServerParams"
        url = f"http://localhost:{port}/sse"
    elif transport == "streamable-http":
        mcp_toolset = "StreamableHTTPServerParams"
        url = f"http://localhost:{port}/mcp"
    else:
        raise ValueError("Invalid transport protocol specified. Use 'sse' or 'streamable-http'.")

    template = f"""
from google.adk.agents import Agent
from google.adk.models.lite_llm import LiteLlm
from google.adk.tools.mcp_tool.mcp_session_manager import SseServerParams, StreamableHTTPServerParams
from dp.agent.adapter.adk import CalculationMCPToolset

toolset = CalculationMCPToolset(
    connection_params={mcp_toolset}(
        url="{url}",
    )
)

model = LiteLlm(
    model={repr(llm_model)},
    api_key={repr(llm_api_key)},
    base_url={repr(llm_base_url)}
)
root_agent = Agent(
    name='agent',
    model=model,
    description=(
        "ils4gas agent"
    ),
    instruction="You are an agent that can use tools to solve problems",
    tools=[toolset]
)
    """
    
    os.makedirs("myagent", exist_ok=True)
    with open("myagent/agent.py", "w") as f:
        f.write(template)
    with open("myagent/__init__.py", "w") as f:
        f.write("from . import agent")
        
    print("Google ADK template generated, please type below command to run the agent:")
    print("    adk web --host 0.0.0.0 --port 50002\n")