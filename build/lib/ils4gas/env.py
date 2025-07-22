import os
import json
import time

ENVS = {
    "ILS4GAS_WORK_PATH": "/tmp/ils4gas",

    # connection settings
    "ILS4GAS_TRANSPORT": "sse",  # sse, streamable-http
    "ILS4GAS_HOST": "localhost",
    "ILS4GAS_PORT": "50001", 
    
    "LLM_MODEL": "",
    "LLM_API_KEY": "",
    "LLM_BASE_URL": "",
    
    "_comments":{
        "ILS4GAS_WORK_PATH": "The working directory, where all temporary files will be stored.",
        "ILS4GAS_TRANSPORT": "The transport protocol, can be 'sse' or 'streamable-http'.",
        "ILS4GAS_HOST": "The host address for the server.",
        "ILS4GAS_PORT": "The port number for the server.",
        "LLM_MODEL": "The model to use for the LLM, e.g., 'gpt-3.5-turbo'.",
        "LLM_API_KEY": "The API key for the LLM service.",
        "LLM_BASE_URL": "The base URL for the LLM service, if applicable."
    }
}

def set_envs(transport_input=None, port_input=None, host_input=None):
    """
    Set environment variables for AbacusAgent.
    
    Args:
        transport_input (str, optional): The transport protocol to use. Defaults to None.
        port_input (int, optional): The port number to run the MCP server on. Defaults to None.
        host_input (str, optional): The host address to run the MCP server on. Defaults to None.
    
    Returns:
        dict: The environment variables that have been set.
    """
    envjson_file = os.path.expanduser("~/.ils4gas/env.json")
    if os.path.isfile(envjson_file):
        envjson = json.load(open(envjson_file, "r"))
    else:
        envjson = {}
        
    update_envjson = False    
    for key, value in ENVS.items():
        if key not in envjson:
            envjson[key] = value
            update_envjson = True
    
    if transport_input is not None:
        envjson["ILS4GAS_TRANSPORT"] = str(transport_input)
    if port_input is not None:
        envjson["ILS4GAS_PORT"] = str(port_input)
    if host_input is not None:
        envjson["ILS4GAS_HOST"] = str(host_input)
        
    for key, value in envjson.items():
        os.environ[key] = str(value)
    
    if update_envjson:
        # write envjson to ~/.abacusagent/env.json
        os.makedirs(os.path.dirname(envjson_file), exist_ok=True)
        json.dump(
            envjson,
            open(envjson_file, "w"),
            indent=4
        )
    return envjson
    
def create_workpath():
    """
    Create the working directory for AbacusAgent, and change the current working directory to it.
    
    Returns:
        str: The path to the working directory.
    """
    work_path = os.environ.get("ILS4GAS_WORK_PATH", "/tmp/ils4gas") + f"/{time.strftime('%Y%m%d%H%M%S')}"
    os.makedirs(work_path, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(work_path)
    print(f"Changed working directory to: {work_path}")
    # write the environment variables to a file
    json.dump({
        k: os.environ.get(k) for k in ENVS.keys()
    }.update({"ILS4GAS_START_PATH": cwd}), 
        open("env.json", "w"), indent=4)
    
    return work_path    

def get_llm_config():
    """
    Get the LLM configuration from environment variables.
    
    Returns:
        tuple: Contains LLM model, API key, and base URL.
    """
    llm_model = os.environ.get("LLM_MODEL", "")
    llm_api_key = os.environ.get("LLM_API_KEY", "")
    llm_base_url = os.environ.get("LLM_BASE_URL", "")
    
    
    if not llm_model:
        llm_model = input("Enter LLM model (e.g., 'gpt-3.5-turbo'): ").strip()
        print("Note: you can also set the LLM_MODEL in ~/.ils4gas/env.json to avoid this prompt in the future.")
    if not llm_api_key:
        llm_api_key = input("Enter LLM API key: ").strip()
        print("Note: you can also set the LLM_API_KEY in ~/.ils4gas/env.json to avoid this prompt in the future.")
    if not llm_base_url:
        llm_base_url = input("Enter LLM base URL (if applicable): ").strip()
        print("Note: you can also set the LLM_BASE_URL in ~/.ils4gas/env.json to avoid this prompt in the future.")
    
    return llm_model, llm_api_key, llm_base_url


def load_tools():
    """
    Load all tools from the ils4gas package.
    """
    from pathlib import Path
    import importlib
    import inspect
    
    module_dir = Path(__file__).parent / "modules"
    
    function_names = []
    
    for py_file in module_dir.glob("*.py"):
        if py_file.name.startswith("_") or py_file.stem in ["utils", "comm"]: 
            continue  # skipt __init__.py and utils.py
        
        module_name = f"ils4gas.modules.{py_file.stem}"
        try:
            module = importlib.import_module(module_name)
            print(f"✅ Successfully loaded: {module_name}")

            for name, obj in inspect.getmembers(module_name, inspect.isfunction):
                function_names.append(name)

        except Exception as e:
            print(f"⚠️ Failed to load {module_name}: {str(e)}")
    return function_names

