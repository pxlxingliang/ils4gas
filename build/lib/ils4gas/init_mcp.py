import os

port = os.environ.get("ILS4GAS_PORT", "50001")
host = os.environ.get("ILS4GAS_HOST", "0.0.0.0")
model = os.environ.get("ILS4GAS_MODEL", "fastmcp")
if model == "fastmcp":
    from mcp.server.fastmcp import FastMCP
    mcp = FastMCP("ILS4GAS", port=port, host=host)
elif model == "test": # For unit test of models
    class MCP:
        def tool(self):
            def decorator(func):
                return func
            return decorator
    mcp = MCP()
else:
    print("Please set the environment variable ILS4GAS_MODEL to fastmcp or test.")
    raise ValueError("Invalid ILS4GAS_MODEL. Please set it to fastmcp or test.")

