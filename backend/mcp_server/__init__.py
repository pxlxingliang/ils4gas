import os

port = int(os.environ.get("ILS4GAS_MCP_PORT", "50001"))
host = os.environ.get("ILS4GAS_MCP_HOST", "0.0.0.0")

from mcp.server.fastmcp import FastMCP

mcp = FastMCP("ILS4GAS", port=port, host=host)
