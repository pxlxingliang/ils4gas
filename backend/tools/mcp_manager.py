import asyncio
from typing import Dict, List, Optional, Any

from loguru import logger

from backend.tools.mcp_adapter import (
    MCPServerConfig,
    MCPTransport,
    create_transport,
)

logger = logger.bind(module="mcp_manager")


class MCPServerManager:
    def __init__(self):
        self._transports: Dict[str, MCPTransport] = {}
        self._tools: Dict[str, Dict] = {}
        self._tool_server_map: Dict[str, str] = {}
        self._configs: Dict[str, MCPServerConfig] = {}

    def register_configs(self, configs: List[MCPServerConfig]) -> None:
        for cfg in configs:
            self._configs[cfg.name] = cfg

    async def connect_all(self) -> None:
        for cfg in self._configs.values():
            if cfg.auto_connect:
                try:
                    await self.connect_server(cfg)
                except Exception as e:
                    logger.warning(f"Failed to connect MCP server '{cfg.name}': {e}")

    async def connect_server(self, config: MCPServerConfig) -> None:
        if config.name in self._transports:
            await self.disconnect_server(config.name)

        transport = create_transport(config)
        await transport.connect()
        self._transports[config.name] = transport
        await self._discover_tools(config.name, transport, config)

    async def disconnect_server(self, name: str) -> None:
        transport = self._transports.pop(name, None)
        if not transport:
            return
        self._remove_tools(name)
        await transport.disconnect()
        logger.info(f"MCP server '{name}' disconnected")

    async def disconnect_all(self) -> None:
        for name in list(self._transports.keys()):
            try:
                await self.disconnect_server(name)
            except Exception as e:
                logger.warning(f"Error disconnecting '{name}': {e}")

    async def _discover_tools(
        self, server_name: str, transport: MCPTransport, config: MCPServerConfig
    ) -> None:
        try:
            tools = await transport.list_tools()
        except Exception as e:
            logger.warning(f"Failed to list tools from '{server_name}': {e}")
            return

        for tool in tools:
            name = tool["name"]
            if config.tool_denylist and name in config.tool_denylist:
                continue
            if config.tool_allowlist and name not in config.tool_allowlist:
                continue

            full_name = f"{server_name}__{name}"
            self._tools[full_name] = {
                "name": full_name,
                "original_name": name,
                "server": server_name,
                "description": tool.get("description", ""),
                "inputSchema": tool.get("inputSchema", {}),
            }
            self._tool_server_map[full_name] = server_name
            logger.debug(f"  tool: {full_name}")

        logger.info(
            f"Discovered {len(tools)} tools from MCP server '{server_name}'"
        )

    def _remove_tools(self, server_name: str) -> None:
        to_remove = [
            name
            for name, server in self._tool_server_map.items()
            if server == server_name
        ]
        for name in to_remove:
            self._tools.pop(name, None)
            self._tool_server_map.pop(name, None)

    def list_servers(self) -> List[Dict]:
        return [
            {
                "name": cfg.name,
                "description": cfg.description,
                "transport": cfg.transport,
                "connected": cfg.name in self._transports,
                "auto_connect": cfg.auto_connect,
            }
            for cfg in self._configs.values()
        ]

    def list_tools(self, server_name: Optional[str] = None) -> List[Dict]:
        if server_name:
            return [
                t
                for t in self._tools.values()
                if t["server"] == server_name
            ]
        return list(self._tools.values())

    def get_tool(self, name: str) -> Optional[Dict]:
        if name in self._tools:
            return self._tools[name]
        for full_name, tool in self._tools.items():
            if full_name.endswith(f"__{name}"):
                return tool
        return None

    async def execute_tool(
        self, name: str, arguments: Dict
    ) -> str:
        tool = self.get_tool(name)
        if not tool:
            return f"Tool '{name}' not found"

        server_name = tool["server"]
        transport = self._transports.get(server_name)
        if not transport:
            return f"MCP server '{server_name}' is not connected"

        try:
            return await asyncio.wait_for(
                transport.call_tool(tool["original_name"], arguments),
                timeout=self._configs[server_name].timeout_ms / 1000,
            )
        except asyncio.TimeoutError:
            return f"Tool '{name}' execution timed out"
        except Exception as e:
            return f"Tool execution error: {e}"
