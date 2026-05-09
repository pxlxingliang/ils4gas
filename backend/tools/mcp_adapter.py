import asyncio
import json
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from loguru import logger

logger = logger.bind(module="mcp")


@dataclass
class MCPServerConfig:
    name: str
    description: str = ""
    transport: str = "stdio"
    command: Optional[str] = None
    args: List[str] = field(default_factory=list)
    url: Optional[str] = None
    env: Dict[str, str] = field(default_factory=dict)
    auto_connect: bool = True
    tool_allowlist: Optional[List[str]] = None
    tool_denylist: Optional[List[str]] = None
    timeout_ms: int = 30000

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MCPServerConfig":
        return cls(
            name=data["name"],
            description=data.get("description", ""),
            transport=data.get("transport", "stdio"),
            command=data.get("command"),
            args=data.get("args", []),
            url=data.get("url"),
            env=data.get("env", {}),
            auto_connect=data.get("auto_connect", True),
            tool_allowlist=data.get("tool_allowlist"),
            tool_denylist=data.get("tool_denylist"),
            timeout_ms=data.get("timeout_ms", 30000),
        )


class MCPTransport(ABC):
    def __init__(self, config: MCPServerConfig):
        self.config = config
        self._session: Any = None
        self._context: Any = None
        self._connected = False

    @abstractmethod
    async def connect(self) -> None:
        ...

    @abstractmethod
    async def disconnect(self) -> None:
        ...

    async def _init_session(self, read, write) -> None:
        from mcp import ClientSession

        self._session = ClientSession(read, write)
        await self._session.initialize()

    async def list_tools(self) -> List[Dict]:
        if not self._session:
            raise RuntimeError(f"MCP server '{self.config.name}' not connected")
        result = await self._session.list_tools()
        return [
            {
                "name": t.name,
                "description": t.description or "",
                "inputSchema": t.inputSchema if hasattr(t, "inputSchema") else {},
            }
            for t in result.tools
        ]

    async def call_tool(self, name: str, arguments: Dict) -> str:
        if not self._session:
            raise RuntimeError(f"MCP server '{self.config.name}' not connected")
        result = await self._session.call_tool(name, arguments)
        if result.isError:
            return f"Error: {result.content}"
        texts = []
        for c in result.content:
            if hasattr(c, "text"):
                texts.append(c.text)
            else:
                texts.append(str(c))
        return "\n".join(texts) if texts else "(no output)"

    async def list_resources(self) -> List[Dict]:
        if not self._session:
            return []
        result = await self._session.list_resources()
        return [
            {"uri": str(r.uri), "name": r.name, "description": r.description or ""}
            for r in result.resources
        ]

    @property
    def connected(self) -> bool:
        return self._connected


class StdioTransport(MCPTransport):
    async def connect(self) -> None:
        from mcp.client.stdio import stdio_client, StdioServerParameters

        env = {**os.environ, **self.config.env}
        params = StdioServerParameters(
            command=self.config.command,
            args=self.config.args,
            env=env,
        )
        self._context = stdio_client(params)
        read, write = await self._context.__aenter__()
        await self._init_session(read, write)
        self._connected = True
        logger.info(f"MCP stdio server '{self.config.name}' connected")

    async def disconnect(self) -> None:
        self._connected = False
        if self._context:
            try:
                await self._context.__aexit__(None, None, None)
            except Exception:
                pass
            self._context = None
        self._session = None


class SSETransport(MCPTransport):
    async def connect(self) -> None:
        from mcp.client.sse import sse_client

        self._context = sse_client(self.config.url)
        read, write = await self._context.__aenter__()
        await self._init_session(read, write)
        self._connected = True
        logger.info(f"MCP SSE server '{self.config.name}' connected at {self.config.url}")

    async def disconnect(self) -> None:
        self._connected = False
        if self._context:
            try:
                await self._context.__aexit__(None, None, None)
            except Exception:
                pass
            self._context = None
        self._session = None


class HTTPTransport(MCPTransport):
    async def connect(self) -> None:
        from mcp.client.streamable_http import streamablehttp_client

        self._context = streamablehttp_client(self.config.url)
        read, write = await self._context.__aenter__()
        await self._init_session(read, write)
        self._connected = True
        logger.info(f"MCP HTTP server '{self.config.name}' connected at {self.config.url}")

    async def disconnect(self) -> None:
        self._connected = False
        if self._context:
            try:
                await self._context.__aexit__(None, None, None)
            except Exception:
                pass
            self._context = None
        self._session = None


def create_transport(config: MCPServerConfig) -> MCPTransport:
    if config.transport == "stdio":
        return StdioTransport(config)
    elif config.transport == "sse":
        return SSETransport(config)
    elif config.transport == "http":
        return HTTPTransport(config)
    raise ValueError(f"Unsupported transport: {config.transport}")
