from typing import Dict, List, Optional

from backend.tools.mcp_adapter import MCPServerConfig
from backend.tools.mcp_manager import MCPServerManager


class MCPService:
    def __init__(self, mcp_manager: MCPServerManager):
        self._manager = mcp_manager

    @classmethod
    def from_config_entries(
        cls, entries: List[Dict], manager: MCPServerManager | None = None
    ) -> "MCPService":
        mgr = manager or MCPServerManager()
        configs = [
            MCPServerConfig.from_dict(e) for e in entries if e.get("name")
        ]
        mgr.register_configs(configs)
        return cls(mgr)

    @property
    def manager(self) -> MCPServerManager:
        return self._manager

    async def connect_all(self) -> None:
        await self._manager.connect_all()

    async def connect_server(self, name: str) -> Optional[str]:
        cfg = self._manager._configs.get(name)
        if not cfg:
            return f"Server '{name}' not found in config"
        try:
            await self._manager.connect_server(cfg)
            return None
        except Exception as e:
            return str(e)

    async def disconnect_server(self, name: str) -> None:
        await self._manager.disconnect_server(name)

    async def disconnect_all(self) -> None:
        await self._manager.disconnect_all()

    def list_servers(self) -> List[Dict]:
        return self._manager.list_servers()

    def list_all_tools(self) -> List[Dict]:
        return self._manager.list_tools()

    def get_tool(self, name: str) -> Optional[Dict]:
        return self._manager.get_tool(name)

    async def execute_tool(self, name: str, arguments: Dict) -> str:
        return await self._manager.execute_tool(name, arguments)
