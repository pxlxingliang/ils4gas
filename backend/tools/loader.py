import importlib
from pathlib import Path
from typing import Optional

from loguru import logger

from backend.tools.registry import ToolRegistry, ToolInfo

logger = logger.bind(module="tool_loader")


def load_tools_from_mcp_modules(modules_dir: Optional[Path] = None) -> ToolRegistry:
    registry = ToolRegistry()

    if modules_dir is None:
        modules_dir = Path(__file__).parent.parent / "mcp_server" / "modules"

    if not modules_dir.exists():
        logger.warning(f"Modules directory not found: {modules_dir}")
        return registry

    for py_file in sorted(modules_dir.glob("*.py")):
        if py_file.name.startswith("_"):
            continue

        module_name = f"backend.mcp_server.modules.{py_file.stem}"
        try:
            importlib.import_module(module_name)
            logger.debug(f"Imported module: {module_name}")
        except Exception as e:
            logger.warning(f"Failed to import {module_name}: {e}")

    try:
        from backend.mcp_server import mcp
        for tool in mcp._tool_manager.list_tools():
            info = ToolInfo(
                name=tool.name,
                description=tool.description or "",
                parameters=tool.parameters,
                fn=tool.fn,
            )
            registry.register(info)
            logger.info(f"Registered tool: {tool.name}")
    except Exception as e:
        logger.warning(f"Failed to extract tools from MCP: {e}")

    return registry
