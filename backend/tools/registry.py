from typing import Callable, Dict, List, Optional

from backend.tools.base import Tool


class ToolInfo(Tool):
    def __init__(self, name: str, description: str, parameters: Dict, fn: Callable):
        self._name = name
        self._description = description
        self._parameters = parameters
        self._fn = fn

    @property
    def name(self) -> str:
        return self._name

    @property
    def description(self) -> str:
        return self._description

    @property
    def parameters(self) -> Dict:
        return self._parameters

    def execute(self, **kwargs):
        return self._fn(**kwargs)


class ToolRegistry:
    def __init__(self):
        self._tools: Dict[str, ToolInfo] = {}

    def register(self, tool: ToolInfo):
        self._tools[tool.name] = tool

    def get(self, name: str) -> Optional[ToolInfo]:
        return self._tools.get(name)

    def list_all(self) -> List[ToolInfo]:
        return list(self._tools.values())

    def to_openai_tools(self) -> List[Dict]:
        return [t.to_openai_schema() for t in self._tools.values()]

    def __len__(self) -> int:
        return len(self._tools)

    def __bool__(self) -> bool:
        return len(self._tools) > 0
