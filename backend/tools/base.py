import json
from abc import ABC, abstractmethod
from typing import Any, Dict


class Tool(ABC):
    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    @abstractmethod
    def description(self) -> str: ...

    @property
    @abstractmethod
    def parameters(self) -> Dict[str, Any]: ...

    def call(self, **kwargs) -> str:
        result = self.execute(**kwargs)
        if isinstance(result, (dict, list)):
            return json.dumps(result, ensure_ascii=False, default=str)
        return str(result)

    @abstractmethod
    def execute(self, **kwargs) -> Any:
        ...

    def to_openai_schema(self) -> Dict:
        return {
            "type": "function",
            "function": {
                "name": self.name,
                "description": self.description.strip(),
                "parameters": self._clean_schema(self.parameters),
            },
        }

    @staticmethod
    def _clean_schema(schema: Dict) -> Dict:
        out = {"type": schema.get("type", "object")}
        if "title" in schema:
            out["title"] = schema["title"]
        if "required" in schema:
            out["required"] = schema["required"]
        if "properties" in schema:
            props = {}
            for k, v in schema["properties"].items():
                cleaned = {"type": v.get("type", "string")}
                if "description" in v:
                    cleaned["description"] = v["description"]
                if "title" in v:
                    cleaned["title"] = v["title"]
                if "enum" in v:
                    cleaned["enum"] = v["enum"]
                if "items" in v:
                    cleaned["items"] = v["items"]
                if "default" in v:
                    cleaned["default"] = v["default"]
                props[k] = cleaned
            out["properties"] = props
        return out
