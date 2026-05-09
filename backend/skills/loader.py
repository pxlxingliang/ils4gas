import importlib.util
import sys
from pathlib import Path
from typing import Any, List, Optional, Tuple

from backend.skills import SkillMeta
from backend.tools.base import Tool


class SkillLoader:
    @staticmethod
    def parse_skill_md(path: Path) -> Tuple[SkillMeta, str]:
        if not path.exists():
            raise FileNotFoundError(f"SKILL.md not found: {path}")

        content = path.read_text(encoding="utf-8")

        if content.startswith("---"):
            parts = content.split("---", 2)
            if len(parts) >= 3:
                import yaml
                front_matter = yaml.safe_load(parts[1]) or {}
                prompt = parts[2].strip()
            else:
                front_matter = {}
                prompt = content
        else:
            front_matter = {}
            prompt = content

        meta = SkillMeta(
            name=front_matter.get("name", path.parent.name),
            version=str(front_matter.get("version", "1.0")),
            description=front_matter.get("description", ""),
            author=front_matter.get("author", ""),
            tags=front_matter.get("tags", []),
            dependencies=front_matter.get("dependencies", {}),
            trigger_keywords=front_matter.get("trigger_keywords", []),
        )
        return meta, prompt

    @staticmethod
    def load_skill_py(path: Path) -> Optional[object]:
        if not path.exists():
            return None
        module_name = f"skill_{path.parent.name}"
        if module_name in sys.modules:
            return sys.modules[module_name]
        spec = importlib.util.spec_from_file_location(module_name, path)
        if spec is None or spec.loader is None:
            return None
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        return module

    @staticmethod
    def load_tools(path: Path, skill_name: str) -> List[Tool]:
        if not path.exists():
            return []
        module_name = f"skill_{skill_name}_tools"
        if module_name in sys.modules:
            return _extract_tools(sys.modules[module_name])
        spec = importlib.util.spec_from_file_location(module_name, path)
        if spec is None or spec.loader is None:
            return []
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        return _extract_tools(module)


def _extract_tools(module: Any) -> List[Tool]:
    tools: List[Tool] = []
    for attr_name in dir(module):
        obj = getattr(module, attr_name)
        if isinstance(obj, type) and issubclass(obj, Tool) and obj is not Tool:
            tools.append(obj())
        elif isinstance(obj, Tool):
            tools.append(obj)
    return tools
