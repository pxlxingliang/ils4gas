from dataclasses import dataclass, field
from typing import Dict, List, Optional

from backend.tools.base import Tool


@dataclass
class SkillMeta:
    name: str
    version: str = "1.0"
    description: str = ""
    author: str = ""
    tags: List[str] = field(default_factory=list)
    dependencies: Dict[str, List[str]] = field(default_factory=dict)
    trigger_keywords: List[str] = field(default_factory=list)


@dataclass
class Skill:
    meta: SkillMeta
    prompt: str = ""
    module: Optional[object] = None
    tools: List[Tool] = field(default_factory=list)

    def get_system_prompt(self) -> str:
        return (
            f'<skill name="{self.meta.name}">\n'
            f"{self.prompt}\n"
            f"</skill>"
        )
