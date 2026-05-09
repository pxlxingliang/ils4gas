from typing import Any, Dict, Optional

from backend.skills import Skill
from backend.tools.registry import ToolRegistry


class SkillExecutor:
    def __init__(self, tool_registry: Optional[ToolRegistry] = None):
        self.tool_registry = tool_registry

    def execute(self, skill: Skill, context: Dict[str, Any]) -> Optional[str]:
        if skill.module and hasattr(skill.module, "execute"):
            return skill.module.execute(context)

        if skill.tools and self.tool_registry:
            for t in skill.tools:
                self.tool_registry.register(t)

        return None

    def inject_prompt(self, skill: Skill, system_prompt: str) -> str:
        skill_prompt = skill.get_system_prompt()
        return system_prompt + "\n\n" + skill_prompt
