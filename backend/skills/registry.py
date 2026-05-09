from pathlib import Path
from typing import Dict, List, Optional

from backend.skills import Skill, SkillMeta

SKILLS_DIR = Path("~/.ils4gas/skills").expanduser()


class SkillRegistry:
    def __init__(self, skills_dirs: Optional[List[str]] = None):
        self._skills: Dict[str, Skill] = {}
        self._dirs = [Path(d).expanduser() for d in (skills_dirs or [])]
        if not self._dirs:
            self._dirs = [SKILLS_DIR]

    def discover(self) -> List[str]:
        discovered: List[str] = []
        for base_dir in self._dirs:
            if not base_dir.exists():
                continue
            for skill_dir in base_dir.iterdir():
                if skill_dir.is_dir() and (skill_dir / "SKILL.md").exists():
                    discovered.append(skill_dir.name)
        return sorted(discovered)

    def load(self, name: str) -> Optional[Skill]:
        if name in self._skills:
            return self._skills[name]

        for base_dir in self._dirs:
            skill_dir = base_dir / name
            if not skill_dir.exists() or not (skill_dir / "SKILL.md").exists():
                continue

            try:
                skill = _load_skill_from_dir(skill_dir)
                self._skills[name] = skill
                return skill
            except Exception:
                return None
        return None

    def load_all(self) -> List[Skill]:
        discovered = self.discover()
        return [self.load(name) for name in discovered if self.load(name)]

    def get(self, name: str) -> Optional[Skill]:
        return self._skills.get(name)

    def match_by_keyword(self, user_input: str) -> List[str]:
        matched: List[str] = []
        for name, skill in self._skills.items():
            for kw in skill.meta.trigger_keywords:
                if kw.lower() in user_input.lower():
                    matched.append(name)
                    break
        return matched

    def get_active_prompts(self) -> str:
        parts: list[str] = []
        for skill in self._skills.values():
            parts.append(skill.get_system_prompt())
        return "\n".join(parts)


def _load_skill_from_dir(skill_dir: Path) -> Skill:
    from backend.skills.loader import SkillLoader

    loader = SkillLoader()
    meta, prompt = loader.parse_skill_md(skill_dir / "SKILL.md")
    module = loader.load_skill_py(skill_dir / "skill.py")
    tools = loader.load_tools(skill_dir / "tools.py", meta.name)

    return Skill(meta=meta, prompt=prompt, module=module, tools=tools)
