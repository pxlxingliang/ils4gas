from pathlib import Path
from typing import List, Tuple


class SkillValidator:
    REQUIRED_FIELDS = ["name", "description"]

    @staticmethod
    def validate_skill_dir(skill_dir: Path) -> Tuple[bool, List[str]]:
        errors: List[str] = []

        if not skill_dir.exists():
            return False, [f"Directory not found: {skill_dir}"]

        skill_md = skill_dir / "SKILL.md"
        if not skill_md.exists():
            errors.append("Missing SKILL.md")

        if errors:
            return False, errors
        return True, []

    @staticmethod
    def validate_skill_md(path: Path) -> Tuple[bool, List[str]]:
        errors: List[str] = []
        if not path.exists():
            return False, ["SKILL.md not found"]

        content = path.read_text(encoding="utf-8")
        if not content.startswith("---"):
            errors.append("SKILL.md must start with YAML front matter (---)")

        if errors:
            return False, errors
        return True, []
