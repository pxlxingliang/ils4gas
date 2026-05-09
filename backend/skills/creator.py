from typing import List, Optional

from backend.skills import Skill, SkillMeta


class SkillCreator:
    @staticmethod
    async def create_from_conversation(
        conversation: List[dict], skill_name: str, llm=None
    ) -> Optional[Skill]:
        raise NotImplementedError("SkillCreator is not yet implemented")

    @staticmethod
    def create_manual(
        name: str,
        prompt: str,
        description: str = "",
        trigger_keywords: Optional[List[str]] = None,
    ) -> Skill:
        meta = SkillMeta(
            name=name,
            description=description,
            trigger_keywords=trigger_keywords or [],
        )
        return Skill(meta=meta, prompt=prompt)
