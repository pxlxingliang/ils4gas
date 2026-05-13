from fastapi import APIRouter

from backend.tools.builtin.skill import get_skill_registry

router = APIRouter(prefix="/api/v1/skills", tags=["skills"])


@router.get("")
async def list_skills():
    registry = get_skill_registry()
    return {"skills": registry.get_summaries()}
