from fastapi import APIRouter

router = APIRouter(prefix="/api/v1/skills", tags=["skills"])


@router.get("")
async def list_skills():
    return {"skills": []}
