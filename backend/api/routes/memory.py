from fastapi import APIRouter

router = APIRouter(prefix="/api/v1/memory", tags=["memory"])


@router.get("/long-term")
async def list_memories():
    return {"memories": []}


@router.post("/search")
async def search_memories():
    return {"results": []}
