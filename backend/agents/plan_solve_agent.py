from typing import AsyncIterator, Dict, List

from backend.core.agent import BaseAgent
from backend.core.events import AgentEvent


class PlanSolveAgent(BaseAgent):
    async def run(self, messages: List[Dict]) -> str:
        raise NotImplementedError("PlanSolveAgent is not yet implemented")

    async def stream_run(self, messages: List[Dict]) -> AsyncIterator[AgentEvent]:
        raise NotImplementedError("PlanSolveAgent is not yet implemented")
