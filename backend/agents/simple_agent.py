from typing import AsyncIterator, Dict, List

from backend.core.agent import BaseAgent, AgentState
from backend.core.events import AgentEvent, AgentEventType


class SimpleAgent(BaseAgent):
    async def run(self, messages: List[Dict]) -> str:
        self._state = AgentState.RUNNING
        try:
            result = await self.llm.ainvoke(messages)
            return result
        finally:
            self._state = AgentState.IDLE

    async def stream_run(self, messages: List[Dict]) -> AsyncIterator[AgentEvent]:
        self._state = AgentState.RUNNING
        try:
            accumulated = ""
            async for chunk in self.llm.astream_invoke(messages):
                if self._state == AgentState.PAUSED:
                    break
                accumulated += chunk
                yield AgentEvent(
                    type=AgentEventType.CONTENT_CHUNK,
                    data={"text": chunk},
                )
            yield AgentEvent(
                type=AgentEventType.DONE,
                data={"full_text": accumulated},
            )
        except Exception as e:
            self._state = AgentState.ERROR
            yield AgentEvent(
                type=AgentEventType.ERROR,
                data={"message": str(e)},
            )
        finally:
            if self._state != AgentState.ERROR:
                self._state = AgentState.IDLE
