from typing import AsyncIterator, Dict, List

from backend.core.agent import BaseAgent, AgentState
from backend.core.events import AgentEvent, AgentEventType
from backend.core.token_counter import count_tokens


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
        provider = self.llm.get_provider()
        model_name = provider.model_name

        system_messages = [{"role": "system", "content": self.system_prompt}]
        full_messages = system_messages + messages
        prompt_tokens = count_tokens(full_messages, model_name)

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
            done_data: dict = {
                "full_text": accumulated,
                "prompt_tokens": prompt_tokens,
            }
            usage = provider.last_usage
            if usage:
                done_data["usage"] = usage.to_dict()
            yield AgentEvent(
                type=AgentEventType.DONE,
                data=done_data,
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
