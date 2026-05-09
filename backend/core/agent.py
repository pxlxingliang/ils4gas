from abc import ABC, abstractmethod
from enum import Enum
from typing import AsyncIterator, Dict, List, Optional

from backend.core.events import AgentEvent
from backend.core.exceptions import AgentError


class AgentState(str, Enum):
    IDLE = "idle"
    RUNNING = "running"
    PAUSED = "paused"
    ERROR = "error"


class BaseAgent(ABC):
    def __init__(
        self,
        llm_service,
        tool_registry=None,
        system_prompt: Optional[str] = None,
    ):
        self.llm = llm_service
        self.tools = tool_registry
        self.system_prompt = system_prompt or "You are a helpful assistant."
        self._state = AgentState.IDLE

    @property
    def state(self) -> AgentState:
        return self._state

    @abstractmethod
    async def run(self, messages: List[Dict]) -> str:
        ...

    @abstractmethod
    async def stream_run(self, messages: List[Dict]) -> AsyncIterator[AgentEvent]:
        ...

    def cancel(self) -> None:
        self._state = AgentState.IDLE

    def pause(self) -> None:
        if self._state == AgentState.RUNNING:
            self._state = AgentState.PAUSED

    def resume(self) -> None:
        if self._state == AgentState.PAUSED:
            self._state = AgentState.RUNNING
