from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict


class AgentEventType(str, Enum):
    CONTENT_CHUNK = "content_chunk"
    REASONING_CHUNK = "reasoning_chunk"
    TOOL_CALL_START = "tool_call_start"
    TOOL_CALL_END = "tool_call_end"
    DONE = "done"
    ERROR = "error"
    CANCELLED = "cancelled"


@dataclass
class AgentEvent:
    type: AgentEventType
    data: Dict[str, Any] = field(default_factory=dict)

    def to_sse(self) -> str:
        import json

        payload = {"type": self.type.value}
        payload.update(self.data)
        return f"data: {json.dumps(payload, ensure_ascii=False)}\n\n"

    def to_dict(self) -> Dict[str, Any]:
        payload = {"type": self.type.value}
        payload.update(self.data)
        return payload
