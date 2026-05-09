from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any


@dataclass
class Message:
    role: str  # user | assistant | system | tool
    content: str
    tool_calls: Optional[List[Dict[str, Any]]] = None
    tool_call_id: Optional[str] = None

    def to_openai_format(self) -> Dict[str, Any]:
        msg: Dict[str, Any] = {"role": self.role, "content": self.content}
        if self.tool_calls:
            msg["tool_calls"] = self.tool_calls
        if self.tool_call_id:
            msg["tool_call_id"] = self.tool_call_id
        return msg
