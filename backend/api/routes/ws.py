import json
from fastapi import APIRouter, WebSocket, WebSocketDisconnect

from backend.agents.react_agent import ReActAgent
from backend.api.deps import get_llm_service, get_session_service
from backend.core.context import WorkspaceContext
from backend.skills.registry import SkillRegistry

router = APIRouter()


def _make_agent():
    llm = get_llm_service()
    system_prompt = WorkspaceContext().build_system_prompt()

    registry = SkillRegistry()
    registry.load_all()
    skill_prompts = registry.get_active_prompts()
    if skill_prompts:
        system_prompt = system_prompt + "\n\n" + skill_prompts

    return ReActAgent(
        llm_service=llm,
        tool_registry=llm.tool_registry,
        system_prompt=system_prompt,
    )


@router.websocket("/api/v1/ws/chat")
async def ws_chat(websocket: WebSocket):
    await websocket.accept()
    sess = get_session_service()
    agent = _make_agent()

    try:
        while True:
            raw = await websocket.receive_text()
            data = json.loads(raw)

            msg_type = data.get("type", "")
            if msg_type == "cancel":
                agent.cancel()
                await websocket.send_json({"type": "done", "cancelled": True})
                continue

            content = data.get("content", "")
            session_id = data.get("session_id", "default")

            session = sess.get_session(session_id)
            if not session:
                session = sess.create_session()
                session_id = session["id"]

            messages = [
                {"role": msg["role"], "content": msg["content"]}
                for msg in sess.get_messages(session_id)
            ]
            messages.append({"role": "user", "content": content})
            sess.add_message(session_id, "user", content)

            text_accum: list[str] = []

            async for event in agent.stream_run(messages):
                await websocket.send_json(event.to_dict())
                if event.data.get("text") and event.type.value == "content_chunk":
                    text_accum.append(event.data["text"])

            full_text = "".join(text_accum)
            if full_text:
                sess.add_message(session_id, "assistant", full_text)

            await websocket.send_json({"type": "done", "session_id": session_id})

    except WebSocketDisconnect:
        pass
