import asyncio
import json
from fastapi import APIRouter, WebSocket, WebSocketDisconnect

from backend.agents.react_agent import ReActAgent
from backend.api.deps import get_llm_service, get_session_service
from backend.core.context import WorkspaceContext

router = APIRouter()


def _make_agent():
    llm = get_llm_service()
    system_prompt = WorkspaceContext().build_system_prompt()

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
    _stream_task = None

    async def _run_stream(session_id: str, messages: list):
        text_accum: list[str] = []
        done_data: dict = {}
        try:
            async for event in agent.stream_run(messages):
                await websocket.send_json(event.to_dict())
                if event.data.get("text") and event.type.value == "content_chunk":
                    text_accum.append(event.data["text"])
                elif event.type.value == "done":
                    done_data = dict(event.data)
        except asyncio.CancelledError:
            pass
        except Exception:
            pass

        full_text = "".join(text_accum)
        if full_text and not any(tc.get("function", {}).get("name") for tc in (done_data.get("tool_calls") or [])):
            full_text += "\n\n*(Output truncated by user)*"

        if full_text:
            metadata = None
            usage = done_data.get("usage")
            if usage:
                metadata = {"token_usage": usage}
            sess.add_message(session_id, "assistant", full_text, metadata=metadata)

        final_done: dict = {"type": "done", "session_id": session_id}
        if "prompt_tokens" in done_data:
            final_done["prompt_tokens"] = done_data["prompt_tokens"]
        if "usage" in done_data:
            final_done["usage"] = done_data["usage"]
        try:
            await websocket.send_json(final_done)
        except Exception:
            pass

    try:
        while True:
            raw = await websocket.receive_text()
            data = json.loads(raw)

            msg_type = data.get("type", "")
            if msg_type == "cancel":
                agent.cancel()
                if _stream_task and not _stream_task.done():
                    _stream_task.cancel()
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

            _stream_task = asyncio.create_task(_run_stream(session_id, messages))
            await _stream_task

    except WebSocketDisconnect:
        if _stream_task and not _stream_task.done():
            _stream_task.cancel()
