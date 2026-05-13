import asyncio
import json
from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

from backend.agents.react_agent import ReActAgent
from backend.api.deps import get_llm_service, get_session_service
from backend.core.context import WorkspaceContext
from backend.services.title_service import generate_title

router = APIRouter(prefix="/api/v1/chat", tags=["chat"])


class ChatRequest(BaseModel):
    content: str


def _build_history(session_id: str) -> list:
    sess = get_session_service()
    messages = []
    for msg in sess.get_messages(session_id):
        role = msg["role"]
        if role == "tool":
            continue
        if role == "assistant" and msg.get("tool_calls") and not msg.get("content"):
            continue
        entry = {"role": role, "content": msg["content"]}
        if msg.get("tool_calls"):
            entry["tool_calls"] = msg["tool_calls"]
        if msg.get("tool_call_id"):
            entry["tool_call_id"] = msg["tool_call_id"]
        messages.append(entry)
    return messages


def _make_agent():
    llm = get_llm_service()
    system_prompt = WorkspaceContext().build_system_prompt()

    return ReActAgent(
        llm_service=llm,
        tool_registry=llm.tool_registry,
        system_prompt=system_prompt,
    )


async def _maybe_generate_title(session_id: str, user_content: str):
    sess = get_session_service()
    session = sess.get_session(session_id)
    if not session or session.get("title") != "New Chat":
        return
    try:
        llm = get_llm_service()
        title = await generate_title(user_content, llm)
        if title:
            sess.update_session(session_id, title=title)
    except Exception:
        pass


@router.post("/{session_id}/send")
async def send_message(session_id: str, req: ChatRequest):
    sess = get_session_service()
    if not sess.get_session(session_id):
        raise HTTPException(status_code=404, detail="Session not found")

    messages = _build_history(session_id)
    messages.append({"role": "user", "content": req.content})
    sess.add_message(session_id, "user", req.content)

    agent = _make_agent()
    response = await agent.run(messages)
    sess.add_message(session_id, "assistant", response)
    asyncio.create_task(_maybe_generate_title(session_id, req.content))
    return {"content": response, "session_id": session_id}


@router.post("/{session_id}/stream")
async def stream_message(session_id: str, req: ChatRequest):
    sess = get_session_service()
    if not sess.get_session(session_id):
        raise HTTPException(status_code=404, detail="Session not found")

    messages = _build_history(session_id)
    messages.append({"role": "user", "content": req.content})
    sess.add_message(session_id, "user", req.content)

    async def generate():
        try:
            agent = _make_agent()
            text_accum: list[str] = []
            tool_calls: list[dict] = []

            async for event in agent.stream_run(messages):
                yield event.to_sse()
                if event.type.value == "content_chunk" and event.data.get("text"):
                    text_accum.append(event.data["text"])
                elif event.type.value == "tool_call_start":
                    tool_calls.append({
                        "id": event.data.get("tool_call_id", ""),
                        "type": "function",
                        "function": {
                            "name": event.data.get("tool_name", ""),
                            "arguments": event.data.get("args", ""),
                        },
                    })
                elif event.type.value == "tool_call_end":
                    if tool_calls:
                        tool_calls[-1]["_result"] = event.data.get("result", "")

            full_text = "".join(text_accum)
            if full_text:
                saved_calls = [
                    {k: v for k, v in tc.items() if k != "_result"}
                    for tc in tool_calls
                ] if tool_calls else None
                sess.add_message(
                    session_id, "assistant", full_text,
                    tool_calls=saved_calls,
                )
                for tc in tool_calls:
                    if "_result" in tc:
                        sess.add_message(
                            session_id, "tool", tc["_result"],
                            tool_call_id=tc["id"],
                        )

            yield "data: [DONE]\n\n"
            asyncio.create_task(_maybe_generate_title(session_id, req.content))

        except Exception as e:
            yield f"data: {json.dumps({'type': 'error', 'message': str(e)}, ensure_ascii=False)}\n\n"
            yield "data: [DONE]\n\n"

    return StreamingResponse(
        generate(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",
        },
    )


@router.get("/{session_id}/history")
async def get_history(session_id: str):
    sess = get_session_service()
    if not sess.get_session(session_id):
        raise HTTPException(status_code=404, detail="Session not found")
    return {"messages": sess.get_messages(session_id)}
