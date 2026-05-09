# Bug Reference

Common pitfalls and solutions encountered during ILS4GAS development.

---

## 1. `async def` without `await` → coroutine leak

**Symptom**: `'coroutine' object is not subscriptable`, tool results are `<coroutine object ...>` strings.

**Root cause**: A method declared `async def` but containing no `await` and called without `await`. Python returns a coroutine instead of executing.

```python
# BAD — returns coroutine, never executes
async def _execute_tool_sync(self, name, args) -> str:
    return tool_info.call(**json.loads(args))

result = self._execute_tool_sync("tool", '{}')  # coroutine!
display = result[:500]  # TypeError: 'coroutine' object is not subscriptable

# GOOD
def _execute_tool_sync(self, name, args) -> str:
    return tool_info.call(**json.loads(args))
```

**Rule**: Only use `async def` when the body actually contains `await`. Call async methods with `await`.

---

## 2. SSE format mismatch after protocol change

**Symptom**: Streaming broken, content only appears after refresh, tool cards missing.

**Root cause**: Backend changed SSE event format but frontend parser was not updated.

```
Before:  {"content": "hello"}
         {"tool_call_start": {"tool_name": "x", "args": "{}"}}

After:   {"type": "content_chunk", "text": "hello"}
         {"type": "tool_call_start", "tool_name": "x", "args": "{}"}
```

**Fix**: When changing event protocols, update every consumer. In this project, frontend SSE parser is the only consumer.

**Rule**: Define the event schema in one place (`core/events.py` `AgentEvent`) and ensure all producers/consumers match. Frontend must handle `parsed.type` to dispatch to correct handler.

---

## 3. API format mismatch → silent React crash

**Symptom**: Web page renders briefly then goes white. No visible error.

**Root cause**: Frontend TypeScript type (`Message.segments: Segment[]`) differs from API response (`{role, content}` — no `segments`). React renders `msg.segments.length` → `TypeError: Cannot read property 'length' of undefined` → uncaught error → white screen.

```typescript
// Frontend type
interface Message { id: string; role: string; segments: Segment[]; }

// API response — NO segments field
{ "id": "msg_1", "role": "user", "content": "hello" }

// Crash
msg.segments.length  // TypeError
```

**Fix**: Validate/transform API data at the boundary before storing in state.

```typescript
const msgs = (d.messages || []).map((m: any) => ({
  id: m.id, role: m.role,
  segments: [{ type: "text", content: m.content || "" }],
}));
```

**Rule**: Never trust API data to match frontend types. Transform at the fetch boundary. Add error boundaries in React to display errors instead of white screen.

---

## 4. OpenAI tool_calls format must be `{id, type, function: {name, arguments}}`

**Symptom**: First message with tools works. Second message returns HTTP 400: `missing messages.tool_calls.function`.

**Root cause**: Tool calls saved in simplified flat format, sent back to API in next turn.

```json
// BAD — missing id, type, function wrapper
[{ "tool_name": "search", "args": "{}" }]

// REQUIRED format
[{ "id": "call_1", "type": "function", "function": { "name": "search", "arguments": "{}" } }]
```

**Fix**: Save `tool_calls` exactly as received from OpenAI streaming response. Store complete `{id, type, function}` objects.

**Rule**: When persisting tool call data for replay, preserve exact OpenAI format. Never reshape or flatten.

---

## 5. Incomplete tool call sequences in history → API validation failure

**Symptom**: Similar to #4 — 400 error on multi-turn conversations.

**Root cause**: Assistant messages with `content: null` and `tool_calls` are saved to history, but the corresponding `role: "tool"` response messages may be missing or filtered. OpenAI requires: if assistant has `tool_calls`, there MUST be following `tool` messages with matching `tool_call_id`.

**Fix**: Two approaches:
- **A)** Save complete sequences including tool messages, then load them all for API.
- **B)** Filter history to exclude tool-only messages when building API request.

This project uses **B** (simpler, no benefit to replaying tool calls):

```python
def _build_history(session_id):
    for msg in messages:
        if msg["role"] == "tool":
            continue  # skip tool results
        if msg["role"] == "assistant" and msg.get("tool_calls") and not msg.get("content"):
            continue  # skip assistant messages that only contain tool calls
        ...
```

**Rule**: When sending message history to OpenAI, ensure no orphaned `tool_calls`. Either save complete tool sequences or filter them out.

---

## 6. Textual captures mouse → terminal text selection requires Shift

**Symptom**: Cannot select text with mouse drag in TUI. Right-click copies only the word under cursor.

**Root cause**: Textual enables terminal mouse mode (SGR tracking). All mouse events go to Textual for widget interaction (scrolling, clicking). The terminal cannot handle drag-based text selection unless Shift is held.

**Fix**: This is by-design for Textual. Options:
- Hold **Shift** while dragging to select text (terminal native selection)
- Show this tip prominently (status bar, `/help`)

**Rule**: Textual apps cannot support default mouse text selection. Always inform users about Shift+drag.

---

## 7. FastAPI service singletons: `assert` vs lazy init

**Symptom**: `GET /api/v1/sessions` returns 500 AssertionError → white screen on first load.

**Root cause**: `deps.py` used `assert service is not None` to access services initialized in `lifespan`. If lifespan hasn't run (TestClient, early request), assertion fails with no useful error message.

```python
# BAD — crashes with AssertionError
def get_session_service():
    assert session_service is not None
    return session_service

# GOOD — lazy initialization
def get_session_service():
    global session_service
    if session_service is None:
        session_service = SessionService()
    return session_service
```

**Rule**: Service singletons accessed by API routes should be lazy-initialized, not asserted. This handles cold starts, tests, and lifespan failures gracefully.

---

## 8. FastAPI SPA catch-all route ordering

**Symptom**: `/api/*` routes return HTML (index.html) instead of JSON.

**Root cause**: SPA catch-all `/{full_path:path}` registered before API routes, intercepting all requests.

**Fix**: Register API routers first, then SPA catch-all last. FastAPI matches routes in registration order.

```python
# 1. API routes first
app.include_router(chat.router)
app.include_router(llm.router)
...

# 2. SPA catch-all LAST
@app.get("/{full_path:path}")
async def spa(full_path: str):
    return index.html
```

**Rule**: Always register catch-all/wildcard routes AFTER all specific routes.

---

## 9. Duplicate tool execution from generator + sync call

**Symptom**: Tool called twice per LLM tool invocation.

**Root cause**: Tool executed inside an async generator for yielding events, then called again synchronously to get the return value for the message loop.

```python
# BAD — execute_tool called twice
async for event in self._execute_tool_stream(name, args):  # executes once
    yield event
result = self._execute_tool_sync(name, args)  # executes again

# GOOD — execute once, then emit event
result = self._execute_tool_sync(name, args)
yield self._make_tool_event(name, result)
```

**Rule**: When a function has side effects (tool execution), call it exactly once. Separate computation from event emission.

---

## 10. Session JSON field consistency

**Symptom**: Messages missing `tool_call_id` after persistence.

**Root cause**: `add_message()` didn't accept `tool_call_id` parameter. Tool response messages couldn't be saved with their call ID.

**Fix**: Add `tool_call_id` parameter to `add_message()` and include it in the stored message dict. Backward-compatible (defaults to `None`).

**Rule**: When adding new fields to persistence, always use optional parameters with `None` default to avoid breaking existing calls.

---

## Quick Checklist

- [ ] `async def` methods all contain `await`? Called with `await`?
- [ ] SSE/event format changed? Frontend parser updated?
- [ ] API format changed? Frontend type mapping updated?
- [ ] Tool calls stored in OpenAI `{id, type, function}` format?
- [ ] History filtered for incomplete tool call sequences?
- [ ] Textual app: Shift+drag tip visible to users?
- [ ] Service singletons: lazy init or lifespan guaranteed?
- [ ] SPA catch-all: registered AFTER API routes?
- [ ] Side effects: called exactly once per logical operation?
- [ ] New persistence fields: use `Optional` defaults?
