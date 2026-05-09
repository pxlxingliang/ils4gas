# ILS4GAS Project Progress

Last updated: 2026-05-09 (Session 2)

## Overall Status: ~60% complete

Core chat, MCP tools, web UI, TUI, agent framework, skills, and workspace context are working.
Memory system, advanced agents, A2A, evaluation, built-in tools, and middleware are not yet implemented.

---

## Completed Modules

### core/ — Core Framework Layer (90%)

| File | Lines | Status | Notes |
|------|-------|--------|-------|
| `agent.py` | 49 | ✓ | `BaseAgent` abstract, `AgentState` state machine, cancel/pause/resume |
| `llm.py` | 88 | ✓ | `LLMProvider` + `OpenAICompatibleProvider` (sync/async/streaming) |
| `message.py` | 18 | ✓ | `Message` dataclass with `to_openai_format()` |
| `config.py` | 42 | ✓ | Loads `~/.ils4gas/config.json`, resolves `${ENV:VAR}` |
| `context.py` | 73 | ✓ | `WorkspaceContext` — AGENT/PERSONA/MEMORY.md → system prompt |
| `events.py` | 31 | ✓ | `AgentEventType` enum + `AgentEvent` dataclass (sse/json output) |
| `exceptions.py` | 10 | ✓ | `AgentError`, `ToolNotFoundError`, `LLMError` |
| `port.py` | 13 | ✓ | `find_free_port()` |
| `env.py` | 50 | ✓ | `EnvManager` — ILS4GAS_* env var defaults + persistence |

### agents/ — Agent Implementations (60%)

| File | Lines | Status | Notes |
|------|-------|--------|-------|
| `simple_agent.py` | 40 | ✓ | Pure LLM pass-through, no tool calling |
| `react_agent.py` | 205 | ✓ | ReAct pattern: Thought→Action loop, full tool calling, streaming events |
| `reflection_agent.py` | 12 | ○ Stub | `raise NotImplementedError` |
| `plan_solve_agent.py` | 12 | ○ Stub | `raise NotImplementedError` |
| `a2a_agent.py` | 0 | ✗ Empty | A2A orchestrator entry point |

### services/ — Business Service Layer (70%)

| File | Lines | Status | Notes |
|------|-------|--------|-------|
| `llm_service.py` | 95 | ✓ | Provider management, model switching, OpenAI tool schemas |
| `session_service.py` | 133 | ✓ | JSON file persistence (`~/.ils4gas/data/sessions/`) |
| `mcp_service.py` | 55 | ✓ | MCP server connect/disconnect/tool execution |
| `title_service.py` | 25 | ✓ | Auto-generate session titles via LLM on first message |
| `memory_service.py` | 0 | ✗ | |
| `skill_service.py` | 0 | ✗ | |
| `workspace_service.py` | 0 | ✗ | |
| `scheduler_service.py` | 0 | ✗ | |

### skills/ — Skills System (80%)

| File | Lines | Status | Notes |
|------|-------|--------|-------|
| `__init__.py` | 30 | ✓ | `Skill`, `SkillMeta` dataclasses |
| `registry.py` | 65 | ✓ | Discover, load, load_all, match_by_keyword, get_active_prompts |
| `loader.py` | 75 | ✓ | Parse SKILL.md YAML, dynamic skill.py/tools.py import |
| `executor.py` | 22 | ✓ | Execute Skill, inject prompt |
| `validator.py` | 28 | ✓ | Validate skill dir and SKILL.md |
| `creator.py` | 20 | ✓ | `create_manual()`, `create_from_conversation()` is stub |
| `builtin/code_review.py` | 0 | ✗ | |
| `builtin/data_analysis.py` | 0 | ✗ | |

### tools/ — Tool System (55%)

| File | Lines | Status | Notes |
|------|-------|--------|-------|
| `base.py` | 62 | ✓ | `Tool` ABC (name, description, parameters, call, execute, to_openai_schema) |
| `registry.py` | 44 | ✓ | `ToolInfo(Tool)` + `ToolRegistry` |
| `loader.py` | 48 | ✓ | `load_tools_from_mcp_modules()` — auto-discover MCP tool modules |
| `mcp_adapter.py` | 180+ | ✓ | MCP transport adapters (stdio/SSE/HTTP) |
| `mcp_manager.py` | 148 | ✓ | Multi-server lifecycle, lazy loading, tool discovery |
| `chain.py` | 0 | ✗ | Tool pipeline orchestration |
| `async_executor.py` | 0 | ✗ | Timeout/retry/concurrency |
| `permission.py` | 0 | ✗ | Tool permission grading |
| `builtin/` | 0 | ✗ | file_tools, web_tools, code_tools, system_tools |

### memory/ — Memory System (0%)

All 6 files empty: `short_term.py`, `long_term.py`, `vector_store.py`, `rag.py`, `summarizer.py`, `embeddings.py`

### a2a/ — Multi-Agent Communication (0%)

All 5 files empty: `protocol.py`, `router.py`, `orchestrator.py`, `client.py`, `agents/a2a_agent.py`

### evaluation/ — Evaluation System (0%)

All 4 files empty: `metrics.py`, `benchmarks.py`, `tracer.py`, `reporter.py`

### api/ — API Routes (90%)

| File | Status | Notes |
|------|--------|-------|
| `routes/chat.py` | ✓ | SSE streaming + sync send; ReActAgent + WorkspaceContext + Skills + auto-title |
| `routes/ws.py` | ✓ | WebSocket chat; shared Agent architecture |
| `routes/session.py` | ✓ | Full CRUD |
| `routes/llm.py` | ✓ | Model list, switch, current |
| `routes/mcp.py` | ✓ | MCP server/tool management |
| `routes/skill.py` | ✓ | Skill list endpoint |
| `routes/memory.py` | ✓ | Memory endpoints |
| `deps.py` | ✓ | Lazy-initialized service singletons (no crash on cold start) |
| `middleware/` | ✗ | All 4 files empty (auth, ratelimit, logging, input_filter) |

### mcp_server/ — Built-in MCP Server (100%)

| File | Status | Notes |
|------|--------|-------|
| `server.py` | ✓ | stdio/SSE/HTTP transport CLI |
| `modules/materials.py` | ✓ | 5 tools: generate_bulk_structure, generate_molecule_structure, search_properties, search_dft_feature, generate_novel_ions |
| `util/comm.py` | ✓ | Work path generator, command runner, core counter |

### tui/ — Text UI (100%)

| File | Status | Notes |
|------|--------|-------|
| `app.py` | ✓ | Textual App, CSS layout (history height: 1fr, bottom bar dock), Ctrl+Q quit |
| `cli.py` | ✓ | CLI entry |
| `screen.py` | ✓ | ChatScreen: ReActAgent + WorkspaceContext + Skills + auto-title, model switching, commands |
| `widgets/history.py` | ✓ | Scrollable MessageHistory with role-colored prefixes |
| `widgets/status_bar.py` | ✓ | Model/tools/help status |

### frontend/ — React SPA (95%)

| Component | Status | Notes |
|------|--------|-------|
| `App.tsx` | ✓ | Shell layout (sidebar + main) |
| `ChatPanel` | ✓ | Message list, SSE streaming (AgentEvent format), tool call cards, auto-title on first msg |
| `SessionList` | ✓ | Sidebar session CRUD |
| `ModelSelector` | ✓ | Model dropdown switcher |
| `ToolCallCard` | ✓ | Expandable tool call display |
| `MarkdownView` | ✓ | Markdown + syntax highlighting |
| `LocaleToggle` | ✓ | zh-CN / en-US |
| `store` | ✓ | Zustand-like global state |
| `i18n` | ✓ | zh-CN.json, en-US.json, useLocale |
| `theme` | ✓ | Light/dark CSS variables, ils4gas_theme localStorage |

---

## Bug Fixes (Session 2)

| Issue | Root Cause | Fix |
|-------|-----------|-----|
| TUI `'coroutine' object is not subscriptable` | `_execute_tool_sync` declared `async def` (no await) → returned coroutine | Changed to `def` |
| Tool call 400 error on next turn | tool_calls saved in wrong format (flat instead of OpenAI `{id, type, function}`); history loaded without filtering incomplete tool sequences | Fixed format + skip tool-only messages in history |
| TUI input box mid-screen | Vertical layout distributed space equally | CSS `#history { height: 1fr; }` |
| Web white screen | `deps.py` used `assert` → 500 crash; SPA only served `/` and `/chat` | Lazy init + catch-all route `/{full_path:path}` |
| Web content not streaming + tool cards missing | Frontend SSE parser used old format after AgentEvent refactor introduced new `{type, text/tool_name/args}` format | Updated to switch-case on `parsed.type` |
| All sessions titled "New Chat" | No title generation | `title_service.py` + auto-generate on first user message |
| React crash after history load | API returns `{role, content}` but frontend expects `{segments: [...]}` | Convert in `loadHistory()` |

---

## Architecture

All three interfaces share a single code path:

```
Web SSE / WebSocket / TUI
        │
        ▼
    ReActAgent (react_agent.py)
        ├── WorkspaceContext → ~/.ils4gas/workspace/{AGENT,PERSONA,MEMORY}.md
        ├── SkillRegistry    → ~/.ils4gas/skills/{name}/SKILL.md + skill.py
        └── ToolRegistry     → mcp_server/modules/*.py (MCP tools)
                │
                ▼
        LLMService → OpenAI-compatible API
```

### AgentEvent SSE Protocol

| Event Type | JSON Format |
|-----------|-------------|
| content_chunk | `{"type": "content_chunk", "text": "..."}` |
| reasoning_chunk | `{"type": "reasoning_chunk", "text": "..."}` |
| tool_call_start | `{"type": "tool_call_start", "tool_call_id": "...", "tool_name": "...", "args": "{...}"}` |
| tool_call_end | `{"type": "tool_call_end", "tool_name": "...", "result": "..."}` |
| done | `{"type": "done", "full_text": "..."}` |
| error | `{"type": "error", "message": "..."}` |

---

## User Directory Layout

```
~/.ils4gas/
├── config.json                  # LLM providers, models, MCP servers
├── env.json                     # Environment variable defaults
├── data/
│   ├── sessions.json            # Session index [{id, title, ...}]
│   └── sessions/
│       └── sess_*.json          # Per-session messages
├── workspace/
│   ├── AGENT.md                 # Agent behavior config
│   ├── PERSONA.md               # Agent persona
│   └── MEMORY.md                # Long-term memory
└── skills/
    └── code_review/
        ├── SKILL.md             # Skill metadata + prompt (YAML front matter)
        └── skill.py             # Skill execution logic
```

---

## Next Steps (Priority Order)

### P1 — Memory System
- `memory/short_term.py` — Token-aware sliding window
- `memory/long_term.py` — SQLite persistent memory
- Wire into ReActAgent for context-aware responses

### P2 — Built-in Tools
- `tools/builtin/file_tools.py` — File system operations
- `tools/builtin/web_tools.py` — HTTP requests, web scraping

### P3 — Advanced Agents
- `agents/reflection_agent.py` — Self-evaluation loop
- `agents/plan_solve_agent.py` — Plan-then-execute

### P4 — Long-term
- `a2a/` — Multi-agent orchestration
- `evaluation/` — Benchmarks and metrics
- `api/middleware/` — Auth, rate limiting, logging
