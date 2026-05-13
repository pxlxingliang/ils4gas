# Project Structure

Code lives in the project directory. Runtime data lives in user home.

## Root — Project Directory (version controlled)

| File | Purpose |
|------|---------|
| `pyproject.toml` | Python project metadata, build config, entry point (`ils4gas` command) |
| `requirements.txt` | Python backend dependencies |
| `config.example.json` | Config template — copied to `~/.ils4gas/config.json` on first launch |
| `Dockerfile` | Container image build |
| `docker-compose.yml` | Dev environment (backend + frontend + ChromaDB) |
| `README.md` | Project overview and quick start guide |
| `.gitignore` | Git ignore rules |

---

## `backend/` — Python Backend

### `backend/core/` — Core Framework Layer

Lowest-level abstractions with no external dependencies.

| File | Purpose |
|------|---------|
| `agent.py` | Agent abstract base class, event bus, interrupt/resume |
| `llm.py` | `LLMProvider` abstract + `OpenAICompatibleProvider` implementation |
| `message.py` | Message model (user / assistant / system / tool) |
| `config.py` | Config loader — reads `~/.ils4gas/config.json` + resolves `${ENV:VAR}` |
| `context.py` | Workspace context loader (AGENT.md, PERSONA.md, etc.) |
| `events.py` | `AgentEvent` enum + `EventBus` |
| `exceptions.py` | Custom exception types |

### `backend/agents/` — Agent Implementations

| File | Purpose |
|------|---------|
| `simple_agent.py` | Basic chat agent, no tool calling |
| `react_agent.py` | ReAct pattern: Thought → Action loop |
| `reflection_agent.py` | Self-reflection optimization agent |
| `a2a_agent.py` | A2A multi-agent orchestrator entry point |

### `backend/tools/` — Tool System

| File | Purpose |
|------|---------|
| `base.py` | `Tool` base class (sync/async unified interface) |
| `registry.py` | Tool registration, discovery, namespace isolation |
| `mcp_adapter.py` | MCP transport adapters (stdio / SSE / HTTP) |
| `mcp_manager.py` | Multi-server lifecycle management, lazy loading |
| `loader.py` | Load MCP module tools into ToolRegistry |
| `chain.py` | Tool pipeline orchestration |
| `async_executor.py` | Timeout / retry / concurrency control executor |
| `permission.py` | `ToolPermission` grading + permission manager |
| `builtin/read.py` | Read files with offset/limit |
| `builtin/glob.py` | File pattern matching |
| `builtin/grep.py` | Content search with regex |
| `builtin/bash.py` | Shell command execution |
| `builtin/write.py` | Create/overwrite files |
| `builtin/edit.py` | Exact string replacement |
| `builtin/webfetch.py` | Fetch web content |
| `builtin/websearch.py` | Web search (DuckDuckGo) |
| `builtin/skill.py` | Skill registration, discovery, on-demand loading (`load_skill` tool) |

### `backend/memory/` — Memory & Knowledge Layer

| File | Purpose |
|------|---------|
| `short_term.py` | Sliding window short-term memory (token-aware) |
| `long_term.py` | SQLite long-term memory persistence |
| `vector_store.py` | ChromaDB vector store wrapper |
| `rag.py` | RAG retrieval: semantic search + document indexing |
| `summarizer.py` | Conversation summarizer (LLM compresses long chats) |
| `embeddings.py` | Embedding model management (local / cloud modes) |

### `backend/a2a/` — A2A Multi-agent Communication

| File | Purpose |
|------|---------|
| `protocol.py` | A2A message format + task definitions |
| `router.py` | Message routing (from_agent → to_agent) |
| `orchestrator.py` | Task decomposition → capability matching → dispatch |
| `client.py` | Cross-process / cross-network A2A client |

### `backend/evaluation/` — Evaluation System

| File | Purpose |
|------|---------|
| `metrics.py` | Agent run metrics collection and calculation |
| `benchmarks.py` | BFCL / GAIA benchmark adapters |
| `tracer.py` | Execution tracing and logging |
| `reporter.py` | Evaluation report generation |

### `backend/services/` — Business Service Layer

| File | Purpose |
|------|---------|
| `llm_service.py` | LLM management: reads `~/.ils4gas/config.json`, creates providers |
| `session_service.py` | Session CRUD (SQLite, stores at `~/.ils4gas/data/sessions.db`) |
| `mcp_service.py` | MCP server connect / disconnect / tool discovery |
| `memory_service.py` | Memory read / write / search |
| `workspace_service.py` | Context file read / write (`~/.ils4gas/workspace/`) |
| `scheduler_service.py` | APScheduler cron jobs |

### `backend/api/` — FastAPI Route Layer

| File | Purpose |
|------|---------|
| `routes/chat.py` | Message send, streaming, cancel / pause / resume |
| `routes/llm.py` | Model list, switch, info |
| `routes/session.py` | Session CRUD |
| `routes/mcp.py` | MCP server and tool management |
| `routes/skill.py` | Skill management |
| `routes/memory.py` | Memory management |
| `routes/ws.py` | WebSocket streaming + SSE fallback |
| `middleware/auth.py` | Bearer token authentication |
| `middleware/ratelimit.py` | Sliding window rate limiting |
| `middleware/logging.py` | Request logging + Request ID tracing |
| `middleware/input_filter.py` | Input length / zero-width char filtering |
| `deps.py` | FastAPI dependency injection (service singletons) |

`backend/main.py` — FastAPI app creation, middleware registration, route mounting, lifespan

---

## `frontend/` — React Frontend

| File | Purpose |
|------|---------|
| `index.html` | SPA entry HTML |
| `package.json` | npm dependencies and scripts |
| `vite.config.ts` | Vite build config (proxy /api to backend) |
| `tsconfig.json` | TypeScript config |
| `src/main.tsx` | React mount entry |
| `src/App.tsx` | Root component (AppShell layout) |
| `src/store/index.ts` | Zustand global state |
| `src/locales/zh-CN.json` | Chinese locale dictionary |
| `src/locales/en-US.json` | English locale dictionary |
| `src/locales/index.ts` | Locale type definitions |
| `src/hooks/useLocale.ts` | i18n context + `useLocale()` + toggle |
| `src/hooks/useWebSocket.ts` | WebSocket connection + SSE fallback |
| `src/hooks/useStream.ts` | Streaming message state management |
| `src/hooks/useSessions.ts` | Session list state management |
| `src/styles/theme.css` | CSS variables theme (light / dark) |
| `src/styles/global.css` | Global base styles |
| `src/components/ChatPanel/` | Main chat area (message list + input) |
| `src/components/SessionList/` | Sidebar session list |
| `src/components/ModelSelector/` | Model dropdown switcher |
| `src/components/ToolCallCard/` | Collapsible tool call card (with diff view) |
| `src/components/MarkdownView/` | Markdown + syntax highlighting |
| `src/components/common/LocaleToggle.tsx` | Language toggle button (中 / EN) |

---

## `tests/` — Tests

| File | Scope |
|------|-------|
| `test_agent.py` | Agent base class, SimpleAgent |
| `test_react_agent.py` | ReAct think-act loop |
| `test_llm.py` | Provider abstraction, model switching |
| `test_tools.py` | Tool registration, executor, permissions |
| `test_mcp.py` | MCP adapters, server management |
| `test_skills.py` | Skill load, execute, create |
| `test_memory.py` | Memory read/write, RAG search |
| `test_evaluation.py` | Metrics, LLM-assisted evaluation |

---

## `docs/` — Documentation

| File | Purpose |
|------|---------|
| `DESIGN.md` | Technical design document |
| `STRUCTURE.md` | This file — directory and file descriptions |

---

## `~/.ils4gas/` — User Config (NOT version controlled)

| File | Purpose |
|------|---------|
| `config.json` | Single merged config: LLM providers, models, MCP servers, service settings |

First-launch behavior: copies `config.example.json` from the installed package to this location if `config.json` does not exist.

---

## `~/.ils4gas/` — User Data (NOT version controlled)

| Path | Purpose |
|------|---------|
| `.env` | API keys and secrets (never committed) |
| `data/sessions.db` | SQLite session + message persistence |
| `data/memory.db` | SQLite long-term memory persistence |
| `data/vectors/` | ChromaDB vector data |
| `workspace/AGENT.md` | Agent behavior configuration |
| `workspace/PERSONA.md` | Agent persona definition |
| `workspace/MEMORY.md` | Long-term memory file |
| `skills/` | User-installed skills |
| `logs/` | Runtime logs |

First-launch behavior: creates full directory tree and default files if they do not exist.
