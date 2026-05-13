# ILS4GAS 智能体框架 - 技术设计文档

## 1. 项目概述

### 1.1 项目背景

本项目旨在开发一个支持以下核心功能的智能体工具：
- **MCP工具调用**：集成已有的MCP工具函数，在合适的时机自动调用
- **多LLM支持**：支持多种大语言模型，并能在前端动态切换
- **会话管理**：支持多会话管理，类似opencode/openclaw的交互体验
- **远程访问**：支持在服务器部署，通过远程网页进行交互
- **美观界面**：提供简洁易用的Web前端界面

### 1.2 设计理念

本框架借鉴了 **hello-agents** 项目的核心设计思想：

1. **轻量级与教学友好**：代码简洁，易于理解和扩展
2. **统一工具抽象**："一切皆为工具"（包括MCP、记忆、RAG等）
3. **基于标准API**：兼容OpenAI接口，参考opencode的配置方式
4. **渐进式扩展**：模块化设计，按需添加功能
5. **完全控制**：开发者对每一行代码都有完全掌控

### 1.3 项目约定

1. **文件命名**：项目中所有文件和目录统一使用英文命名（如 `agent.py`、`session_service.py`）
2. **代码语言**：`src/` 下所有功能代码（含注释）使用英文编写
3. **前端国际化**：前端界面同时支持中文（默认）和英文，通过一个按钮切换。参照 6.6 节 i18n 方案设计
4. **文档**：项目核心文档使用英文命名（如 `DESIGN.md`、`README.md`）

### 1.4 数据目录模型

项目遵循 XDG 规范，将**代码**与**运行时数据**分离：

```
项目目录 (源代码)                用户目录 (运行时数据)
/your/project/                  ~/
├── backend/         安装到      ├── .ils4gas/
├── frontend/        site-       │   └── config.json        # 合并：模型+MCP+服务器配置
├── docs/            packages    └── .ils4gas/
├── tests/                           ├── data/
├── pyproject.toml                   │   ├── sessions.db     # 会话+消息
└── README.md                        │   ├── memory.db      # 长期记忆
                                     │   └── vectors/       # ChromaDB 向量
                                     ├── workspace/
                                     │   ├── AGENT.md       # Agent 行为配置
                                     │   ├── PERSONA.md     # Agent 人设
                                     │   └── MEMORY.md     # 长期记忆文件
                                     ├── skills/           # 用户安装的 Skill
                                     ├── logs/             # 运行日志
                                     └── .env              # API Key 等秘钥
```

**规则：**
- 项目目录仅含**源代码、构建配置、文档、测试**，不提交任何运行时配置或用户数据
- 首次启动时自动在 `~/.ils4gas/` 生成默认配置文件
- 首次启动时自动在 `~/.ils4gas/` 创建完整数据目录树
- `.env` 和 `config.json` 的 API Key 仅存在于用户目录，永不进入版本控制

---

## 2. 技术架构设计

### 2.1 整体架构图

```
┌───────────────────────────────────────────────────────────────────┐
│                          前端界面 (SPA)                            │
│  ┌──────────┐ ┌──────────┐ ┌──────────────┐ ┌──────────────────┐ │
│  │ 会话管理  │ │ 模型切换  │ │  聊天界面     │ │ 工具/Skill可视化 │ │
│  │ (侧边栏)  │ │ (顶栏)   │ │ (流式渲染)    │ │ (折叠卡片+diff)  │ │
│  └──────────┘ └──────────┘ └──────────────┘ └──────────────────┘ │
└─────────────────────┬─────────────────────────────────────────────┘
                      │ HTTP REST + WebSocket (SSE fallback)
                      ▼
┌───────────────────────────────────────────────────────────────────┐
│                     FastAPI 服务层                                 │
│  Middleware: CORS | Auth | RateLimit | RequestID | Logging        │
│  /api/v1/chat/* | /api/v1/sessions/* | /api/v1/llm/*             │
│  /api/v1/skills/* | /api/v1/mcp/* | /api/v1/memory/* | /ws/chat  │
└─────────────────────┬─────────────────────────────────────────────┘
                      │
                      ▼
┌───────────────────────────────────────────────────────────────────┐
│                        业务服务层                                  │
│  LLMService | SessionService | MCPService                          │
│  MemoryService | A2AService | WorkspaceService                    │
└──────┬──────────┬──────────┬──────────┬───────────────────────────┘
       │          │          │          │
       ▼          ▼          ▼          ▼
┌───────────┐ ┌──────────┐ ┌───────────┐ ┌──────────────────┐
│ 核心框架层 │ │ 工具系统层 │ │ Skill系统  │ │  记忆与知识层     │
│ Agent基类  │ │ Tool基类  │ │ Skill注册  │ │ 短期记忆(会话)    │
│ LLM接口   │ │ MCP适配器  │ │ Skill加载  │ │ 长期记忆(向量DB)  │
│ Message   │ │ 工具注册表 │ │ Skill执行  │ │ RAG检索增强      │
│ 上下文工程 │ │ 异步执行器 │ │ Skill自动创建│ │ 知识图谱(可选)   │
└─────┬─────┘ └────┬─────┘ └─────┬─────┘ └────────┬─────────┘
      │            │             │                │
      └────────────┴──────┬──────┴────────────────┘
                          │
                          ▼
┌───────────────────────────────────────────────────────────────────┐
│                    数据存储层（用户目录）                            │
│  ~/.ils4gas/data/sessions.db | vectors/ | memory.db              │
│  ~/.ils4gas/config.json | skills/ | workspace/            │
└───────────────────────────────────────────────────────────────────┘

                    ┌──────────────────────┐
                    │   A2A 通信层(可选)    │
                    │  Agent间消息传递      │
                    │  任务分发与负载均衡    │
                    └──────────────────────┘
```

### 2.2 技术栈选型

| 层级 | 技术选型 | 说明 |
|------|---------|------|
| **后端框架** | FastAPI | 异步高性能，自动生成API文档 |
| **智能体核心** | 自研（基于hello-agents理念） | 轻量级，完全可控 |
| **LLM接口** | OpenAI Python SDK + LiteLLM(可选) | 兼容所有OpenAI接口，可扩展100+提供商 |
| **MCP集成** | mcp (官方Python SDK) | 支持 stdio/SSE/HTTP 三种传输协议，支持工具懒加载 |
| **数据存储** | SQLite + ChromaDB (用户目录) | 运行时数据在 `~/.ils4gas/data/`，配置在 `~/.ils4gas/` |
| **前端界面** | React/Vue + Vite | 组件化开发，支持虚拟滚动、流式渲染 |
| **实时通信** | WebSocket (主) + SSE (fallback) | 流式响应，自动降级 |
| **配置管理** | Pydantic + Pydantic-Settings | 数据验证、环境变量解析、配置热加载 |
| **任务调度** | APScheduler | 定时任务、Cron 表达式支持 |
| **向量嵌入** | sentence-transformers / OpenAI Embeddings | 文本向量化，支持本地和云端两种模式 |
| **日志系统** | loguru | 结构化日志，支持请求追踪 |

---

## 3. 项目结构

```
ils4gas/
├── backend/
│   ├── core/                    # 核心框架层
│   │   ├── __init__.py
│   │   ├── agent.py
│   │   ├── llm.py
│   │   ├── message.py
│   │   ├── config.py
│   │   ├── context.py
│   │   ├── exceptions.py
│   │   └── events.py
│   │
│   ├── agents/                  # Agent实现层
│   │   ├── __init__.py
│   │   ├── simple_agent.py
│   │   ├── react_agent.py
│   │   ├── reflection_agent.py
│   │   └── a2a_agent.py
│   │
│   ├── tools/                   # 工具系统层
│   │   ├── __init__.py
│   │   ├── base.py
│   │   ├── registry.py
│   │   ├── mcp_adapter.py
│   │   ├── mcp_manager.py
│   │   ├── chain.py
│   │   ├── async_executor.py
│   │   ├── permission.py
│   │   ├── loader.py             # MCP 模块工具加载
│   │   └── builtin/              # 内置工具（含 Skill）
│   │       ├── __init__.py
│   │       ├── read.py
│   │       ├── glob.py
│   │       ├── grep.py
│   │       ├── bash.py
│   │       ├── write.py
│   │       ├── edit.py
│   │       ├── webfetch.py
│   │       ├── websearch.py
│   │       └── skill.py          # Skill 加载工具
│   │
│   ├── memory/                  # 记忆与知识层
│   │   ├── __init__.py
│   │   ├── short_term.py
│   │   ├── long_term.py
│   │   ├── vector_store.py
│   │   ├── rag.py
│   │   ├── summarizer.py
│   │   └── embeddings.py
│   │
│   ├── a2a/                     # A2A 多智能体通信层
│   │   ├── __init__.py
│   │   ├── protocol.py
│   │   ├── router.py
│   │   ├── orchestrator.py
│   │   └── client.py
│   │
│   ├── evaluation/              # Agent评估层
│   │   ├── __init__.py
│   │   ├── metrics.py
│   │   ├── benchmarks.py
│   │   ├── tracer.py
│   │   └── reporter.py
│   │
│   ├── services/               # 业务服务层
│   │   ├── __init__.py
│   │   ├── llm_service.py
│   │   ├── session_service.py
│   │   ├── mcp_service.py
│   │   ├── memory_service.py
│   │   ├── workspace_service.py
│   │   └── scheduler_service.py
│   │
│   ├── api/                    # API层
│   │   ├── __init__.py
│   │   ├── routes/
│   │   │   ├── chat.py
│   │   │   ├── llm.py
│   │   │   ├── session.py
│   │   │   ├── mcp.py
│   │   │   ├── skill.py
│   │   │   ├── memory.py
│   │   │   └── ws.py
│   │   ├── middleware/
│   │   │   ├── auth.py
│   │   │   ├── ratelimit.py
│   │   │   └── logging.py
│   │   └── deps.py
│   │
│   └── main.py
│
├── config.example.json          # 配置文件模板（首次启动复制到 ~/.ils4gas/）
│
├── frontend/                   # 前端界面
│   ├── src/
│   │   ├── components/
│   │   │   ├── ChatPanel/
│   │   │   ├── SessionList/
│   │   │   ├── ModelSelector/
│   │   │   ├── ToolCallCard/
│   │   │   ├── MarkdownView/
│   │   │   └── common/
│   │   ├── hooks/
│   │   │   ├── useWebSocket.ts
│   │   │   ├── useStream.ts
│   │   │   ├── useSessions.ts
│   │   │   └── useLocale.ts
│   │   ├── store/
│   │   │   └── index.ts
│   │   ├── locales/
│   │   │   ├── zh-CN.json
│   │   │   ├── en-US.json
│   │   │   └── index.ts
│   │   ├── styles/
│   │   │   ├── theme.css
│   │   │   └── global.css
│   │   ├── App.tsx
│   │   └── main.tsx
│   ├── index.html
│   ├── package.json
│   ├── vite.config.ts
│   └── tsconfig.json
│
├── tests/                     # 测试代码
│   ├── test_agent.py
│   ├── test_react_agent.py
│   ├── test_llm.py
│   ├── test_tools.py
│   ├── test_mcp.py
│   ├── test_skills.py
│   ├── test_memory.py
│   └── test_evaluation.py
│
├── docs/                      # 文档
│   ├── DESIGN.md
│   └── STRUCTURE.md
│
├── pyproject.toml
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
├── README.md
└── .gitignore
```

**用户目录运行时数据（不纳入版本控制）：**

```
~/
├── .config/ils4gas/
│   └── config.json            # 合并：模型+MCP+服务器配置
└── .ils4gas/
    ├── data/
    │   ├── sessions.db        # 会话+消息 SQLite
    │   ├── memory.db         # 长期记忆 SQLite
    │   └── vectors/          # ChromaDB 向量数据
    ├── workspace/
    │   ├── AGENT.md
    │   ├── PERSONA.md
    │   └── MEMORY.md
    ├── skills/               # 用户安装的 Skill
    ├── logs/                 # 运行日志
    └── .env                  # API Key 等秘钥
```

---

## 4. 核心设计

### 4.1 上下文工程（Workspace 配置层）

参考 hello-agents 上下文工程理念，Agent 在启动时会从工作区加载多层配置，构建完整的上下文理解。

**工作区文件结构（用户目录）：**
```
~/.ils4gas/workspace/
├── AGENT.md        # Agent 全局行为配置（工具使用策略、回复风格、安全规则）
├── PERSONA.md      # Agent 人设定义（角色、语气、专业领域）
├── MEMORY.md       # 长期记忆文件（用户偏好、历史决策、关键信息）
├── TOOLS.md        # 工具使用说明（可选，覆盖默认工具描述）
└── SKILLS.md       # Skill 使用指南（可选）
```

**上下文加载流程：** (详见 `backend/core/context.py`)

```
WorkspaceContext(workspace_path)
  ├── _files: 映射到 AGENT.md / PERSONA.md / MEMORY.md / TOOLS.md / SKILLS.md
  │
  ├── load_context() → {key: content}   // 读取所有存在的文件
  ├── build_system_prompt() → str       // 合并 persona + agent + memory + tools 为系统提示词
  ├── update_memory(content)            // 追加行到 MEMORY.md
  └── write_file(name, content)         // 写入指定的上下文文件
```

**上下文加载时机：**
- Agent 初始化时加载一次作为 system prompt
- 用户可通过对话指令动态修改（如 "请记住我喜欢简洁的回答" → 写入 MEMORY.md）
- 支持多用户隔离（每个用户有独立的 workspace 子目录）

### 4.2 LLM配置设计（多Provider抽象）

**Provider 抽象层设计：** (详见 `backend/core/llm.py`)

```
LLMProvider (ABC)
  ├── invoke(messages) → str           // 同步调用
  ├── stream_invoke(messages) → stream  // 同步流式调用
  ├── ainvoke(messages) → str          // 异步调用
  ├── astream_invoke(messages) → stream // 异步流式调用
  ├── model_name: str                  // 模型名称
  ├── provider_name: str               // 厂商名称
  └── context_limit: int               // 上下文长度限制

OpenAICompatibleProvider(LLMProvider)
  ├── __init__(api_key, base_url, model_name, context_limit)
  │     // 创建 OpenAI + AsyncOpenAI client
  ├── invoke / stream_invoke          // 通过 OpenAI client 调用
  └── ainvoke / astream_invoke        // 通过 AsyncOpenAI client 调用
```

**~/.ils4gas/config.json 结构：** (详见 `backend/core/config.py` 的 `load_config()`)

```
{
  currentModel: "provider/model_id",     // 当前使用的模型
  providers: {
    "<provider_key>": {
      name, type,                        // 厂商名称、类型(openai-compatible等)
      options: { baseURL, apiKey },      // 连接参数，apiKey 支持 ${ENV:VAR}
      models: {
        "<model_key>": {
          name,                           // 显示名称
          limit: { context, output },     // token 限制
          supports_tools,                 // 是否支持 function calling
          supports_reasoning              // 是否支持推理
        }
      }
    }
  }
}
```

**环境变量解析机制：**
- `${ENV:VOLC_API_KEY}` 表示从环境变量 `VOLC_API_KEY` 读取值
- 在加载配置时自动解析并替换

### 4.3 会话管理设计（SQLite持久化）

**会话表结构（~/.ils4gas/data/sessions.db）：**

```sql
-- 会话表
CREATE TABLE sessions (
    id TEXT PRIMARY KEY,
    title TEXT NOT NULL DEFAULT '新对话',
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    model_provider TEXT NOT NULL,
    model_name TEXT NOT NULL,
    metadata TEXT DEFAULT '{}',  -- JSON: token_count, mcp_servers_used, total_tool_calls
    is_active INTEGER DEFAULT 1
);

-- 消息表
CREATE TABLE messages (
    id TEXT PRIMARY KEY,
    session_id TEXT NOT NULL,
    role TEXT NOT NULL,  -- user / assistant / system / tool
    content TEXT NOT NULL,
    tool_calls TEXT,     -- JSON array: [{tool_name, args, result}]
    timestamp TEXT NOT NULL,
    FOREIGN KEY (session_id) REFERENCES sessions(id) ON DELETE CASCADE
);

-- 索引
CREATE INDEX idx_messages_session ON messages(session_id, timestamp);
CREATE INDEX idx_sessions_updated ON sessions(updated_at DESC);
```

**会话服务设计：** (详见 `backend/services/session_service.py`)

```
SessionService(db_path)
  ├── _init_db() → 创建 sessions + messages 表及索引
  ├── create_session(provider, model, title) → session_id
  ├── list_sessions(limit) → [{id, title, ...}]
  ├── add_message(session_id, role, content, tool_calls) → msg_id
  ├── get_messages(session_id, limit) → [{role, content, tool_calls, ...}]
  └── switch_model(session_id, provider, model) → 更新会话的模型配置
```

### 4.4 MCP工具集成方案（增强版）

MCP (Model Context Protocol) 是本框架的核心工具集成协议。支持三种传输方式，多服务器管理，以及工具懒加载。

**MCP服务器配置结构：** (详见 `~/.ils4gas/config.json` 的 `mcp_servers` 字段，解析逻辑见 `backend/tools/mcp_adapter.py` 的 `MCPServerConfig`)

```
mcp_servers: [
  {
    name, description, transport,     // 名称、描述、传输协议(stdio|sse|http)
    command, args, env,               // stdio 模式: 启动命令 + 环境变量
    url,                               // sse/http 模式: 服务地址
    auto_connect,                      // 是否自动连接
    tool_allowlist, tool_denylist,     // 工具白名单/黑名单
    timeout_ms                         // 超时
  }
]
```

**MCP传输层抽象：** (详见 `backend/tools/mcp_adapter.py`)

```
MCPTransport (ABC)
  ├── connect() / disconnect()        // 连接/断开
  ├── list_tools() → [schema]         // 发现工具
  ├── call_tool(name, args) → result  // 调用工具
  └── list_resources() / read_resource(uri)

StdioTransport    — 子进程 stdio 通信，JSON-RPC 协议
SSETransport      — HTTP Server-Sent Events，端点发现
HTTPTransport     — RESTful POST 调用
```

**MCP服务器管理器：** (详见 `backend/tools/mcp_manager.py`)

```
MCPServerManager
  ├── connect_server(config)           // 根据 transport 创建适配器、连接、发现工具
  ├── get_tool(name) → ToolAdapter    // 命名空间隔离: server_name__tool_name
  ├── list_all_tools() → [{name, server, description}]
  └── disconnect_all()
```

### 4.5 Skills 系统设计

Skills 是框架的核心扩展机制，每个 Skill 是一个包含元信息和指令文档的自包含目录。Skill 被实现为一个内置工具（`load_skill`），与 read/write/bash 等工具处于同一层级。

**核心设计理念：按需加载** — Agent 启动时只暴露 skill 摘要列表，LLM 需要时通过 `load_skill` 工具调用加载完整内容。

**Skill 目录结构：**

```
~/.ils4gas/skills/
├── code_review/
│   └── SKILL.md          # Skill 元信息(YAML front matter) + 指令文档(markdown body)
├── code_helper/
│   └── SKILL.md
└── material_analysis/
    └── SKILL.md
```

**SKILL.md 格式：**
```markdown
---
name: code_review
description: Review code for bugs, security, and style. Use when the user asks to review or audit code.
---

# Code Review Skill

## Workflow
1. Examine the provided code for issues.
2. Check for: bugs, security vulnerabilities, performance problems.
3. Report findings with severity levels.
```

- YAML front matter 仅含 `name` 和 `description` 两个字段
- body 部分作为 skill 的 `prompt`（LLM 加载时获取的内容）

**数据模型：** (详见 `backend/tools/builtin/skill.py`)

```
SkillMeta  (元信息)
  ├── name
  └── description

Skill  (完整定义)
  ├── meta: SkillMeta          — 从 SKILL.md YAML 解析
  ├── prompt: str              — SKILL.md body 正文
  └── directory: Path          — skill 目录的绝对路径
```

**按需加载流程：**

```
Agent 启动
  └─> register_builtin_tools()
        └─> create_skill_tool()
              ├─> SkillRegistry.get_summaries()     # 只解析 YAML 摘要，不加载全文
              └─> 返回 ToolInfo("load_skill")       # 注册到 ToolRegistry

用户提问
  └─> LLM 判断需要某个 skill
        └─> 调用 load_skill(name="code_review")
              └─> SkillRegistry.load(name)     # 解析 SKILL.md
              └─> 返回格式化内容:
                    ## Skill: code_review
                    Directory: /root/.ils4gas/skills/code_review

                    [SKILL.md 全文]

后续轮次
  └─> skill 内容已在会话历史的 tool message 中，LLM 可回溯，无需重载
```

**`load_skill` 工具设计：** (详见 `backend/tools/builtin/skill.py`)

```
ToolInfo(
  name = "load_skill"
  description = "Load a skill module's full instructions. Available skills:
                   - code_review: Review code for bugs, security, and style
                   - code_helper: 协助进行代码编写、审查和优化
                   - material_analysis: 材料科学分析助手..."
  parameters:
    name: string, enum=[所有已发现 skill 的名称列表], required
  
  execute(name):
    1. 调用 SkillRegistry.load(name) → 返回 Skill 对象
    2. 格式化输出:
         ## Skill: {skill.meta.name}
         Directory: {skill.directory}

         {skill.prompt}
    3. 返回格式化字符串给 LLM
)
```

**核心类与职责：**

| 类 | 文件 | 职责 |
|----|------|------|
| `SkillRegistry` | `tools/builtin/skill.py` | 扫描目录、解析 SKILL.md YAML front matter、按名称加载、缓存、去重 |
| `SkillMeta` | `tools/builtin/skill.py` | 元信息数据类 (name, description) |
| `Skill` | `tools/builtin/skill.py` | 完整 Skill 定义 (meta + prompt + directory) |

**集成点：** `load_skill` 与其他内置工具一样，由 `register_builtin_tools()` 统一注册到 `ToolRegistry`：

```
# backend/tools/builtin/__init__.py
def register_builtin_tools(registry):
    ...
    skill_tool = create_skill_tool()
    if skill_tool:
        registry.register(skill_tool)
```

该注册在 TUI 和 Web 模式下均自动执行（`backend/tui/screen.py` 和 `backend/main.py`）。

**关键设计决策：**

1. **Skill 即 Tool** — `load_skill` 是标准 `ToolInfo`，与其他 builtin 工具同层级，无需独立的服务/加载器/执行器层
2. **按需加载 vs 启动全量注入** — LLM 自主决定何时加载，避免 context 膨胀
3. **Tool 机制复用** — 完全复用 ReActAgent tool-calling 循环，不修改 Agent 核心代码
4. **跨轮次持久化** — skill 内容作为 tool message 保留在对话历史中，后续轮次自然可见
5. **目录路径暴露** — 工具返回 `Directory: /path/to/skill/`，LLM 可通过 read 工具读取 skill 目录下的附加文件
6. **同名冲突：后者覆盖 + 日志警告** — 多个 `skills_dirs` 中出现同名 skill 时，最后扫描的目录生效，同时输出 `logger.warning()`
7. **摘要不解析全文** — `get_summaries()` 只读 YAML front matter，`load()` 才触发完整解析

### 4.6 记忆与RAG系统设计

记忆系统分为三个层次：短期记忆（会话上下文）、长期记忆（持久化存储）、RAG检索增强（向量语义搜索）。

**短期记忆：** (详见 `backend/memory/short_term.py`)

```
ShortTermMemory(max_turns, max_tokens)
  ├── add(role, content, token_estimate)  // 追加消息，自动裁剪超限
  ├── _trim()                             // 超出 token 上限时弹出最早消息
  ├── get_messages() → [dict]             // 获取滑动窗口内消息
  └── clear()                             // 清空
```

**长期记忆：** (详见 `backend/memory/long_term.py`)

```
LongTermMemory(db_path)
  ├── _init_db() → 创建 memories + conversation_summaries 表
  ├── save(user_id, content, category, importance)
  ├── search(user_id, query, limit) → [dict]
  └── save_summary(session_id, user_id, summary)
```

**RAG检索增强：** (详见 `backend/memory/rag.py`)

```
RAGManager(persist_dir)
  ├── get_or_create_collection(name) → ChromaDB collection
  ├── index_messages(session_id, messages, embedding_fn)  // 消息向量化
  ├── search(session_id, query, embedding_fn) → [dict]    // 语义搜索
  └── index_documents(collection_name, documents, embedding_fn)

ConversationSummarizer(llm)
  └── summarize(messages, max_length) → str  // LLM 压缩长对话为摘要
```

### 4.7 A2A 多智能体通信设计

A2A (Agent-to-Agent) 协议用于多智能体之间的任务分发与协作。

**A2A协议消息格式：** (详见 `backend/a2a/protocol.py`)

```
A2AMessage
  ├── message_id, from_agent, to_agent
  ├── message_type: "task" | "result" | "query" | "notify"
  ├── content, context, priority

A2ATask
  ├── task_id, description
  ├── required_skills, required_tools
  ├── status ("pending"), assigned_agent, result
```

**多智能体编排器：** (详见 `backend/a2a/orchestrator.py`)

```
AgentOrchestrator
  ├── register_agent(name, agent, capabilities)
  ├── decompose_and_dispatch(task_description, llm) → [A2ATask]
  │     // LLM 分解任务 → 按 capability 匹配 Agent → 分发
  └── execute_tasks(tasks) → {task_id: result}
        // 并发执行已分发的子任务，聚合结果
```

### 4.8 Agent核心设计

**Agent基类：** (详见 `backend/core/agent.py`)

```
BaseAgent (ABC)
  ├── 属性: name, llm, system_prompt, tools (ToolRegistry)
  ├── 状态机: IDLE → RUNNING → PAUSED / ERROR → IDLE
  │
  ├── run(messages) → str                    // 抽象：同步运行
  ├── stream_run(messages) → AsyncIterator   // 抽象：流式运行
  ├── cancel()                               // 中断当前执行
  ├── pause() / resume()                     // 暂停/恢复（保存检查点）
  └── _execute_tool_sync(name, args) → str   // 同步执行工具
```

**SimpleAgent：** (详见 `backend/agents/simple_agent.py`)

```
SimpleAgent(BaseAgent)
  // 纯 LLM 透传，无工具调用
  run(messages) → str        // 构建 system prompt + history → LLM invoke → 返回
  stream_run(messages) → y   // 同上但流式
```

**ReActAgent：** (详见 `backend/agents/react_agent.py`)

```
ReActAgent(BaseAgent)       // 实际使用的 Agent，基于 OpenAI tool-calling
  run(messages):
    loop:
      调用 LLM (带 tools schema)
      if LLM 返回 tool_calls:
        执行每个 tool_call → 追加 assistant/tool 消息 → 继续循环
      else:
        返回文本内容 → 结束
  stream_run(messages):
    同上，但流式产出 AgentEvent (content_chunk | tool_call_start | tool_call_end | done)
```

**ReflectionAgent / PlanSolveAgent：** (详见 `backend/agents/reflection_agent.py`, `plan_solve_agent.py`)
当前为 stub（`raise NotImplementedError`），留待后续实现。

**事件系统：** (详见 `backend/core/events.py`)

```
AgentEventType enum: CONTENT_CHUNK | REASONING_CHUNK | TOOL_CALL_START | TOOL_CALL_END | DONE | ERROR
AgentEvent: {type, data} → to_dict() | to_sse() | to_json()
```

### 4.9 异步工具执行器设计

**工具权限分级：** (详见 `backend/tools/permission.py`)

```
ToolPermission: READ_ONLY | WRITE | EXECUTE | NETWORK
  // 四级权限：只读 < 写入 < 网络 < 执行
```

**异步执行器：** (详见 `backend/tools/async_executor.py`)

```
AsyncToolExecutor
  ├── register_tool_config(name, timeout, retries, permission)
  ├── execute(tool_name, tool_fn, args) → result
  │     // 带超时控制 + 自动重试 + 信号量并发限制
  └── cancel(tool_name)  // 取消正在执行的任务

当前状态: async_executor.py + permission.py 尚未实现（空文件）
```

### 4.10 Agent评估体系设计

**运行指标：** (详见 `backend/evaluation/metrics.py`)

```
AgentMetrics: task_id, success, steps_taken, tools_called, total_time, llm_calls, total_tokens, errors
```

**评估器：** (详见 `backend/evaluation/`)

```
AgentEvaluator
  ├── evaluate(metrics) → {task_completion, efficiency, tool_usage, error_rate}
  └── evaluate_answer_quality(question, answer, expected) → {accuracy, completeness, ...}
        // 使用 LLM 对回答质量进行 1-10 评分

当前状态: evaluation/ 目录下所有文件为空，尚未实现
```

### 4.11 LLM服务设计

**LLMService：** (详见 `backend/services/llm_service.py`)

```
LLMService(config_path, tool_registry)
  ├── _load_config() → dict         // 读取 ~/.ils4gas/config.json，解析 ${ENV:VAR}
  ├── get_provider(model_id) → LLMProvider  // 创建/缓存 Provider 实例
  ├── get_openai_tools() → [dict]   // 返回 ToolRegistry 中所有工具的 OpenAI schema
  ├── switch_model(model_id)         // 切换当前模型
  ├── list_models() → [{id, name, provider, limit}]
  ├── get_current_model() → dict
  ├── invoke / stream_invoke / ainvoke / astream_invoke(messages)  // 委托到 Provider
  └── tool_registry: ToolRegistry    // 共享的工具注册表
```

---

## 5. API设计

### 5.1 RESTful API

```
# LLM管理
GET    /api/v1/llm/providers            # 获取所有LLM提供商
GET    /api/v1/llm/models               # 获取所有模型
POST   /api/v1/llm/switch              # 切换当前模型
GET    /api/v1/llm/current             # 获取当前使用的模型
GET    /api/v1/llm/models/{id}/info    # 获取特定模型信息

# 会话管理
GET    /api/v1/sessions                 # 获取会话列表
POST   /api/v1/sessions                 # 创建新会话
GET    /api/v1/sessions/{id}            # 获取特定会话（含消息）
PUT    /api/v1/sessions/{id}            # 更新会话（标题、模型等）
DELETE /api/v1/sessions/{id}            # 删除会话
POST   /api/v1/sessions/{id}/switch-model  # 切换会话模型

# 聊天交互
POST   /api/v1/chat/{session_id}/send          # 发送消息（非流式）
POST   /api/v1/chat/{session_id}/stream        # 发送消息（SSE流式）
GET    /api/v1/chat/{session_id}/history       # 获取聊天历史
POST   /api/v1/chat/{session_id}/cancel        # 取消当前生成
POST   /api/v1/chat/{session_id}/pause         # 暂停执行（保存检查点）
POST   /api/v1/chat/{session_id}/resume        # 从检查点恢复

# MCP工具管理
GET    /api/v1/mcp/servers              # 获取MCP服务器列表
POST   /api/v1/mcp/servers/connect      # 连接MCP服务器
POST   /api/v1/mcp/servers/disconnect   # 断开MCP服务器
GET    /api/v1/mcp/tools                # 获取所有可用MCP工具
GET    /api/v1/mcp/tools/{name}         # 获取特定工具信息
POST   /api/v1/mcp/tools/{name}/execute # 直接执行工具（调试用）

# Skill管理
GET    /api/v1/skills                   # 获取所有Skill列表

# 记忆管理
GET    /api/v1/memory/long-term         # 获取长期记忆
POST   /api/v1/memory/long-term         # 添加长期记忆
DELETE /api/v1/memory/long-term/{id}    # 删除长期记忆
POST   /api/v1/memory/search            # 语义搜索历史对话
POST   /api/v1/memory/summarize         # 生成对话摘要

# Agent评估
POST   /api/v1/eval/run                 # 运行评估任务
GET    /api/v1/eval/results/{task_id}   # 获取评估结果
GET    /api/v1/eval/metrics             # 获取全局指标统计

# 工作区/上下文
GET    /api/v1/workspace/context        # 获取当前工作区上下文
PUT    /api/v1/workspace/file/{name}    # 更新上下文文件（AGENT.md等）
```

### 5.2 WebSocket API & SSE

框架同时支持 WebSocket（首选）和 SSE（降级 fallback）。

**WebSocket API:**
```
WebSocket: /api/v1/ws/chat?session_id={id}&token={token}

# 客户端 → 服务器
{
  "type": "message",
  "content": "你好"
}

# 客户端 → 服务器（取消生成）
{
  "type": "cancel"
}

# 服务器 → 客户端（推理/思考过程）
{
  "type": "thinking",
  "content": "让我分析一下这个问题..."
}

# 服务器 → 客户端（流式内容）
{
  "type": "content",
  "content": "好的，我来帮你"
}

# 服务器 → 客户端（工具调用开始）
{
  "type": "tool_call_start",
  "tool_name": "read_file",
  "tool_args": {"file_path": "main.py"},
  "display_args": {"file_path": "main.py"}
}

# 服务器 → 客户端（工具调用完成）
{
  "type": "tool_call_end",
  "tool_name": "read_file",
  "result": "import os\n...",
  "duration_ms": 234,
  "is_error": false
}

# 服务器 → 客户端（完成）
{
  "type": "done",
  "message_id": "msg_xyz789",
  "total_tokens": 1234,
  "tool_calls": 3,
  "metrics": {"steps": 3, "time_s": 5.2}
}

# 服务器 → 客户端（错误）
{
  "type": "error",
  "code": "RATE_LIMITED",
  "message": "请求过于频繁，请稍后重试"
}
```

**SSE 降级方案（当 WebSocket 不可用）：**
```
POST /api/v1/chat/{session_id}/stream

Response (text/event-stream):
event: thinking
data: {"content": "让我思考一下..."}

event: content
data: {"content": "你好"}

event: tool_call_start
data: {"tool_name": "read_file", "args": {...}}

event: tool_call_end
data: {"tool_name": "read_file", "result": "...", "duration_ms": 234}

event: done
data: {"message_id": "msg_xyz", "total_tokens": 1234}
```

---

## 6. 前端界面设计

### 6.1 技术选型

| 技术 | 选型 | 说明 |
|------|------|------|
| **框架** | React 18 + TypeScript | 组件化开发，类型安全 |
| **构建** | Vite | 快速HMR，ESM原生支持 |
| **状态管理** | Zustand | 轻量级，无模板代码 |
| **Markdown渲染** | react-markdown + rehype-highlight | 支持代码高亮 |
| **代码高亮** | Shiki | 服务端级别的高亮质量 |
| **CSS方案** | CSS Modules + CSS Variables | 组件隔离 + 主题变量 |
| **HTTP客户端** | fetch + EventSource | 原生API，无额外依赖 |
| **虚拟滚动** | @tanstack/react-virtual | 大量消息高性能渲染 |

### 6.2 组件架构

```
App
├── AppShell (布局容器)
│   ├── Sidebar (侧边栏)
│   │   ├── Logo
│   │   ├── NewChatButton
│   │   ├── SessionList
│   │   │   └── SessionItem (可重命名、删除)
│   │   └── SettingsButton
│   │
│   ├── TopBar (顶部栏)
│   │   ├── ModelSelector (模型下拉切换)
│   │   └── ThemeToggle (亮色/暗色切换)
│   │
│   └── MainPanel (主区域)
│       ├── ChatArea (聊天区域)
│       │   ├── MessageList (虚拟滚动消息列表)
│       │   │   ├── UserMessage (用户消息气泡)
│       │   │   ├── AssistantMessage (AI消息)
│       │   │   │   ├── ThinkingBlock (推理过程，可折叠)
│       │   │   │   ├── MarkdownView (Markdown渲染)
│       │   │   │   └── ToolCallCard (工具调用卡片)
│       │   │   │       ├── ToolCallArgs (参数JSON展开)
│       │   │   │       ├── ToolCallDiff (文件编辑diff视图)
│       │   │   │       └── ToolCallResult (执行结果)
│       │   │   └── SystemMessage (系统通知)
│       │   └── ScrollAnchor (自动滚动锚点)
│       │
│       └── InputArea (输入区)
│           ├── AttachButton (附件上传)
│           ├── TextInput (多行文本输入)
│           ├── SendButton / StopButton (发送/停止)
│           └── TokenCounter (Token计数)
│
├── SkillPanel (Skill管理抽屉，可滑出)
│   ├── SkillList
│   └── SkillEditor
│
└── EvalPanel (评估面板)
    ├── EvalForm
    └── EvalResults
```

### 6.3 流式渲染设计

```typescript
// hooks/useStream.ts
interface StreamState {
  content: string;
  thinking: string;
  toolCalls: ToolCall[];
  isStreaming: boolean;
  isThinking: boolean;
}

function useChatStream(sessionId: string): StreamState & {
  send: (message: string) => void;
  cancel: () => void;
} {
  const ws = useRef<WebSocket | null>(null);
  const [state, setState] = useState<StreamState>(initialState);
  
  const connect = useCallback(() => {
    // 首选 WebSocket
    try {
      ws.current = new WebSocket(`ws://${host}/api/v1/ws/chat?session_id=${sessionId}`);
      ws.current.onmessage = handleMessage;
      ws.current.onerror = () => fallbackToSSE();  // 降级到SSE
    } catch {
      fallbackToSSE();
    }
  }, [sessionId]);
  
  const fallbackToSSE = () => {
    // SSE fallback
    const source = new EventSource(`/api/v1/chat/${sessionId}/stream`);
    source.addEventListener('content', ...);
    source.addEventListener('tool_call_start', ...);
    source.addEventListener('done', ...);
  };
  
  return { ...state, send, cancel };
}
```

### 6.4 主题系统（亮色/暗色）

```css
/* styles/theme.css */
:root {
  /* 亮色主题 */
  --bg-primary: #ffffff;
  --bg-secondary: #f7f7f8;
  --bg-tertiary: #ececf1;
  --text-primary: #1a1a2e;
  --text-secondary: #6b7280;
  --border-color: #e5e7eb;
  --accent: #6366f1;
  --accent-hover: #4f46e5;
  --user-msg-bg: #6366f1;
  --user-msg-text: #ffffff;
  --assistant-msg-bg: #f7f7f8;
  --tool-card-bg: #f0fdf4;
  --tool-card-border: #86efac;
  --code-bg: #1e1e2e;
  --error-text: #ef4444;
}

[data-theme="dark"] {
  /* 暗色主题 */
  --bg-primary: #0f0f23;
  --bg-secondary: #1a1a2e;
  --bg-tertiary: #25253e;
  --text-primary: #e2e8f0;
  --text-secondary: #94a3b8;
  --border-color: #2d2d50;
  --accent: #818cf8;
  --accent-hover: #6366f1;
  --user-msg-bg: #6366f1;
  --user-msg-text: #ffffff;
  --assistant-msg-bg: #1a1a2e;
  --tool-card-bg: #0f2f1a;
  --tool-card-border: #22c55e;
  --code-bg: #0d0d1a;
  --error-text: #f87171;
}
```

### 6.5 工具栏可视化设计

工具调用在前端以可折叠卡片形式展示，针对文件编辑类工具提供 diff 视图：

```typescript
// components/ToolCallCard/ToolCallDiff.tsx
interface DiffViewProps {
  toolName: string;
  args: Record<string, any>;
  result: string;
  isError: boolean;
  durationMs: number;
}

function ToolCallCard({ toolName, args, result, isError, durationMs }: DiffViewProps) {
  const [expanded, setExpanded] = useState(false);
  
  const isFileEdit = ['edit_file', 'write_file', 'multi_edit_file'].includes(toolName);
  
  return (
    <div className={cn(styles.card, isError && styles.cardError)}>
      <div className={styles.header} onClick={() => setExpanded(!expanded)}>
        <span className={styles.icon}>{isFileEdit ? '📝' : '🔧'}</span>
        <span className={styles.name}>{toolName}</span>
        <span className={styles.duration}>{durationMs}ms</span>
        <span className={styles.chevron}>{expanded ? '▾' : '▸'}</span>
      </div>
      
      {expanded && (
        <div className={styles.body}>
          {isFileEdit ? (
            <DiffView oldString={args.old_string} newString={args.new_string} />
          ) : (
            <pre className={styles.args}>{JSON.stringify(args, null, 2)}</pre>
          )}
          <div className={styles.result}>
            {isError ? '❌ ' : '✅ '}{truncate(result, 500)}
          </div>
        </div>
      )}
    </div>
  );
}
```


### 6.6 国际化（i18n）方案

框架采用轻量级 i18n 方案，仅依赖 React Context + JSON 字典，无需引入第三方 i18n 库。

**设计原则：**
- 默认中文、一键切换英文
- 所有文案集中在 locale JSON 文件中，组件内不使用硬编码字符串
- 支持占位符变量插值（如 `{name}`、`{count}`）

**目录结构：**
```
frontend/src/
├── locales/
│   ├── zh-CN.json    # 中文文案（默认）
│   ├── en-US.json    # 英文文案
│   └── index.ts      # 导出类型定义
├── hooks/
│   └── useLocale.ts  # 语言切换 Hook
└── components/
    └── common/
        └── LocaleToggle.tsx  # 语言切换按钮
```

**文案字典（locales/zh-CN.json）：**
```json
{
  "app": { "title": "ILS4GAS 智能体", "version": "v2.0" },
  "sidebar": {
    "newChat": "新建对话",
    "searchPlaceholder": "搜索会话...",
    "noSessions": "暂无会话"
  },
  "chat": {
    "inputPlaceholder": "输入消息... (Enter 发送, Shift+Enter 换行)",
    "send": "发送",
    "stop": "停止生成",
    "regenerate": "重新生成",
    "copy": "复制",
    "copied": "已复制",
    "thinking": "思考中..."
  },
  "model": {
    "selector": "切换模型",
    "context": "上下文长度",
    "output": "最大输出"
  },
  "tools": {
    "callTitle": "工具调用",
    "executing": "执行中...",
    "result": "结果",
    "duration": "耗时",
    "error": "执行失败"
  },
  "skills": {
    "title": "技能管理",
    "load": "加载技能",
    "unload": "卸载",
    "create": "从对话创建",
    "noSkills": "暂无可用技能"
  },
  "memory": {
    "title": "记忆管理",
    "save": "保存记忆",
    "search": "搜索记忆",
    "empty": "暂无记忆"
  },
  "settings": {
    "title": "设置",
    "theme": "主题",
    "light": "亮色",
    "dark": "暗色",
    "language": "语言",
    "chinese": "中文",
    "english": "English"
  }
}
```

**文案字典（locales/en-US.json）：**
```json
{
  "app": { "title": "ILS4GAS", "version": "v2.0" },
  "sidebar": {
    "newChat": "New Chat",
    "searchPlaceholder": "Search sessions...",
    "noSessions": "No sessions"
  },
  "chat": {
    "inputPlaceholder": "Type a message... (Enter to send, Shift+Enter for new line)",
    "send": "Send",
    "stop": "Stop",
    "regenerate": "Regenerate",
    "copy": "Copy",
    "copied": "Copied",
    "thinking": "Thinking..."
  },
  "model": {
    "selector": "Switch Model",
    "context": "Context Length",
    "output": "Max Output"
  },
  "tools": {
    "callTitle": "Tool Call",
    "executing": "Executing...",
    "result": "Result",
    "duration": "Duration",
    "error": "Execution Failed"
  },
  "skills": {
    "title": "Skills",
    "load": "Load Skill",
    "unload": "Unload",
    "create": "Create from chat",
    "noSkills": "No skills available"
  },
  "memory": {
    "title": "Memory",
    "save": "Save Memory",
    "search": "Search Memory",
    "empty": "No memories"
  },
  "settings": {
    "title": "Settings",
    "theme": "Theme",
    "light": "Light",
    "dark": "Dark",
    "language": "Language",
    "chinese": "中文",
    "english": "English"
  }
}
```

**类型定义（locales/index.ts）：**
```typescript
export type Locale = 'zh-CN' | 'en-US';

export interface LocaleDict {
  app: { title: string; version: string };
  sidebar: { newChat: string; searchPlaceholder: string; noSessions: string };
  chat: { inputPlaceholder: string; send: string; stop: string; regenerate: string; copy: string; copied: string; thinking: string };
  model: { selector: string; context: string; output: string };
  tools: { callTitle: string; executing: string; result: string; duration: string; error: string };
  skills: { title: string; load: string; unload: string; create: string; noSkills: string };
  memory: { title: string; save: string; search: string; empty: string };
  settings: { title: string; theme: string; light: string; dark: string; language: string; chinese: string; english: string };
}
```

**i18n Context + Hook（hooks/useLocale.ts）：**
```typescript
import React, { createContext, useContext, useState, useCallback, useEffect } from 'react';
import { Locale, LocaleDict } from '../locales';
import zhCN from '../locales/zh-CN.json';
import enUS from '../locales/en-US.json';

const localeMap: Record<Locale, LocaleDict> = {
  'zh-CN': zhCN as LocaleDict,
  'en-US': enUS as LocaleDict,
};

interface LocaleContextValue {
  locale: Locale;
  t: LocaleDict;
  setLocale: (l: Locale) => void;
  toggleLocale: () => void;
}

const LocaleContext = createContext<LocaleContextValue>(null!);

const STORAGE_KEY = 'ils4gas_locale';

export function LocaleProvider({ children }: { children: React.ReactNode }) {
  const [locale, setLocaleState] = useState<Locale>(() => {
    const saved = localStorage.getItem(STORAGE_KEY);
    return (saved as Locale) || 'zh-CN';  // 默认中文
  });

  const setLocale = useCallback((l: Locale) => {
    setLocaleState(l);
    localStorage.setItem(STORAGE_KEY, l);
  }, []);

  const toggleLocale = useCallback(() => {
    setLocale(locale === 'zh-CN' ? 'en-US' : 'zh-CN');
  }, [locale, setLocale]);

  const value: LocaleContextValue = {
    locale,
    t: localeMap[locale],
    setLocale,
    toggleLocale,
  };

  return <LocaleContext.Provider value={value}>{children}</LocaleContext.Provider>;
}

export function useLocale() {
  return useContext(LocaleContext);
}

// 模板插值工具函数
export function interpolate(text: string, vars: Record<string, string | number>): string {
  return text.replace(/\{(\w+)\}/g, (_, key) => String(vars[key] ?? `{${key}}`));
}
```

**语言切换按钮（components/common/LocaleToggle.tsx）：**
```typescript
import { useLocale } from '../../hooks/useLocale';

export function LocaleToggle() {
  const { locale, toggleLocale } = useLocale();

  return (
    <button
      onClick={toggleLocale}
      className={styles.toggle}
      title={locale === 'zh-CN' ? 'Switch to English' : '切换中文'}
    >
      {locale === 'zh-CN' ? 'EN' : '中'}
    </button>
  );
}
```

**组件使用示例：**
```typescript
// 在任何组件中：
import { useLocale, interpolate } from '../../hooks/useLocale';

function ChatInput() {
  const { t } = useLocale();

  return (
    <div>
      <textarea placeholder={t.chat.inputPlaceholder} />
      <button>{t.chat.send}</button>
    </div>
  );
}

// 带变量的插值：
function ToolCard({ toolName, durationMs }: { toolName: string; durationMs: number }) {
  const { t } = useLocale();
  return (
    <span>{interpolate(t.tools.duration, { name: toolName, ms: durationMs })}</span>
  );
}
```

### 7.1 服务器配置

**.env 环境变量：**
```env
# 服务器绑定
HOST=127.0.0.1
PORT=8789

# 认证令牌（必须设置！）
GATEWAY_TOKEN=your-strong-secret-token-here

# 速率限制
RATE_LIMIT_ENABLED=true
RATE_LIMIT_PER_MINUTE=30

# 安全配置
ALLOWED_ORIGINS=http://localhost:5173,https://your-domain.com
MAX_REQUEST_SIZE_MB=10
TOOL_EXECUTION_TIMEOUT_S=30

# 其他
DEBUG=false
LOG_LEVEL=info
```

### 7.2 多层安全防护

**1. 认证层（Auth Middleware）：**
```python
# api/middleware/auth.py
class AuthMiddleware:
    """认证中间件 - Bearer Token / Query Param 两种方式"""
    
    EXEMPT_PATHS = {"/health", "/api/health", "/chat", "/static"}
    
    async def __call__(self, request: Request, call_next):
        if request.url.path in self.EXEMPT_PATHS:
            return await call_next(request)
        
        token = None
        auth = request.headers.get("Authorization", "")
        if auth.startswith("Bearer "):
            token = auth[7:]
        else:
            token = request.query_params.get("token")
        
        if not token or not self._verify_token(token):
            return JSONResponse(status_code=401, content={"detail": "Unauthorized"})
        
        return await call_next(request)
```

**2. 速率限制（Rate Limit Middleware）：**
```python
# api/middleware/ratelimit.py
from collections import defaultdict
import time

class RateLimitMiddleware:
    """基于滑动窗口的速率限制"""
    
    def __init__(self, max_requests: int = 30, window_seconds: int = 60):
        self.max_requests = max_requests
        self.window = window_seconds
        self._clients: dict[str, list] = defaultdict(list)
    
    async def __call__(self, request: Request, call_next):
        client_ip = request.client.host
        now = time.time()
        
        # 清理过期记录
        self._clients[client_ip] = [
            t for t in self._clients[client_ip] 
            if now - t < self.window
        ]
        
        if len(self._clients[client_ip]) >= self.max_requests:
            return JSONResponse(
                status_code=429,
                content={"detail": "Too many requests. 请稍后重试。"}
            )
        
        self._clients[client_ip].append(now)
        return await call_next(request)
```

**3. 工具权限分级：**
```python
# tools/permission.py
from enum import Enum

class ToolPermission(Enum):
    READ_ONLY = "read_only"      # 只读：文件读取、搜索
    WRITE = "write"             # 写入：文件编辑、创建
    NETWORK = "network"         # 网络：HTTP请求
    EXECUTE = "execute"         # 执行：代码运行、系统命令

class ToolPermissionManager:
    """工具权限管理器"""
    
    def __init__(self):
        self._tool_permissions: dict[str, ToolPermission] = {}
        self._session_permissions: dict[str, set[ToolPermission]] = {}
    
    def register_tool(self, tool_name: str, permission: ToolPermission):
        self._tool_permissions[tool_name] = permission
    
    def grant_session(self, session_id: str, permissions: set[ToolPermission]):
        """授予会话工具权限"""
        self._session_permissions[session_id] = permissions
    
    def can_use_tool(self, session_id: str, tool_name: str) -> bool:
        """检查会话是否有权限使用工具"""
        if session_id not in self._session_permissions:
            return False
        required = self._tool_permissions.get(tool_name, ToolPermission.READ_ONLY)
        return required in self._session_permissions[session_id]
    
    def require_confirmation(self, tool_name: str) -> bool:
        """高风险工具需要用户确认"""
        permission = self._tool_permissions.get(tool_name, ToolPermission.READ_ONLY)
        return permission in (ToolPermission.EXECUTE, ToolPermission.WRITE)
```

**4. 输入过滤：**
```python
# api/middleware/input_filter.py
import re

class InputFilter:
    """输入内容安全检查"""
    
    MAX_LENGTH = 100000  # 最大输入长度
    
    @staticmethod
    def sanitize(text: str) -> str:
        """过滤潜在危险内容"""
        if len(text) > InputFilter.MAX_LENGTH:
            raise ValueError(f"输入超过最大长度限制 ({InputFilter.MAX_LENGTH})")
        # 移除零宽字符（防止隐藏注入）
        text = re.sub(r'[\u200b-\u200f\u2028-\u202f]', '', text)
        return text
```

### 7.3 远程部署安全

1. **生产环境建议**：Nginx 反向代理 + Let's Encrypt SSL
2. **容器隔离**：Docker 部署，工具执行在独立容器中
3. **日志审计**：记录所有工具调用和模型交互
4. **定期更新**：保持依赖包和系统的安全更新
5. **最小权限原则**：运行时使用非 root 用户

### 7.4 Nginx反向代理配置

```nginx
server {
    listen 443 ssl;
    server_name your-domain.com;
    
    ssl_certificate /etc/letsencrypt/live/your-domain.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/your-domain.com/privkey.pem;
    
    # 速率限制
    limit_req_zone $binary_remote_addr zone=api_limit:10m rate=30r/m;
    
    location / {
        limit_req zone=api_limit burst=10 nodelay;
        proxy_pass http://localhost:8789;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # WebSocket支持
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        
        # 请求大小限制
        client_max_body_size 10M;
    }
}
```

---

## 8. 实施计划

### 阶段一：核心框架 MVP（第1-2周）

1. **项目初始化**
   - 创建目录结构，初始化虚拟环境
   - 配置 `pyproject.toml`、`requirements.txt`
   - 创建 `config.example.json` 模板文件

2. **核心框架层**
   - `core/llm.py` - Provider抽象 + OpenAI兼容实现
   - `core/message.py` - 消息系统
   - `core/config.py` - 配置管理（热加载）
   - `core/context.py` - 上下文工程（AGENT.md等）
   - `core/agent.py` - Agent基类（含中断/恢复/事件总线）

3. **Agent实现**
   - `agents/simple_agent.py` - 基础对话Agent
   - `agents/react_agent.py` - ReAct推理-行动Agent
   - `agents/reflection_agent.py` - 反思优化Agent

4. **服务层基础**
   - `services/llm_service.py` - LLM管理
   - `services/session_service.py` - SQLite会话管理

### 阶段二：工具与记忆系统（第3-4周）

5. **工具系统**
   - `tools/base.py` - 工具基类（同步/异步统一）
   - `tools/registry.py` - 工具注册与发现
   - `tools/mcp_adapter.py` - MCP适配器（stdio/SSE/HTTP）
   - `tools/mcp_manager.py` - MCP服务器生命周期管理
   - `tools/async_executor.py` - 异步执行器（超时/重试/并发）
   - `tools/permission.py` - 工具权限分级
   - `tools/builtin/` - 内置工具集

6. **Skills 系统**
   - `tools/builtin/skill.py` - Skill 注册、发现、按需加载

7. **记忆系统**
   - `memory/short_term.py` - 短期记忆（滑动窗口）
   - `memory/long_term.py` - 长期记忆（SQLite）
   - `memory/vector_store.py` - 向量存储（ChromaDB）
   - `memory/rag.py` - RAG检索增强

### 阶段三：API与前端（第5-7周）

8. **API层**
   - `api/routes/chat.py` - 聊天路由（含流式）
   - `api/routes/session.py` - 会话管理
   - `api/routes/llm.py` - LLM管理
   - `api/routes/mcp.py` - MCP工具管理
   - `api/routes/skill.py` - Skill管理
   - `api/routes/memory.py` - 记忆管理
   - `api/routes/ws.py` - WebSocket + SSE

9. **中间件**
   - `api/middleware/auth.py` - 认证
   - `api/middleware/ratelimit.py` - 速率限制
   - `api/middleware/input_filter.py` - 输入过滤

10. **前端开发**
    - 项目搭建：Vite + React + TypeScript
    - 核心组件：AppShell、ChatArea、MessageList、InputArea
    - 流式渲染：WebSocket + SSE fallback
    - 工具可视化：ToolCallCard（diff视图、折叠展开）
    - 主题系统：CSS Variables 亮色/暗色切换
    - Skill面板、评估面板

### 阶段四：高级功能与评估（第8-9周）

11. **A2A通信**
    - `a2a/protocol.py` - A2A消息协议
    - `a2a/orchestrator.py` - 多智能体编排
    - `a2a/router.py` - 任务路由

12. **评估体系**
    - `evaluation/metrics.py` - 评估指标
    - `evaluation/benchmarks.py` - 基准测试

13. **测试**
    - 单元测试：Agent、Tools、Memory、Skills
    - 集成测试：API端点、WebSocket通信
    - 端到端测试：完整对话流程

### 阶段五：部署与文档（第10周）

14. **部署准备**
    - Dockerfile + docker-compose.yml
    - Nginx配置模板
    - 健康检查端点

15. **文档完善**
    - README.md、API文档（FastAPI自动生成）
    - 用户手册、开发指南
    - Skill编写指南

---

## 9. 关键依赖

**requirements.txt:**
```
# Web框架
fastapi>=0.115.0
uvicorn[standard]>=0.32.0

# OpenAI兼容接口
openai>=1.60.0

# MCP支持
mcp>=1.4.0

# 数据处理
pydantic>=2.10.0
pydantic-settings>=2.7.0

# 环境变量
python-dotenv>=1.0.0

# 向量数据库
chromadb>=0.5.0

# LLM扩展（可选，支持100+提供商）
litellm>=1.55.0

# 嵌入模型（可选，本地模式）
sentence-transformers>=3.3.0

# 日志
loguru>=0.7.0

# 定时任务
apscheduler>=3.10.0

# HTTP客户端（MCP SSE传输）
aiohttp>=3.10.0

# YAML解析（Skill配置）
pyyaml>=6.0.0

# 开发依赖
pytest>=8.0.0
pytest-asyncio>=0.24.0
httpx>=0.28.0
black>=24.0.0
ruff>=0.7.0
```

**package.json (前端):**
```json
{
  "dependencies": {
    "react": "^18.3.0",
    "react-dom": "^18.3.0",
    "react-markdown": "^9.0.0",
    "rehype-highlight": "^7.0.0",
    "zustand": "^5.0.0",
    "@tanstack/react-virtual": "^3.10.0",
    "shiki": "^1.20.0"
  },
  "devDependencies": {
    "typescript": "^5.6.0",
    "vite": "^6.0.0",
    "@types/react": "^18.3.0",
    "@types/react-dom": "^18.3.0"
  }
}
```

---

## 10. 总结

本技术设计文档详细描述了一个功能完备的智能体框架开发方案，核心特点包括：

✅ **轻量级架构** - 借鉴hello-agents的设计理念，代码简洁易维护  
✅ **多LLM支持** - Provider抽象层，支持动态切换，切换时自动刷新Agent缓存  
✅ **MCP工具深度集成** - 支持stdio/SSE/HTTP三种传输协议，多服务器管理，工具懒加载与权限分级  
✅ **Skills技能系统** - 可复用的专业化技能模块，支持自动创建（从对话中学习）  
✅ **记忆与RAG** - 三级记忆体系（短期/长期/向量语义搜索），对话摘要自动压缩  
✅ **上下文工程** - Workspace配置层（AGENT.md/PERSONA.md/MEMORY.md），多用户隔离  
✅ **多种Agent范式** - SimpleAgent / ReActAgent / ReflectionAgent，支持中断/暂停/恢复  
✅ **A2A多智能体** - Agent间任务分发与协作编排  
✅ **评估体系** - 内置指标计算与LLM辅助质量评估  
✅ **会话管理** - SQLite持久化，支持并发读写和模型切换  
✅ **前端组件化** - React+TypeScript，流式渲染（WebSocket+SSE fallback），暗色主题  
✅ **多层安全** - Bearer认证、速率限制、工具权限分级、输入过滤  
✅ **远程部署** - Docker容器化，Nginx反向代理，HTTPS支持  
✅ **完全可控** - 每一行代码都由开发者掌控，无黑盒  

通过这个方案，你将拥有一个功能完整、易于扩展的智能体框架，可以根据实际需求灵活调整和优化。

---

**文档版本：** v2.0  
**创建日期：** 2024-01-15  
**最后更新：** 2025-04-30  
**作者：** AI Assistant
