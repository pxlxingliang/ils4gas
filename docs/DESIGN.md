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
│  LLMService | SessionService | MCPService | SkillService          │
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
│   │   └── builtin/
│   │       ├── file_tools.py
│   │       ├── web_tools.py
│   │       ├── code_tools.py
│   │       └── system_tools.py
│   │
│   ├── skills/                  # Skills系统层
│   │   ├── __init__.py
│   │   ├── registry.py
│   │   ├── loader.py
│   │   ├── executor.py
│   │   ├── creator.py
│   │   └── validator.py
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
│   │   ├── skill_service.py
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

**上下文加载流程：**
```python
# core/context.py
from pathlib import Path
from typing import Optional, Dict

class WorkspaceContext:
    """工作区上下文管理器"""
    
    def __init__(self, workspace_path: Path):
        self.workspace_path = workspace_path
        self._files = {
            "agent": workspace_path / "AGENT.md",
            "persona": workspace_path / "PERSONA.md",
            "memory": workspace_path / "MEMORY.md",
            "tools": workspace_path / "TOOLS.md",
            "skills": workspace_path / "SKILLS.md",
        }
    
    def load_context(self) -> Dict[str, Optional[str]]:
        """加载所有上下文文件"""
        context = {}
        for key, path in self._files.items():
            if path.exists():
                context[key] = path.read_text(encoding="utf-8")
            else:
                context[key] = None
        return context
    
    def build_system_prompt(self) -> str:
        """构建系统提示词（合并上下文文件）"""
        ctx = self.load_context()
        parts = []
        if ctx.get("persona"):
            parts.append(f"## 角色定义\n{ctx['persona']}")
        if ctx.get("agent"):
            parts.append(f"## 行为准则\n{ctx['agent']}")
        if ctx.get("memory"):
            parts.append(f"## 历史记忆\n{ctx['memory']}")
        if ctx.get("tools"):
            parts.append(f"## 工具使用说明\n{ctx['tools']}")
        return "\n\n".join(parts)
    
    def update_memory(self, content: str):
        """追加长期记忆"""
        with open(self._files["memory"], "a", encoding="utf-8") as f:
            f.write(f"\n- {content}")
    
    def write_file(self, name: str, content: str):
        """写入上下文文件"""
        path = self.workspace_path / f"{name.upper()}.md"
        path.write_text(content, encoding="utf-8")
```

**上下文加载时机：**
- Agent 初始化时加载一次作为 system prompt
- 用户可通过对话指令动态修改（如 "请记住我喜欢简洁的回答" → 写入 MEMORY.md）
- 支持多用户隔离（每个用户有独立的 workspace 子目录）

### 4.2 LLM配置设计（多Provider抽象）

**Provider 抽象层设计：**

```python
# core/llm.py
from abc import ABC, abstractmethod
from typing import Iterator, List, Dict, Any, Optional
from openai import OpenAI, AsyncOpenAI

class LLMProvider(ABC):
    """LLM Provider 抽象基类"""
    
    @abstractmethod
    def invoke(self, messages: List[Dict], **kwargs) -> str:
        """同步调用LLM"""
        pass
    
    @abstractmethod
    def stream_invoke(self, messages: List[Dict], **kwargs) -> Iterator[str]:
        """流式调用LLM"""
        pass
    
    @abstractmethod
    async def ainvoke(self, messages: List[Dict], **kwargs) -> str:
        """异步调用LLM"""
        pass
    
    @abstractmethod
    async def astream_invoke(self, messages: List[Dict], **kwargs):
        """异步流式调用LLM"""
        pass
    
    @property
    @abstractmethod
    def model_name(self) -> str:
        pass
    
    @property
    @abstractmethod
    def provider_name(self) -> str:
        pass
    
    @property
    @abstractmethod
    def context_limit(self) -> int:
        pass

class OpenAICompatibleProvider(LLMProvider):
    """OpenAI 兼容接口 Provider（适用于 OpenAI/火山引擎/DeepSeek 等）"""
    
    def __init__(self, api_key: str, base_url: Optional[str], model_name: str, 
                 context_limit: int = 128000):
        self.client = OpenAI(api_key=api_key, base_url=base_url)
        self.async_client = AsyncOpenAI(api_key=api_key, base_url=base_url)
        self._model_name = model_name
        self._context_limit = context_limit
    
    def invoke(self, messages, **kwargs):
        response = self.client.chat.completions.create(
            model=self._model_name, messages=messages, **kwargs
        )
        return response.choices[0].message.content or ""
    
    def stream_invoke(self, messages, **kwargs):
        stream = self.client.chat.completions.create(
            model=self._model_name, messages=messages, stream=True, **kwargs
        )
        for chunk in stream:
            if chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content
    
    async def ainvoke(self, messages, **kwargs):
        response = await self.async_client.chat.completions.create(
            model=self._model_name, messages=messages, **kwargs
        )
        return response.choices[0].message.content or ""
    
    async def astream_invoke(self, messages, **kwargs):
        stream = await self.async_client.chat.completions.create(
            model=self._model_name, messages=messages, stream=True, **kwargs
        )
        async for chunk in stream:
            if chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content
    
    @property
    def model_name(self): return self._model_name
    @property
    def provider_name(self): return "openai-compatible"
    @property
    def context_limit(self): return self._context_limit
```

**~/.ils4gas/config.json（合并配置，包含模型+MCP+服务）：**
```json
{
  "$schema": "./models.schema.json",
  "currentModel": "volcengine/deepseek-v4-pro",
  "providers": {
    "volcengine": {
      "name": "Volcano Engine",
      "type": "openai-compatible",
      "options": {
        "baseURL": "https://ark.cn-beijing.volces.com/api/coding/v3",
        "apiKey": "${ENV:VOLC_API_KEY}"
      },
      "models": {
        "deepseek-v4-pro": {
          "name": "DeepSeek V4 Pro",
          "limit": { "context": 256000, "output": 4096 },
          "modalities": { "input": ["text", "image"], "output": ["text"] },
          "supports_tools": true,
          "supports_reasoning": true
        },
        "glm-4.7": {
          "name": "GLM-4.7",
          "limit": { "context": 200000, "output": 4096 },
          "supports_tools": true
        },
        "deepseek-v3.2": {
          "name": "DeepSeek V3.2",
          "limit": { "context": 128000, "output": 4096 },
          "supports_tools": true
        }
      }
    },
    "openai": {
      "name": "OpenAI",
      "type": "openai",
      "options": { "apiKey": "${ENV:OPENAI_API_KEY}" },
      "models": {
        "gpt-4.1": {
          "name": "GPT-4.1",
          "limit": { "context": 1000000, "output": 32768 },
          "supports_tools": true,
          "supports_reasoning": false
        }
      }
    },
    "anthropic": {
      "name": "Anthropic",
      "type": "anthropic",
      "options": { "apiKey": "${ENV:ANTHROPIC_API_KEY}" },
      "models": {
        "claude-opus-4-20250514": {
          "name": "Claude Opus 4",
          "limit": { "context": 200000, "output": 32000 },
          "supports_tools": true,
          "supports_reasoning": true
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

**会话服务设计（SQLite适配）：**
```python
# services/session_service.py
import sqlite3
import json
from datetime import datetime, timezone
from typing import List, Optional
from uuid import uuid4

class SessionService:
    """会话管理服务（SQLite后端，支持并发读写）"""
    
    def __init__(self, db_path: str = "~/.ils4gas/data/sessions.db"):
        self.db_path = db_path
        self._init_db()
    
    def _get_conn(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA journal_mode=WAL")   # 支持并发读
        conn.execute("PRAGMA foreign_keys=ON")
        return conn
    
    def _init_db(self):
        with self._get_conn() as conn:
            conn.executescript("""
                CREATE TABLE IF NOT EXISTS sessions (...);
                CREATE TABLE IF NOT EXISTS messages (...);
                CREATE INDEX IF NOT EXISTS idx_messages_session ON messages(session_id, timestamp);
                CREATE INDEX IF NOT EXISTS idx_sessions_updated ON sessions(updated_at DESC);
            """)
    
    def create_session(self, model_provider: str, model_name: str, 
                       title: str = "新对话") -> dict:
        session_id = f"sess_{uuid4().hex[:12]}"
        now = datetime.now(timezone.utc).isoformat()
        with self._get_conn() as conn:
            conn.execute(
                "INSERT INTO sessions VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (session_id, title, now, now, model_provider, model_name, "{}", 1)
            )
        return {"id": session_id, "title": title, "created_at": now}
    
    def list_sessions(self, limit: int = 50) -> List[dict]:
        with self._get_conn() as conn:
            rows = conn.execute(
                "SELECT * FROM sessions WHERE is_active=1 ORDER BY updated_at DESC LIMIT ?",
                (limit,)
            ).fetchall()
        return [dict(r) for r in rows]
    
    def add_message(self, session_id: str, role: str, content: str, 
                    tool_calls: list = None) -> dict:
        msg_id = f"msg_{uuid4().hex[:12]}"
        now = datetime.now(timezone.utc).isoformat()
        with self._get_conn() as conn:
            conn.execute(
                "INSERT INTO messages VALUES (?, ?, ?, ?, ?, ?)",
                (msg_id, session_id, role, content, 
                 json.dumps(tool_calls) if tool_calls else None, now)
            )
            conn.execute(
                "UPDATE sessions SET updated_at=? WHERE id=?",
                (now, session_id)
            )
        return {"id": msg_id, "role": role, "content": content}

    def get_messages(self, session_id: str, limit: int = 100) -> List[dict]:
        with self._get_conn() as conn:
            rows = conn.execute(
                "SELECT * FROM messages WHERE session_id=? ORDER BY timestamp ASC LIMIT ?",
                (session_id, limit)
            ).fetchall()
        return [dict(r) for r in rows]
    
    def switch_model(self, session_id: str, provider: str, model: str):
        """切换会话的模型配置"""
        with self._get_conn() as conn:
            conn.execute(
                "UPDATE sessions SET model_provider=?, model_name=?, updated_at=? WHERE id=?",
                (provider, model, datetime.now(timezone.utc).isoformat(), session_id)
            )
```

### 4.4 MCP工具集成方案（增强版）

MCP (Model Context Protocol) 是本框架的核心工具集成协议。支持三种传输方式，多服务器管理，以及工具懒加载。

**MCP服务器配置（在 ~/.ils4gas/config.json 中）：**
```json
{
  "mcp_servers": [
    {
      "name": "filesystem",
      "description": "本地文件系统操作",
      "transport": "stdio",
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-filesystem", "/workspace"],
      "env": {},
      "auto_connect": true,
      "tool_allowlist": null,
      "tool_denylist": ["rm", "delete"],
      "timeout_ms": 30000
    },
    {
      "name": "web_search",
      "description": "网页搜索服务",
      "transport": "sse",
      "url": "http://localhost:3001/sse",
      "auto_connect": true,
      "timeout_ms": 15000
    },
    {
      "name": "custom_python",
      "description": "自定义Python工具服务",
      "transport": "stdio",
      "command": "python",
      "args": ["my_mcp_server.py"],
      "auto_connect": false,
      "timeout_ms": 60000
    }
  ]
}
```

**MCP适配器设计（支持三种传输协议）：**

```python
# tools/mcp_adapter.py
import asyncio
import json
from typing import Any, Dict, List, Optional, Callable
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

@dataclass
class MCPServerConfig:
    """MCP服务器配置"""
    name: str
    description: str = ""
    transport: str = "stdio"  # stdio | sse | http
    command: Optional[str] = None     # stdio 模式
    args: List[str] = field(default_factory=list)
    url: Optional[str] = None         # sse/http 模式
    env: Dict[str, str] = field(default_factory=dict)
    auto_connect: bool = True
    tool_allowlist: Optional[List[str]] = None  # 允许的工具列表
    tool_denylist: Optional[List[str]] = None   # 禁止的工具列表
    timeout_ms: int = 30000

class MCPTransport(ABC):
    """MCP传输层抽象"""
    
    @abstractmethod
    async def connect(self) -> None: ...
    
    @abstractmethod
    async def disconnect(self) -> None: ...
    
    @abstractmethod
    async def list_tools(self) -> List[dict]: ...
    
    @abstractmethod
    async def call_tool(self, tool_name: str, arguments: dict) -> Any: ...
    
    @abstractmethod
    async def list_resources(self) -> List[dict]: ...
    
    @abstractmethod
    async def read_resource(self, uri: str) -> Any: ...

class StdioTransport(MCPTransport):
    """stdio传输（子进程通信）"""
    
    async def connect(self) -> None:
        """启动子进程并建立MCP会话"""
        self.process = await asyncio.create_subprocess_exec(
            self.config.command, *self.config.args,
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            env={**os.environ, **self.config.env}
        )
        # 发送 initialize 请求，建立 MCP 会话...
    
    async def call_tool(self, tool_name: str, arguments: dict) -> Any:
        """通过stdio发送工具调用请求"""
        request = json.dumps({
            "jsonrpc": "2.0",
            "method": "tools/call",
            "params": {"name": tool_name, "arguments": arguments},
            "id": 1
        })
        self.process.stdin.write((request + "\n").encode())
        self.process.stdin.flush()
        # 读取响应...
        return result

class SSETransport(MCPTransport):
    """SSE传输（HTTP Server-Sent Events）"""
    
    async def connect(self) -> None:
        """建立SSE连接"""
        self.session = aiohttp.ClientSession()
        self.event_stream = await self.session.get(
            f"{self.config.url}/sse",
            headers={"Accept": "text/event-stream"}
        )
        # 处理端点发现...

class HTTPTransport(MCPTransport):
    """HTTP传输（RESTful）"""
    
    async def call_tool(self, tool_name: str, arguments: dict) -> Any:
        response = await self.session.post(
            f"{self.config.url}/tools/call",
            json={"name": tool_name, "arguments": arguments}
        )
        return await response.json()

class MCPToolAdapter:
    """MCP工具适配器（将MCP工具统一为标准工具接口）"""
    
    def __init__(self, transport: MCPTransport, tool_schema: dict):
        self._transport = transport
        self._schema = tool_schema
    
    @property
    def name(self) -> str:
        return self._schema["name"]
    
    @property
    def description(self) -> str:
        return self._schema.get("description", "")
    
    @property
    def parameters(self) -> dict:
        return self._schema.get("inputSchema", {})
    
    async def run(self, **kwargs) -> str:
        """异步执行MCP工具（带超时控制）"""
        try:
            result = await asyncio.wait_for(
                self._transport.call_tool(self.name, kwargs),
                timeout=self._transport.config.timeout_ms / 1000
            )
            return json.dumps(result, ensure_ascii=False)
        except asyncio.TimeoutError:
            return f"工具 '{self.name}' 执行超时"
        except Exception as e:
            return f"工具执行失败: {str(e)}"
```

**MCP服务器生命周期管理（支持懒加载）：**

```python
# tools/mcp_manager.py
class MCPServerManager:
    """管理多个MCP服务器的连接与工具发现"""
    
    def __init__(self):
        self._servers: Dict[str, MCPTransport] = {}
        self._tools: Dict[str, MCPToolAdapter] = {}
        self._tool_to_server: Dict[str, str] = {}  # tool_name → server_name
    
    async def connect_server(self, config: MCPServerConfig) -> None:
        """连接MCP服务器并发现工具（懒加载：仅在使用时连接）"""
        transport = self._create_transport(config)
        await transport.connect()
        self._servers[config.name] = transport
        
        # 发现工具并注册
        tools = await transport.list_tools()
        for tool_schema in tools:
            name = tool_schema["name"]
            # 命名空间隔离：server_name__tool_name
            full_name = f"{config.name}__{name}"
            
            # 权限过滤
            if config.tool_denylist and name in config.tool_denylist:
                continue
            if config.tool_allowlist and name not in config.tool_allowlist:
                continue
            
            adapter = MCPToolAdapter(transport, tool_schema)
            self._tools[full_name] = adapter
            self._tool_to_server[full_name] = config.name
    
    def _create_transport(self, config: MCPServerConfig) -> MCPTransport:
        if config.transport == "stdio":
            return StdioTransport(config)
        elif config.transport == "sse":
            return SSETransport(config)
        elif config.transport == "http":
            return HTTPTransport(config)
        raise ValueError(f"Unsupported transport: {config.transport}")
    
    def get_tool(self, name: str) -> Optional[MCPToolAdapter]:
        """获取工具（支持 server__tool 格式或直接工具名）"""
        if name in self._tools:
            return self._tools[name]
        # 尝试匹配未加前缀的工具名
        for full_name, tool in self._tools.items():
            if full_name.endswith(f"__{name}"):
                return tool
        return None
    
    def list_all_tools(self) -> List[dict]:
        """列出所有可用工具及其来源服务器"""
        return [
            {"name": name, "server": self._tool_to_server.get(name), 
             "description": tool.description}
            for name, tool in self._tools.items()
        ]
    
    async def disconnect_all(self) -> None:
        """断开所有MCP服务器连接"""
        for server in self._servers.values():
            await server.disconnect()
        self._servers.clear()
        self._tools.clear()
```

### 4.5 Skills 系统设计

Skills 是框架的核心扩展机制，允许 Agent 动态加载和执行专业化技能模块。每个 Skill 是一个包含元信息、提示词、工具依赖和执行逻辑的自包含单元。

**Skill 结构定义：**

```
~/.ils4gas/skills/
├── code_review/
│   ├── SKILL.md          # Skill 元信息和提示词
│   ├── skill.py          # Skill 执行逻辑
│   └── tools.py          # Skill 专属工具（可选）
├── data_analysis/
│   ├── SKILL.md
│   └── skill.py
└── pdf_reader/
    ├── SKILL.md
    └── skill.py
```

**SKILL.md 格式：**
```markdown
---
name: code_review
version: "1.0"
description: 对代码进行全面的审查，包括代码质量、安全性、性能分析
author: ILS4GAS Team
tags: [code, review, quality]
dependencies:
  tools: [read_file, search_code]
  python_packages: []
trigger_keywords: [review, 审查, 检查代码, code review]
---

# Code Review Skill

## 工作流程
1. 读取目标文件
2. 分析代码结构
3. 检查常见问题（安全漏洞、性能瓶颈、代码异味）
4. 生成审查报告

## 报告格式
- **严重性**: [高/中/低]
- **问题类型**: [安全/性能/可读性/最佳实践]
- **位置**: [文件名:行号]
- **描述**: [具体问题描述]
- **建议**: [改进建议]
```

**Skill 加载器与执行引擎：**

```python
# skills/registry.py
import yaml
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, field

@dataclass
class SkillMeta:
    """Skill 元信息"""
    name: str
    version: str
    description: str
    author: str = ""
    tags: List[str] = field(default_factory=list)
    dependencies: Dict[str, List[str]] = field(default_factory=dict)
    trigger_keywords: List[str] = field(default_factory=list)

@dataclass
class Skill:
    """Skill 完整定义"""
    meta: SkillMeta
    prompt: str           # SKILL.md 正文（作为 system prompt 注入）
    module: any           # skill.py 模块
    tools: List = field(default_factory=list)  # Skill 专属工具
    
    def get_system_prompt(self) -> str:
        """生成注入到 Agent system prompt 的提示词"""
        return f"""<skill name="{self.meta.name}">
{self.prompt}
</skill>"""

class SkillRegistry:
    """Skill 注册与发现中心"""
    
    def __init__(self, skills_dirs: List[str] = None):
        self._skills: Dict[str, Skill] = {}
        self._skills_dirs = skills_dirs or ["~/.ils4gas/skills/"]
    
    def discover(self) -> List[str]:
        """扫描所有 Skill 目录，发现可用 Skill"""
        discovered = []
        for skills_dir in self._skills_dirs:
            base = Path(skills_dir)
            if not base.exists():
                continue
            for skill_dir in base.iterdir():
                if skill_dir.is_dir() and (skill_dir / "SKILL.md").exists():
                    skill_name = skill_dir.name
                    discovered.append(skill_name)
        return discovered
    
    def load(self, name: str) -> Optional[Skill]:
        """加载指定 Skill"""
        if name in self._skills:
            return self._skills[name]
        
        for skills_dir in self._skills_dirs:
            skill_dir = Path(skills_dir) / name
            if not skill_dir.exists():
                continue
            
            # 解析 SKILL.md
            meta_path = skill_dir / "SKILL.md"
            content = meta_path.read_text(encoding="utf-8")
            meta, prompt = self._parse_markdown(content)
            
            # 动态加载 skill.py
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                f"skill_{name}", skill_dir / "skill.py"
            )
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            # 加载专属工具
            tools_path = skill_dir / "tools.py"
            tools = []
            if tools_path.exists():
                tools = self._load_tools(tools_path, name)
            
            skill = Skill(meta=meta, prompt=prompt, module=module, tools=tools)
            self._skills[name] = skill
            return skill
        return None
    
    def _parse_markdown(self, content: str) -> tuple:
        """解析 Markdown 格式的 Skill 文件"""
        # 分割 YAML front matter 和正文
        parts = content.split("---", 2)
        if len(parts) >= 3:
            meta = yaml.safe_load(parts[1])
            prompt = parts[2].strip()
        else:
            meta = {}
            prompt = content
        return SkillMeta(**meta), prompt
    
    def match_by_keyword(self, user_input: str) -> List[str]:
        """根据用户输入匹配 Skill（关键词匹配）"""
        matched = []
        for name, skill in self._skills.items():
            for keyword in skill.meta.trigger_keywords:
                if keyword.lower() in user_input.lower():
                    matched.append(name)
                    break
        return matched

# skills/creator.py
class SkillCreator:
    """Skill 自动创建器（从对话中学习）"""
    
    def __init__(self, llm, workspace_path: str):
        self.llm = llm
        self.skills_dir = Path(workspace_path) / "config" / "skills"
    
    async def create_from_conversation(self, conversation: List[dict], 
                                        skill_name: str) -> Optional[Skill]:
        """从成功的对话中自动提取并创建 Skill"""
        prompt = f"""基于以下成功完成的对话，创建一个可复用的 Skill。
对话内容：
{json.dumps(conversation, ensure_ascii=False, indent=2)}

请生成：
1. SKILL.md 内容（YAML front matter + 提示词）
2. skill.py 执行逻辑代码

Skill 命名：{skill_name}
"""
        response = await self.llm.ainvoke([{"role": "user", "content": prompt}])
        # 解析 LLM 输出，写入文件，返回 Skill 对象
        return self._save_skill(skill_name, response)
```

        return self._save_skill(skill_name, response)
```

### 4.6 记忆与RAG系统设计

记忆系统分为三个层次：短期记忆（会话上下文）、长期记忆（持久化存储）、RAG检索增强（向量语义搜索）。

**短期记忆（memory/short_term.py）：**
```python
from collections import deque
from typing import List, Dict, Optional

class ShortTermMemory:
    """短期记忆：滑动窗口管理最近N轮对话"""
    
    def __init__(self, max_turns: int = 20, max_tokens: int = 100000):
        self.max_turns = max_turns
        self.max_tokens = max_tokens
        self._messages: deque = deque(maxlen=max_turns * 2)
        self._token_count = 0
    
    def add(self, role: str, content: str, token_estimate: int = None):
        self._messages.append({"role": role, "content": content})
        self._token_count += token_estimate or len(content) // 4
        self._trim()
    
    def _trim(self):
        """当超过Token限制时，自动裁剪最早的消息"""
        while self._token_count > self.max_tokens and len(self._messages) > 2:
            removed = self._messages.popleft()
            self._token_count -= len(removed["content"]) // 4
    
    def get_messages(self) -> List[dict]:
        return list(self._messages)
    
    def clear(self):
        self._messages.clear()
        self._token_count = 0
```

**长期记忆（memory/long_term.py）：**
```python
import sqlite3
from datetime import datetime, timezone
from typing import List, Optional

class LongTermMemory:
    """长期记忆：SQLite持久化用户偏好和历史决策"""
    
    def __init__(self, db_path: str = "~/.ils4gas/data/memory.db"):
        self.db_path = db_path
        self._init_db()
    
    def _init_db(self):
        with self._get_conn() as conn:
            conn.executescript("""
                CREATE TABLE IF NOT EXISTS memories (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    user_id TEXT NOT NULL,
                    content TEXT NOT NULL,
                    category TEXT DEFAULT 'general',
                    importance INTEGER DEFAULT 0,
                    created_at TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS conversation_summaries (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    session_id TEXT NOT NULL,
                    user_id TEXT NOT NULL,
                    summary TEXT NOT NULL,
                    created_at TEXT NOT NULL
                );
                CREATE INDEX IF NOT EXISTS idx_mem_user ON memories(user_id, category);
            """)
    
    def _get_conn(self):
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn
    
    def save(self, user_id: str, content: str, category: str = "general", 
             importance: int = 0):
        now = datetime.now(timezone.utc).isoformat()
        with self._get_conn() as conn:
            conn.execute(
                "INSERT INTO memories VALUES (NULL, ?, ?, ?, ?, ?)",
                (user_id, content, category, importance, now)
            )
    
    def search(self, user_id: str, query: str = "", limit: int = 5) -> List[dict]:
        sql = """SELECT * FROM memories WHERE user_id = ? 
                 AND content LIKE ? ORDER BY importance DESC LIMIT ?"""
        with self._get_conn() as conn:
            rows = conn.execute(sql, (user_id, f"%{query}%", limit)).fetchall()
        return [dict(r) for r in rows]
    
    def save_summary(self, session_id: str, user_id: str, summary: str):
        now = datetime.now(timezone.utc).isoformat()
        with self._get_conn() as conn:
            conn.execute(
                "INSERT INTO conversation_summaries VALUES (NULL, ?, ?, ?, ?)",
                (session_id, user_id, summary, now)
            )
```

**RAG检索增强（memory/rag.py）：**
```python
import chromadb
from chromadb.config import Settings
from typing import List, Dict

class RAGManager:
    """RAG检索增强生成管理器"""
    
    def __init__(self, persist_dir: str = "~/.ils4gas/data/vectors"):
        self.client = chromadb.PersistentClient(
            path=persist_dir,
            settings=Settings(anonymized_telemetry=False)
        )
        self._collections: Dict[str, any] = {}
    
    def get_or_create_collection(self, name: str):
        if name not in self._collections:
            self._collections[name] = self.client.get_or_create_collection(
                name=name, metadata={"hnsw:space": "cosine"}
            )
        return self._collections[name]
    
    async def index_messages(self, session_id: str, messages: List[dict], embedding_fn):
        """将历史消息向量化并存储"""
        collection = self.get_or_create_collection(f"session_{session_id}")
        for i, msg in enumerate(messages):
            embedding = await embedding_fn(msg["content"])
            collection.add(
                embeddings=[embedding],
                documents=[msg["content"]],
                metadatas=[{"role": msg["role"], "idx": i}],
                ids=[f"{session_id}_{i}"]
            )
    
    async def search(self, session_id: str, query: str, 
                     embedding_fn, limit: int = 5) -> List[dict]:
        """语义搜索历史消息"""
        collection = self.get_or_create_collection(f"session_{session_id}")
        query_embedding = await embedding_fn(query)
        results = collection.query(query_embeddings=[query_embedding], n_results=limit)
        return [
            {"content": doc, "metadata": meta}
            for doc, meta in zip(results["documents"][0], results["metadatas"][0])
        ]
    
    def index_documents(self, collection_name: str, documents: List[dict], embedding_fn):
        """索引外部文档（知识库）"""
        collection = self.get_or_create_collection(collection_name)
        for i, doc in enumerate(documents):
            embedding = embedding_fn(doc["content"])
            collection.add(
                embeddings=[embedding],
                documents=[doc["content"]],
                metadatas=[doc.get("metadata", {})],
                ids=[f"{collection_name}_{i}"]
            )

class ConversationSummarizer:
    """对话摘要器 - 自动压缩超长对话"""
    
    def __init__(self, llm):
        self.llm = llm
    
    async def summarize(self, messages: List[dict], max_length: int = 500) -> str:
        """将长对话压缩为摘要"""
        if len(messages) < 5:
            return ""
        prompt = f"""请将以下对话压缩为一句话摘要，保留关键的用户需求和决策：
{json.dumps(messages[-10:], ensure_ascii=False)}

摘要："""
        return await self.llm.ainvoke([{"role": "user", "content": prompt}])
```

### 4.7 A2A 多智能体通信设计

A2A (Agent-to-Agent) 协议用于多智能体之间的任务分发与协作。

**A2A协议消息格式（a2a/protocol.py）：**
```python
from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List
from uuid import uuid4

@dataclass
class A2AMessage:
    """A2A标准消息格式"""
    message_id: str = field(default_factory=lambda: uuid4().hex[:12])
    from_agent: str = ""
    to_agent: str = ""
    message_type: str = "task"  # task | result | query | notify
    content: str = ""
    context: Dict[str, Any] = field(default_factory=dict)
    priority: int = 0

@dataclass
class A2ATask:
    """A2A任务定义"""
    task_id: str = field(default_factory=lambda: uuid4().hex[:12])
    description: str = ""
    required_skills: List[str] = field(default_factory=list)
    required_tools: List[str] = field(default_factory=list)
    status: str = "pending"
    assigned_agent: Optional[str] = None
    result: Optional[str] = None
```

**多智能体编排器（a2a/orchestrator.py）：**
```python
import asyncio
import json
from typing import Dict, List, Optional

class AgentOrchestrator:
    """多智能体编排器：任务分解、分发、结果聚合"""
    
    def __init__(self):
        self._agents: Dict[str, any] = {}
        self._capabilities: Dict[str, List[str]] = {}
    
    def register_agent(self, name: str, agent, capabilities: List[str] = None):
        self._agents[name] = agent
        self._capabilities[name] = capabilities or []
    
    async def decompose_and_dispatch(self, task_description: str, llm) -> List[A2ATask]:
        """使用LLM分解任务并分发给合适的Agent"""
        # 1. LLM分解任务
        prompt = f"""将以下任务分解为独立的子任务，每个子任务描述需要的能力：
任务：{task_description}
输出格式（JSON数组）：[{{"description": "子任务描述", "required_skills": ["skill1"]}}]"""
        
        response = await llm.ainvoke([{"role": "user", "content": prompt}])
        subtasks_raw = json.loads(response)
        
        # 2. 为每个子任务匹配Agent
        tasks = []
        for st in subtasks_raw:
            task = A2ATask(**st)
            task.assigned_agent = self._match_agent(task)
            tasks.append(task)
        return tasks
    
    def _match_agent(self, task: A2ATask) -> Optional[str]:
        """根据能力匹配最合适的Agent"""
        best = None
        best_score = -1
        for name, caps in self._capabilities.items():
            score = sum(1 for s in task.required_skills if s in caps)
            if score > best_score:
                best_score = score
                best = name
        return best if best_score > 0 else None
    
    async def execute_tasks(self, tasks: List[A2ATask]) -> Dict[str, str]:
        """并发执行已分发的任务"""
        async def run_one(task):
            agent = self._agents.get(task.assigned_agent)
            if not agent:
                return task.task_id, None
            result = await agent.run(task.description)
            return task.task_id, result
        
        results = await asyncio.gather(*[run_one(t) for t in tasks])
        return dict(results)
```

### 4.8 Agent核心设计

**Agent基类（含中断/恢复支持）：**

```python
# core/agent.py
from abc import ABC, abstractmethod
from typing import Iterator, Optional, List, Dict, Any
from .message import Message
from .llm import LLMProvider
from .events import EventBus, AgentEvent
from tools.registry import ToolRegistry
from memory.short_term import ShortTermMemory

class Agent(ABC):
    """Agent基类 - 支持中断恢复、事件总线、工具调用"""
    
    def __init__(self, name: str, llm: LLMProvider, 
                 system_prompt: Optional[str] = None,
                 tool_registry: Optional[ToolRegistry] = None,
                 memory: Optional[ShortTermMemory] = None):
        self.name = name
        self.llm = llm
        self.system_prompt = system_prompt
        self.tool_registry = tool_registry
        self.memory = memory or ShortTermMemory()
        self._history: List[Message] = []
        self._event_bus = EventBus()
        self._cancelled = False
        self._paused = False
    
    @abstractmethod
    def run(self, input_text: str, **kwargs) -> str:
        pass
    
    @abstractmethod
    def stream_run(self, input_text: str, **kwargs) -> Iterator[str]:
        pass
    
    def cancel(self):
        """中断当前执行"""
        self._cancelled = True
        self._event_bus.emit(AgentEvent.CANCELLED)
    
    def pause(self):
        """暂停当前执行（保存检查点）"""
        self._paused = True
        checkpoint = self._save_checkpoint()
        self._event_bus.emit(AgentEvent.PAUSED, checkpoint=checkpoint)
        return checkpoint
    
    def resume(self, checkpoint: Dict[str, Any] = None):
        """从检查点恢复执行"""
        if checkpoint:
            self._load_checkpoint(checkpoint)
        self._paused = False
        self._event_bus.emit(AgentEvent.RESUMED)
    
    def _save_checkpoint(self) -> Dict[str, Any]:
        """保存当前执行状态"""
        return {
            "history": [msg.to_dict() for msg in self._history],
            "memory_messages": self.memory.get_messages(),
            "step": getattr(self, '_current_step', 0)
        }
    
    def _load_checkpoint(self, checkpoint: Dict[str, Any]):
        """从检查点恢复状态"""
        self._history = [Message.from_dict(m) for m in checkpoint["history"]]
        for msg in checkpoint.get("memory_messages", []):
            self.memory.add(msg["role"], msg["content"])

    def on_event(self, event: AgentEvent, handler):
        """注册事件处理器"""
        self._event_bus.on(event, handler)

    def add_message(self, role: str, content: str, tool_calls: list = None):
        msg = Message(role=role, content=content, tool_calls=tool_calls)
        self._history.append(msg)
        self.memory.add(role, content)
    
    def get_history(self) -> List[Message]:
        return self._history
```

**SimpleAgent实现：**

```python
# agents/simple_agent.py
from core.agent import Agent
from typing import Optional, Iterator

class SimpleAgent(Agent):
    """简单对话Agent - 无工具调用的基础对话"""
    
    def __init__(self, name: str, llm, system_prompt: Optional[str] = None, 
                 tool_registry=None):
        super().__init__(name, llm, system_prompt, tool_registry)
    
    def run(self, input_text: str, **kwargs) -> str:
        messages = []
        if self.system_prompt:
            messages.append({"role": "system", "content": self.system_prompt})
        for msg in self._history:
            messages.append(msg.to_openai_format())
        messages.append({"role": "user", "content": input_text})
        
        response = self.llm.invoke(messages, **kwargs)
        
        self.add_message("user", input_text)
        self.add_message("assistant", response)
        return response
    
    def stream_run(self, input_text: str, **kwargs) -> Iterator[str]:
        messages = []
        if self.system_prompt:
            messages.append({"role": "system", "content": self.system_prompt})
        for msg in self._history:
            messages.append(msg.to_openai_format())
        messages.append({"role": "user", "content": input_text})
        
        full_response = ""
        for chunk in self.llm.stream_invoke(messages, **kwargs):
            if self._cancelled:
                break
            full_response += chunk
            yield chunk
        
        self.add_message("user", input_text)
        self.add_message("assistant", full_response)

# core/events.py
from enum import Enum
from typing import Callable, Dict, List

class AgentEvent(Enum):
    CANCELLED = "cancelled"
    PAUSED = "paused"
    RESUMED = "resumed"
    TOOL_CALL_START = "tool_call_start"
    TOOL_CALL_END = "tool_call_end"
    STEP_START = "step_start"
    STEP_END = "step_end"
    ERROR = "error"

class EventBus:
    """Agent内部事件总线"""
    def __init__(self):
        self._handlers: Dict[AgentEvent, List[Callable]] = {}
    
    def on(self, event: AgentEvent, handler: Callable):
        if event not in self._handlers:
            self._handlers[event] = []
        self._handlers[event].append(handler)
    
    def emit(self, event: AgentEvent, **kwargs):
        for handler in self._handlers.get(event, []):
            handler(**kwargs)
```

**ReActAgent实现（推理-行动循环）：**

```python
# agents/react_agent.py
import re
from typing import Optional, Iterator
from core.agent import Agent, AgentEvent

REACT_PROMPT = """你是一个具备推理和行动能力的AI助手。

## 可用工具
{tools}

## 工作流程
请严格按照以下格式回应，每次只能执行一个步骤：

Thought: 你的思考过程
Action: 工具名[参数] 或 Finish[最终答案]

## 当前任务
**Question:** {question}

## 执行历史
{history}

现在开始你的推理和行动："""

class ReActAgent(Agent):
    """ReAct模式Agent - 交替进行推理(Thought)和行动(Action)"""
    
    def __init__(self, name: str, llm, system_prompt: Optional[str] = None,
                 tool_registry=None, max_steps: int = 10, 
                 custom_prompt: Optional[str] = None):
        super().__init__(name, llm, system_prompt, tool_registry)
        self.max_steps = max_steps
        self.prompt_template = custom_prompt or REACT_PROMPT
        self._current_step = 0
        self._step_history: List[str] = []
    
    def run(self, input_text: str, **kwargs) -> str:
        self._current_step = 0
        self._step_history = []
        
        while self._current_step < self.max_steps:
            if self._cancelled:
                return "执行已被中断"
            
            self._current_step += 1
            self._event_bus.emit(AgentEvent.STEP_START, step=self._current_step)
            
            # 构建Prompt
            tools_desc = self.tool_registry.get_tools_description() if self.tool_registry else ""
            history_str = "\n".join(self._step_history)
            prompt = self.prompt_template.format(
                tools=tools_desc, question=input_text, history=history_str
            )
            
            # 调用LLM
            messages = [{"role": "user", "content": prompt}]
            response = self.llm.invoke(messages, **kwargs)
            
            # 解析 Thought 和 Action
            thought, action = self._parse_output(response)
            self._step_history.append(f"Thought: {thought}")
            
            # 检查是否完成
            if action and action.startswith("Finish["):
                final_answer = action[7:-1]  # 去除 Finish[ 和 ]
                self.add_message("user", input_text)
                self.add_message("assistant", final_answer)
                return final_answer
            
            # 执行工具调用
            if action and self.tool_registry:
                tool_name, tool_input = self._parse_action(action)
                self._event_bus.emit(AgentEvent.TOOL_CALL_START, 
                                     tool_name=tool_name, args=tool_input)
                observation = self.tool_registry.execute_tool(tool_name, tool_input)
                self._event_bus.emit(AgentEvent.TOOL_CALL_END,
                                     tool_name=tool_name, result=observation)
                self._step_history.append(f"Action: {action}")
                self._step_history.append(f"Observation: {observation}")
            
            self._event_bus.emit(AgentEvent.STEP_END, step=self._current_step)
        
        return "抱歉，我无法在限定步数内完成这个任务。"
    
    def _parse_output(self, text: str) -> tuple:
        """解析LLM输出中的Thought和Action"""
        thought_match = re.search(r"Thought:\s*(.+?)(?:\n|$)", text, re.DOTALL)
        action_match = re.search(r"Action:\s*(.+?)(?:\n|$)", text, re.DOTALL)
        thought = thought_match.group(1).strip() if thought_match else ""
        action = action_match.group(1).strip() if action_match else ""
        return thought, action
    
    def _parse_action(self, action: str) -> tuple:
        """解析 Action: tool_name[params]"""
        match = re.match(r"(\w+)\[(.*)\]", action.strip())
        if match:
            return match.group(1), match.group(2)
        return action, ""

    def stream_run(self, input_text: str, **kwargs) -> Iterator[str]:
        """流式运行（ReAct模式下逐步产出）"""
        result = self.run(input_text, **kwargs)
        yield result
```

**ReflectionAgent实现（自我反思优化）：**

```python
# agents/reflection_agent.py
from typing import Optional, Iterator, List

class ReflectionAgent(Agent):
    """Reflection模式Agent - 生成后自我反思并优化"""
    
    def __init__(self, name: str, llm, system_prompt: Optional[str] = None,
                 tool_registry=None, max_reflections: int = 3):
        super().__init__(name, llm, system_prompt, tool_registry)
        self.max_reflections = max_reflections
    
    def run(self, input_text: str, **kwargs) -> str:
        # 第一步：生成初始回答
        initial_response = self._generate(input_text, **kwargs)
        
        # 多轮自我反思优化
        current_response = initial_response
        for i in range(self.max_reflections):
            if self._cancelled:
                break
            
            # 反思当前回答
            critique = self._reflect(input_text, current_response, **kwargs)
            
            # 如果有改进空间，重新生成
            if "需要改进" in critique or "可以优化" in critique:
                current_response = self._refine(
                    input_text, current_response, critique, **kwargs
                )
            else:
                break  # 无需进一步改进
        
        self.add_message("user", input_text)
        self.add_message("assistant", current_response)
        return current_response
    
    def _generate(self, input_text: str, **kwargs) -> str:
        messages = [{"role": "user", "content": input_text}]
        if self.system_prompt:
            messages.insert(0, {"role": "system", "content": self.system_prompt})
        return self.llm.invoke(messages, **kwargs)
    
    def _reflect(self, question: str, answer: str, **kwargs) -> str:
        """反思回答的质量"""
        prompt = f"""请严格审查以下问答，找出需要改进的地方：
问题：{question}
回答：{answer}

审查要点：
1. 回答是否准确完整？
2. 是否有逻辑漏洞？
3. 表达是否清晰？
4. 是否需要补充细节？

如果回答已经足够好，回复"无需改进"。
否则，指出具体需要改进的地方。"""
        return self.llm.invoke([{"role": "user", "content": prompt}], **kwargs)
    
    def _refine(self, question: str, answer: str, critique: str, **kwargs) -> str:
        """根据反思改进回答"""
        prompt = f"""根据以下反馈改进你的回答：
问题：{question}
原回答：{answer}
反馈：{critique}

请给出改进后的回答："""
        return self.llm.invoke([{"role": "user", "content": prompt}], **kwargs)
    
    def stream_run(self, input_text: str, **kwargs) -> Iterator[str]:
        result = self.run(input_text, **kwargs)
        yield result
```

### 4.9 异步工具执行器设计

```python
# tools/async_executor.py
import asyncio
from typing import Any, Dict, Optional, Set
from dataclasses import dataclass, field
from enum import Enum

class ToolPermission(Enum):
    READ_ONLY = "read_only"     # 只读工具（文件读取、搜索）
    WRITE = "write"            # 写入工具（文件编辑、创建）
    EXECUTE = "execute"        # 执行工具（代码运行、系统命令）
    NETWORK = "network"        # 网络工具（HTTP请求）

@dataclass
class ToolExecutionConfig:
    timeout_seconds: float = 30.0
    max_retries: int = 1
    retry_delay: float = 1.0
    max_concurrent: int = 5
    permission: ToolPermission = ToolPermission.READ_ONLY

class AsyncToolExecutor:
    """异步工具执行器 - 超时控制、重试、并发管理、权限分级"""
    
    def __init__(self, max_concurrent: int = 5):
        self._semaphore = asyncio.Semaphore(max_concurrent)
        self._tool_configs: Dict[str, ToolExecutionConfig] = {}
        self._running_tasks: Dict[str, asyncio.Task] = {}
    
    def register_tool_config(self, tool_name: str, config: ToolExecutionConfig):
        """注册工具的执行配置"""
        self._tool_configs[tool_name] = config
    
    def check_permission(self, tool_name: str, required: ToolPermission) -> bool:
        """检查工具是否有足够权限"""
        config = self._tool_configs.get(tool_name)
        if not config:
            return True  # 无配置的工具默认允许
        permissions_order = {
            ToolPermission.READ_ONLY: 0,
            ToolPermission.WRITE: 1,
            ToolPermission.NETWORK: 2,
            ToolPermission.EXECUTE: 3,
        }
        return permissions_order[config.permission] >= permissions_order[required]
    
    async def execute(self, tool_name: str, tool_fn, args: dict, 
                      session_id: str = None) -> Any:
        """执行工具（带超时、重试、并发控制）"""
        config = self._tool_configs.get(tool_name, ToolExecutionConfig())
        task_id = f"{session_id}:{tool_name}" if session_id else tool_name
        
        async with self._semaphore:  # 并发控制
            last_error = None
            for attempt in range(config.max_retries + 1):
                try:
                    result = await asyncio.wait_for(
                        self._call_tool(tool_fn, args),
                        timeout=config.timeout_seconds
                    )
                    return result
                except asyncio.TimeoutError:
                    last_error = f"工具 '{tool_name}' 执行超时 ({config.timeout_seconds}s)"
                except Exception as e:
                    last_error = str(e)
                    if attempt < config.max_retries:
                        await asyncio.sleep(config.retry_delay)
            
            raise RuntimeError(last_error)
    
    async def _call_tool(self, tool_fn, args: dict) -> Any:
        """实际调用工具函数（支持同步/异步）"""
        result = tool_fn(**args)
        if asyncio.iscoroutine(result):
            return await result
        return result
    
    def cancel(self, tool_name: str, session_id: str = None):
        """取消正在执行的工具"""
        task_id = f"{session_id}:{tool_name}" if session_id else tool_name
        if task_id in self._running_tasks:
            self._running_tasks[task_id].cancel()
```

### 4.10 Agent评估体系设计

框架内置评估能力，支持对Agent的行为进行量化评估。

```python
# evaluation/metrics.py
from dataclasses import dataclass, field
from typing import List, Optional
import time

@dataclass
class AgentMetrics:
    """Agent运行指标"""
    task_id: str = ""
    task_description: str = ""
    success: bool = False
    steps_taken: int = 0
    tools_called: int = 0
    unique_tools_used: List[str] = field(default_factory=list)
    total_time_seconds: float = 0.0
    llm_calls: int = 0
    total_tokens: int = 0
    final_answer: str = ""
    errors: List[str] = field(default_factory=list)

class AgentEvaluator:
    """Agent评估器"""
    
    def __init__(self, llm):
        self.llm = llm
    
    def evaluate(self, metrics: AgentMetrics) -> dict:
        """计算综合评估分数"""
        scores = {
            "task_completion": 1.0 if metrics.success else 0.0,
            "efficiency": self._evaluate_efficiency(metrics),
            "tool_usage": self._evaluate_tool_usage(metrics),
            "error_rate": len(metrics.errors) / max(metrics.llm_calls, 1)
        }
        scores["overall"] = sum(scores.values()) / len(scores)
        return scores
    
    def _evaluate_efficiency(self, m: AgentMetrics) -> float:
        """评估效率（步数/工具调用是否合理）"""
        if m.steps_taken == 0:
            return 0.0
        # 理想情况下步数应与工具调用数接近
        efficiency = m.tools_called / m.steps_taken if m.steps_taken > 0 else 0
        return min(efficiency, 1.0)
    
    def _evaluate_tool_usage(self, m: AgentMetrics) -> float:
        """评估工具使用（是否重复调用同一工具）"""
        if m.tools_called == 0:
            return 0.5
        # 工具使用的多样性
        diversity = len(m.unique_tools_used) / max(m.tools_called, 1)
        return diversity
    
    async def evaluate_answer_quality(self, question: str, answer: str, 
                                       expected: Optional[str] = None) -> dict:
        """使用LLM评估回答质量"""
        prompt = f"""评估以下问答的质量（1-10分）：
问题：{question}
回答：{answer}
{f'预期答案：{expected}' if expected else ''}

评估维度：
1. 准确性（回答是否正确）
2. 完整性（是否覆盖所有要点）
3. 清晰度（表达是否清晰）
4. 相关性（是否切题）

请以JSON格式输出：{{"accuracy": N, "completeness": N, "clarity": N, "relevance": N}}"""
        response = await self.llm.ainvoke([{"role": "user", "content": prompt}])
        import json
        return json.loads(response)
```

### 4.11 LLM服务设计

```python
# services/llm_service.py
import json
import os
import re
from pathlib import Path
from typing import List, Dict, Optional, Callable
from core.llm import LLMProvider, OpenAICompatibleProvider

class LLMService:
    """LLM管理服务 - 基于Provider抽象，支持动态切换与Agent缓存失效"""
    
    def __init__(self, config_path: str = "~/.ils4gas/config.json",
                 on_model_switch: Optional[Callable] = None):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.current_model = self.config.get("currentModel")
        self._provider_cache: Dict[str, LLMProvider] = {}
        self._on_model_switch = on_model_switch  # 模型切换回调（用于刷新Agent缓存）
    
    def _load_config(self) -> dict:
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = json.load(f)
        return self._resolve_env_vars(config)
    
    def _resolve_env_vars(self, config: dict) -> dict:
        config_str = json.dumps(config)
        pattern = r'\$\{ENV:([^}]+)\}'
        for var_name in re.findall(pattern, config_str):
            var_value = os.getenv(var_name, '')
            config_str = config_str.replace(f"${{ENV:{var_name}}}", var_value)
        return json.loads(config_str)
    
    def get_provider(self, model_id: Optional[str] = None) -> LLMProvider:
        """获取Provider实例（带缓存）"""
        model_id = model_id or self.current_model
        provider_name, model_name = model_id.split('/', 1)
        
        cache_key = f"{provider_name}/{model_name}"
        if cache_key in self._provider_cache:
            return self._provider_cache[cache_key]
        
        provider_config = self.config["providers"][provider_name]
        model_config = provider_config["models"][model_name]
        options = provider_config["options"]
        
        provider_type = provider_config.get("type", "openai-compatible")
        
        if provider_type in ("openai", "openai-compatible"):
            provider = OpenAICompatibleProvider(
                api_key=options.get("apiKey", ""),
                base_url=options.get("baseURL"),
                model_name=model_name,
                context_limit=model_config["limit"].get("context", 128000)
            )
        else:
            # 默认回退到OpenAI兼容
            provider = OpenAICompatibleProvider(
                api_key=options.get("apiKey", ""),
                base_url=options.get("baseURL"),
                model_name=model_name
            )
        
        self._provider_cache[cache_key] = provider
        return provider
    
    def switch_model(self, model_id: str):
        """切换模型并通知缓存刷新"""
        old_model = self.current_model
        self.current_model = model_id
        if self._on_model_switch and old_model != model_id:
            self._on_model_switch(model_id)
    
    def get_model_info(self, model_id: Optional[str] = None) -> dict:
        model_id = model_id or self.current_model
        provider_name, model_name = model_id.split('/', 1)
        provider = self.config["providers"][provider_name]
        model = provider["models"][model_name]
        return {
            "id": model_id,
            "name": model["name"],
            "provider": provider_name,
            "limit": model.get("limit", {}),
            "supports_tools": model.get("supports_tools", True),
            "supports_reasoning": model.get("supports_reasoning", False)
        }
    
    def list_models(self) -> List[Dict]:
        models = []
        for provider_name, provider in self.config["providers"].items():
            for model_name, model in provider["models"].items():
                models.append({
                    "id": f"{provider_name}/{model_name}",
                    "name": model["name"],
                    "provider": provider_name,
                    "limit": model.get("limit", {})
                })
        return models
    
    def invalidate_cache(self):
        """刷新Provider缓存（用于API Key变更等场景）"""
        self._provider_cache.clear()
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
POST   /api/v1/skills/load              # 加载指定Skill
POST   /api/v1/skills/unload            # 卸载Skill
POST   /api/v1/skills/create            # 从对话创建新Skill
GET    /api/v1/skills/{name}            # 获取Skill详情
PUT    /api/v1/skills/{name}            # 更新Skill
DELETE /api/v1/skills/{name}            # 删除Skill

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

6. **Skills系统**
   - `skills/registry.py` - Skill注册与发现
   - `skills/loader.py` - Skill加载器
   - `skills/executor.py` - Skill执行引擎
   - `skills/creator.py` - Skill自动创建器

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
