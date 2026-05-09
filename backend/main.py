import getpass
import os
import socket
import sys
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

from backend.api import deps
from backend.core.context import WorkspaceContext
from backend.services.llm_service import LLMService
from backend.services.session_service import SessionService
from backend.services.mcp_service import MCPService
from backend.tools.mcp_adapter import MCPServerConfig
from backend.tools.loader import load_tools_from_mcp_modules


@asynccontextmanager
async def lifespan(app: FastAPI):
    from backend.core.env import EnvManager
    EnvManager.init()

    # Initialize workspace files on first launch
    WorkspaceContext.init_workspace()

    # Load built-in MCP tools for web chat
    tool_registry = load_tools_from_mcp_modules()
    print(f"  tools: {len(tool_registry)} built-in tools loaded")

    deps.llm_service = LLMService(tool_registry=tool_registry)
    deps.session_service = SessionService()

    # MCP servers from config (external servers to connect to)
    mcp_entries = deps.llm_service.config.get("mcp_servers", [])
    deps.mcp_service = MCPService.from_config_entries(mcp_entries)
    await deps.mcp_service.connect_all()

    current = deps.llm_service.get_current_model()
    external_tools = len(deps.mcp_service.list_all_tools())
    port = os.environ.get("ILS4GAS_WEB_PORT", "8789")
    print(f"  model: {current['id']}")
    print(f"  mcp  : {len(mcp_entries)} external servers, {external_tools} external tools")
    print()
    print(f"  Local:   http://localhost:{port}")
    hostname = socket.gethostname()
    print(f"  Remote:  ssh -N -L {port}:localhost:{port} {getpass.getuser()}@{hostname}")
    print("           (set up SSH key first: ssh-copy-id user@server)")
    print(f"           Then open http://localhost:{port} in your local browser")

    yield

    await deps.mcp_service.disconnect_all()


app = FastAPI(title="ILS4GAS", version="0.1.0", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# API routes
from backend.api.routes import chat, llm, session, mcp, skill, memory, ws

app.include_router(chat.router)
app.include_router(llm.router)
app.include_router(session.router)
app.include_router(mcp.router)
app.include_router(skill.router)
app.include_router(memory.router)
app.include_router(ws.router)

# Frontend static files
static_dir = Path(__file__).parent / "static"
if static_dir.exists():
    app.mount("/assets", StaticFiles(directory=static_dir / "assets"), name="assets")

    @app.get("/{full_path:path}", response_class=HTMLResponse)
    async def spa(full_path: str = ""):
        index = static_dir / "index.html"
        if index.exists():
            return HTMLResponse(
                content=index.read_text(encoding="utf-8"),
                headers={"Cache-Control": "no-cache"},
            )
        return HTMLResponse("<h1>Frontend not built. Run: cd frontend && npm install && npm run build</h1>")


def run_web(port: int = None, host: str = None):
    import uvicorn
    from backend.core.port import find_free_port

    p = port or int(os.getenv("ILS4GAS_WEB_PORT", "8789"))
    h = host or os.getenv("ILS4GAS_WEB_HOST", "0.0.0.0")
    actual = find_free_port(p, h)
    if actual != p:
        print(f"  port {p} in use, using {actual}")
    os.environ["ILS4GAS_WEB_PORT"] = str(actual)
    uvicorn.run("backend.main:app", host=h, port=actual, reload=False)


if __name__ == "__main__":
    run_web()
