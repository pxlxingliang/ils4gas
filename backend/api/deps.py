from backend.services.llm_service import LLMService
from backend.services.session_service import SessionService
from backend.services.mcp_service import MCPService

llm_service: LLMService = None  # type: ignore
session_service: SessionService = None  # type: ignore
mcp_service: MCPService = None  # type: ignore


def get_llm_service() -> LLMService:
    global llm_service
    if llm_service is None:
        llm_service = LLMService()
    return llm_service


def get_session_service() -> SessionService:
    global session_service
    if session_service is None:
        session_service = SessionService()
    return session_service


def get_mcp_service() -> MCPService:
    global mcp_service
    if mcp_service is None:
        mcp_service = MCPService.from_config_entries([])
    return mcp_service
