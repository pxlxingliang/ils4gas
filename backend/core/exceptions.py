class AgentError(Exception):
    pass


class ToolNotFoundError(AgentError):
    pass


class LLMError(AgentError):
    pass
