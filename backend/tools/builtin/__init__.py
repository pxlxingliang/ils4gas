from backend.tools.registry import ToolRegistry


def register_builtin_tools() -> ToolRegistry:
    from backend.tools.builtin.read import create_read_tool
    from backend.tools.builtin.glob import create_glob_tool
    from backend.tools.builtin.grep import create_grep_tool
    from backend.tools.builtin.bash import create_bash_tool
    from backend.tools.builtin.write import create_write_tool
    from backend.tools.builtin.edit import create_edit_tool
    from backend.tools.builtin.webfetch import create_webfetch_tool
    from backend.tools.builtin.websearch import create_websearch_tool
    from backend.tools.builtin.skill import create_skill_tool

    registry = ToolRegistry()

    for factory in [
        create_read_tool,
        create_glob_tool,
        create_grep_tool,
        create_bash_tool,
        create_write_tool,
        create_edit_tool,
        create_webfetch_tool,
        create_websearch_tool,
    ]:
        registry.register(factory())

    skill_tool = create_skill_tool()
    if skill_tool:
        registry.register(skill_tool)

    return registry
