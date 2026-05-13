import os
from pathlib import Path

from backend.tools.registry import ToolInfo

DESCRIPTION = """\
Fast file pattern matching tool that works with any codebase size
Supports glob patterns like "**/*.js" or "src/**/*.ts"
Returns matching file paths sorted by modification time
Use this tool when you need to find files by name patterns
When you are doing an open-ended search that may require multiple rounds of \
globbing and grepping, use the Task tool instead
You have the capability to call multiple tools in a single response. \
It is always better to speculatively perform multiple searches as a batch \
that are potentially useful."""

PARAMETERS = {
    "type": "object",
    "properties": {
        "pattern": {
            "type": "string",
            "description": "The glob pattern to match files against",
        },
        "path": {
            "type": "string",
            "description": "The directory to search in. If not specified, the current working directory will be used. IMPORTANT: Omit this field to use the default directory. DO NOT enter 'undefined' or 'null' - simply omit it for the default behavior. Must be a valid directory path if provided.",
        },
    },
    "required": ["pattern"],
}

MAX_RESULTS = 100


def _execute(pattern: str, path: str = "") -> str:
    search_dir = Path(path) if path else Path.cwd()
    if not search_dir.is_absolute():
        search_dir = search_dir.resolve()

    if not search_dir.is_dir():
        return f"Directory not found: {search_dir}"

    files = []
    try:
        for filepath in search_dir.rglob(pattern):
            if filepath.is_file():
                try:
                    mtime = filepath.stat().st_mtime
                except OSError:
                    mtime = 0
                files.append((str(filepath), mtime))
    except Exception as e:
        return f"Error searching files: {e}"

    files.sort(key=lambda x: x[1], reverse=True)

    truncated = len(files) > MAX_RESULTS
    results = files[:MAX_RESULTS]

    if not results:
        return "No files found"

    output = [p for p, _ in results]
    if truncated:
        output.append("")
        output.append("(Results are truncated. Consider using a more specific path or pattern.)")

    return "\n".join(output)


def create_glob_tool() -> ToolInfo:
    return ToolInfo(
        name="glob",
        description=DESCRIPTION,
        parameters=PARAMETERS,
        fn=_execute,
    )
