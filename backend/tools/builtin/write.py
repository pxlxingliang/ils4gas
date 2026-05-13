import os
import difflib
from pathlib import Path

from backend.tools.registry import ToolInfo

DESCRIPTION = """\
Writes a file to the local filesystem.

Usage:
- This tool will overwrite the existing file if there is one at the provided path.
- If this is an existing file, you MUST use the Read tool first to read the file's contents.
- ALWAYS prefer editing existing files in the codebase. NEVER write new files unless explicitly required.
- NEVER proactively create documentation files (*.md) or README files. Only create documentation files if explicitly requested by the User.
- Only use emojis if the user explicitly requests it. Avoid writing emojis to files unless asked."""

PARAMETERS = {
    "type": "object",
    "properties": {
        "filePath": {
            "type": "string",
            "description": "The absolute path to the file to write (must be absolute, not relative)",
        },
        "content": {
            "type": "string",
            "description": "The content to write to the file",
        },
    },
    "required": ["filePath", "content"],
}


def _make_diff(filepath: str, old: str, new: str) -> str:
    diff = difflib.unified_diff(
        old.splitlines(keepends=True),
        new.splitlines(keepends=True),
        fromfile=filepath,
        tofile=filepath,
        lineterm="",
    )
    return "".join(diff)


def _execute(filePath: str, content: str) -> str:
    filepath = filePath
    if not os.path.isabs(filepath):
        filepath = os.path.abspath(filepath)

    parent = os.path.dirname(filepath)
    if parent and not os.path.isdir(parent):
        return f"Parent directory does not exist: {parent}"

    existed = os.path.isfile(filepath)
    old_content = ""
    if existed:
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                old_content = f.read()
        except (OSError, IOError) as e:
            return f"Error reading file: {e}"

    if old_content == content:
        return "No changes: new content is identical to existing file content"

    try:
        os.makedirs(parent, exist_ok=True)
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(content)
    except (OSError, IOError) as e:
        return f"Error writing file: {e}"

    diff = _make_diff(filepath, old_content, content)
    action = "Updated" if existed else "Created"
    return f"{action} file: {filepath}\n\nDiff:\n{diff}" if diff else f"{action} file: {filepath}"


def create_write_tool() -> ToolInfo:
    return ToolInfo(
        name="write",
        description=DESCRIPTION,
        parameters=PARAMETERS,
        fn=_execute,
    )
