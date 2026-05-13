import os
from pathlib import Path

from backend.tools.registry import ToolInfo

DESCRIPTION = """\
Reads a file from the local filesystem. You can access any file directly by using this tool.
Assume this tool is able to read all files on the machine. If the User provides a path to a file \
assume that path is valid. It is okay to read a file that does not exist; an error will be returned.

Usage:
- The file_path parameter must be an absolute path, not a relative path
- By default, it reads up to 2000 lines starting from the beginning of the file
- You can optionally specify an offset and limit (especially handy for long files), \
but it's recommended to read the whole file by not providing these parameters
- Any lines longer than 2000 characters will be truncated
- Results are returned using cat -n format, with line numbers starting at 1
- You have the capability to call multiple tools in a single response. \
It is always better to speculatively read multiple files as a batch that are potentially useful.
- You can read image files and PDFs using this tool."""

PARAMETERS = {
    "type": "object",
    "properties": {
        "file_path": {
            "type": "string",
            "description": "The absolute path to the file to read",
        },
        "offset": {
            "type": "integer",
            "description": "The line number to start reading from (0-based). 0 means start at line 1.",
            "default": 0,
        },
        "limit": {
            "type": "integer",
            "description": "The number of lines to read. Defaults to 2000.",
            "default": 2000,
        },
    },
    "required": ["file_path"],
}

MAX_LINE_LENGTH = 2000
MAX_BYTES = 50 * 1024

BINARY_EXTENSIONS = {
    ".zip", ".tar", ".gz", ".bz2", ".xz", ".7z",
    ".exe", ".dll", ".so", ".dylib", ".class", ".jar", ".war",
    ".doc", ".docx", ".xls", ".xlsx", ".ppt", ".pptx",
    ".odt", ".ods", ".odp",
    ".bin", ".dat", ".obj", ".o", ".a", ".lib", ".wasm",
    ".pyc", ".pyo",
    ".png", ".jpg", ".jpeg", ".gif", ".bmp", ".ico", ".webp",
    ".mp3", ".mp4", ".avi", ".mov", ".mkv",
    ".ttf", ".otf", ".woff", ".woff2",
    ".pdf", ".epub",
}


def _is_binary(filepath: str) -> bool:
    ext = Path(filepath).suffix.lower()
    if ext in BINARY_EXTENSIONS:
        return True

    try:
        with open(filepath, "rb") as f:
            chunk = f.read(4096)
    except (OSError, IOError):
        return True

    if not chunk:
        return False

    non_printable = 0
    total = len(chunk)
    for b in chunk:
        if b == 0:
            return True
        if b < 9 or (13 < b < 32):
            non_printable += 1
    return non_printable / total > 0.3


def _fuzzy_suggest(filepath: str) -> str:
    directory = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    if not os.path.isdir(directory):
        return ""
    try:
        entries = os.listdir(directory)
    except OSError:
        return ""
    lower_base = basename.lower()
    suggestions = []
    for entry in entries:
        lower_entry = entry.lower()
        if lower_entry in lower_base or lower_base in lower_entry:
            suggestions.append(os.path.join(directory, entry))
    if suggestions:
        return "\n\nDid you mean one of these?\n" + "\n".join(suggestions[:3])
    return ""


def _execute(file_path: str, offset: int = 0, limit: int = 2000) -> str:
    filepath = file_path
    if not os.path.isabs(filepath):
        filepath = os.path.abspath(filepath)

    if not os.path.isfile(filepath):
        suggestion = _fuzzy_suggest(filepath)
        if suggestion:
            return f"File not found: {file_path}{suggestion}"
        return f"File not found: {file_path}"

    if _is_binary(filepath):
        ext = Path(filepath).suffix.lower()
        if ext in (".png", ".jpg", ".jpeg", ".gif", ".bmp", ".ico", ".webp"):
            return f"Image file: {file_path} (use Read to view as attachment)"
        if ext == ".pdf":
            return f"PDF file: {file_path} (use Read to view as attachment)"
        return f"Cannot read binary file: {file_path}"

    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            lines = f.read().split("\n")
    except (OSError, IOError) as e:
        return f"Error reading file: {e}"

    total_lines = len(lines)
    offset = max(0, offset)
    limit = max(1, min(limit, 2000))

    raw_lines = []
    byte_count = 0
    truncated_by_bytes = False

    for i in range(offset, min(total_lines, offset + limit)):
        line = lines[i]
        if len(line) > MAX_LINE_LENGTH:
            line = line[:MAX_LINE_LENGTH] + "..."
        size = len(line.encode("utf-8")) + (1 if raw_lines else 0)
        if byte_count + size > MAX_BYTES:
            truncated_by_bytes = True
            break
        raw_lines.append(line)
        byte_count += size

    output_lines = []
    for i, line in enumerate(raw_lines):
        output_lines.append(f"{offset + i + 1:05d}| {line}")

    if not output_lines and offset >= total_lines:
        return f"Offset {offset} is beyond end of file ({total_lines} lines total)"

    output = "\n".join(output_lines)

    last_read = offset + len(raw_lines)
    if truncated_by_bytes:
        output += f"\n\n(Output truncated at {MAX_BYTES} bytes. Use 'offset' parameter to read beyond line {last_read})"
    elif total_lines > last_read:
        output += f"\n\n(File has more lines. Use 'offset' parameter to read beyond line {last_read})"
    else:
        output += f"\n\n(End of file - total {total_lines} lines)"

    return output


def create_read_tool() -> ToolInfo:
    return ToolInfo(
        name="read",
        description=DESCRIPTION,
        parameters=PARAMETERS,
        fn=_execute,
    )
