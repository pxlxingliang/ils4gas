import os
import re
import subprocess
from pathlib import Path

from backend.tools.registry import ToolInfo

DESCRIPTION = """\
Fast content search tool that works with any codebase size
Searches file contents using regular expressions
Supports full regex syntax (eg. "log.*Error", "function\\s+\\w+")
Filter files by pattern with the include parameter (eg. "*.js", "*.{ts,tsx}")
Returns file paths and line numbers with at least one match sorted by modification time
Use this tool when you need to find files containing specific patterns
If you need to identify/count the number of matches within files, use the Bash tool \
with 'rg' (ripgrep) directly. Do NOT use 'grep'.
When you are doing an open-ended search that may require multiple rounds of \
globbing and grepping, use the Task tool instead"""

PARAMETERS = {
    "type": "object",
    "properties": {
        "pattern": {
            "type": "string",
            "description": "The regex pattern to search for in file contents",
        },
        "path": {
            "type": "string",
            "description": "The directory to search in. Defaults to the current working directory.",
        },
        "include": {
            "type": "string",
            "description": "File pattern to include in the search (e.g. \"*.js\", \"*.{ts,tsx}\")",
        },
    },
    "required": ["pattern"],
}

MAX_LINE_LENGTH = 2000
MAX_RESULTS = 100


def _find_ripgrep() -> str:
    try:
        result = subprocess.run(
            ["which", "rg"], capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        pass
    return ""


def _parse_ripgrep_output(output: str, search_dir: str) -> list:
    matches = []
    for line in output.strip().split("\n"):
        if not line:
            continue
        parts = line.split("|", 2)
        if len(parts) < 3:
            continue
        file_path, line_num_str, line_text = parts[0], parts[1], parts[2]
        if not os.path.isabs(file_path):
            file_path = os.path.join(search_dir, file_path)
        try:
            line_num = int(line_num_str)
        except ValueError:
            continue
        try:
            mtime = os.path.getmtime(file_path)
        except OSError:
            continue
        matches.append((file_path, mtime, line_num, line_text))
    return matches


def _search_python(pattern: str, search_dir: str, include: str = "") -> list:
    matches = []
    compiled = re.compile(pattern)
    include_pattern = None

    if include:
        import fnmatch
        include_pattern = re.compile(
            fnmatch.translate(include)
        )

    for root, dirs, files in os.walk(search_dir):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        for filename in files:
            if filename.startswith("."):
                continue
            if include_pattern and not include_pattern.match(filename):
                continue
            filepath = os.path.join(root, filename)
            try:
                mtime = os.path.getmtime(filepath)
                with open(filepath, "r", encoding="utf-8", errors="replace") as f:
                    for line_num, line in enumerate(f, start=1):
                        if compiled.search(line):
                            line = line.rstrip("\n")
                            matches.append((filepath, mtime, line_num, line))
            except (OSError, IOError, UnicodeDecodeError):
                continue

    return matches


def _execute(pattern: str, path: str = "", include: str = "") -> str:
    search_dir = os.path.abspath(path) if path else os.getcwd()
    if not os.path.isdir(search_dir):
        return f"Directory not found: {search_dir}"

    matches = []

    rg_path = _find_ripgrep()
    if rg_path:
        args = [
            rg_path, "-nH", "--hidden", "--follow", "--no-messages",
            "--field-match-separator=|", "--regexp", pattern,
        ]
        if include:
            args.extend(["--glob", include])
        args.append(search_dir)

        try:
            proc = subprocess.run(args, capture_output=True, text=True, timeout=60)
            if proc.returncode == 1 or (proc.returncode == 2 and not proc.stdout.strip()):
                return "No files found"
            if proc.returncode not in (0, 2):
                pass
            else:
                matches = _parse_ripgrep_output(proc.stdout, search_dir)
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            pass

    if not matches:
        matches = _search_python(pattern, search_dir, include)

    if not matches:
        return "No files found"

    matches.sort(key=lambda x: x[1], reverse=True)

    truncated = len(matches) > MAX_RESULTS
    final_matches = matches[:MAX_RESULTS]

    output_lines = [f"Found {len(final_matches)} matches"]

    current_file = ""
    for filepath, _mtime, line_num, line_text in final_matches:
        if current_file != filepath:
            if current_file:
                output_lines.append("")
            current_file = filepath
            output_lines.append(f"{filepath}:")
        display_text = line_text[:MAX_LINE_LENGTH]
        if len(line_text) > MAX_LINE_LENGTH:
            display_text += "..."
        output_lines.append(f"  Line {line_num}: {display_text}")

    if truncated:
        output_lines.append("")
        output_lines.append("(Results are truncated. Consider using a more specific path or pattern.)")

    return "\n".join(output_lines)


def create_grep_tool() -> ToolInfo:
    return ToolInfo(
        name="grep",
        description=DESCRIPTION,
        parameters=PARAMETERS,
        fn=_execute,
    )
