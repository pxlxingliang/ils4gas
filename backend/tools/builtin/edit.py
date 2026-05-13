import os

from backend.tools.registry import ToolInfo

DESCRIPTION = """\
Performs exact string replacements in an existing file.

Usage:
- When editing text, ensure you preserve the exact indentation (tabs/spaces) as it appears BEFORE.
- ALWAYS prefer editing existing files in the codebase. NEVER write new files unless explicitly required.
- Only use emojis if the user explicitly requests it. Avoid adding emojis to files unless asked.
- The edit will FAIL if `oldString` is not found in the file with an error "oldString not found in content".
- The edit will FAIL if `oldString` is found multiple times in the file with an error "Found multiple matches for oldString. Provide more surrounding lines in oldString to identify the correct match." Either provide a larger string with more surrounding context to make it unique or use `replaceAll` to change every instance of `oldString`.
- Use `replaceAll` for replacing and renaming strings across the file. This parameter is useful if you want to rename a variable for instance."""

PARAMETERS = {
    "type": "object",
    "properties": {
        "filePath": {
            "type": "string",
            "description": "The absolute path to the file to modify",
        },
        "oldString": {
            "type": "string",
            "description": "The text to replace",
        },
        "newString": {
            "type": "string",
            "description": "The text to replace it with (must be different from oldString)",
        },
        "replaceAll": {
            "type": "boolean",
            "description": "Replace all occurrences of oldString (default false)",
            "default": False,
        },
    },
    "required": ["filePath", "oldString", "newString"],
}


def _normalize(string: str) -> str:
    return string.replace("\r\n", "\n")


# ── Replacer strategies ──────────────────────────────────

def _simple_replacer(content: str, find: str):
    if find in content:
        yield find


def _line_trimmed_replacer(content: str, find: str):
    original_lines = content.split("\n")
    search_lines = find.split("\n")
    if search_lines and search_lines[-1] == "":
        search_lines.pop()
    for i in range(len(original_lines) - len(search_lines) + 1):
        match = all(
            original_lines[i + j].strip() == search_lines[j].strip()
            for j in range(len(search_lines))
        )
        if match:
            start = sum(len(original_lines[k]) + 1 for k in range(i))
            end = start + sum(len(search_lines[k]) + 1 for k in range(len(search_lines))) - 1
            if end > start:
                yield content[start:end]


def _escape_normalized_replacer(content: str, find: str):
    def unescape(s: str) -> str:
        out = []
        i = 0
        while i < len(s):
            if s[i] == "\\" and i + 1 < len(s):
                nxt = s[i + 1]
                if nxt == "n":
                    out.append("\n")
                elif nxt == "t":
                    out.append("\t")
                elif nxt == "r":
                    out.append("\r")
                elif nxt in ("'", '"', "`", "\\", "$"):
                    out.append(nxt)
                else:
                    out.append(s[i])
                    out.append(nxt)
                i += 2
            else:
                out.append(s[i])
                i += 1
        return "".join(out)

    unescaped = unescape(find)
    if unescaped == find:
        return
    if unescaped in content:
        yield unescaped

    content_lines = content.split("\n")
    find_lines = unescaped.split("\n")
    for i in range(len(content_lines) - len(find_lines) + 1):
        block = "\n".join(content_lines[i:i + len(find_lines)])
        if unescape(block) == unescaped:
            yield block


def _multi_occurrence_replacer(content: str, find: str):
    idx = 0
    while True:
        idx = content.find(find, idx)
        if idx == -1:
            break
        yield find
        idx += len(find)


# ── Replace engine ───────────────────────────────────────

def _replace(content: str, oldString: str, newString: str, replaceAll: bool = False) -> str:
    if oldString == newString:
        raise ValueError("oldString and newString must be different")

    content = _normalize(content)
    oldString = _normalize(oldString)
    newString = _normalize(newString)

    replacers = [
        _simple_replacer,
        _line_trimmed_replacer,
        _escape_normalized_replacer,
        _multi_occurrence_replacer,
    ]

    not_found = True
    for replacer in replacers:
        for search in replacer(content, oldString):
            idx = content.find(search)
            if idx == -1:
                continue
            not_found = False
            if replaceAll:
                return content.replace(search, newString)
            last_idx = content.rfind(search)
            if idx != last_idx:
                continue
            return content[:idx] + newString + content[idx + len(search):]

    if not_found:
        raise ValueError("oldString not found in content")
    raise ValueError(
        "Found multiple matches for oldString. "
        "Provide more surrounding lines in oldString to identify the correct match."
    )


# ── Execute ──────────────────────────────────────────────

def _execute(
    filePath: str,
    oldString: str,
    newString: str,
    replaceAll: bool = False,
) -> str:
    filepath = filePath
    if not os.path.isabs(filepath):
        filepath = os.path.abspath(filepath)

    if not os.path.isfile(filepath):
        return f"File not found: {filepath}"

    if os.path.isdir(filepath):
        return f"Path is a directory, not a file: {filepath}"

    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()
    except (OSError, IOError) as e:
        return f"Error reading file: {e}"

    if oldString == newString:
        return "No changes: oldString and newString are identical"

    try:
        new_content = _replace(content, oldString, newString, replaceAll)
    except ValueError as e:
        return str(e)

    try:
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(new_content)
    except (OSError, IOError) as e:
        return f"Error writing file: {e}"

    return "Edit applied successfully."


def create_edit_tool() -> ToolInfo:
    return ToolInfo(
        name="edit",
        description=DESCRIPTION,
        parameters=PARAMETERS,
        fn=_execute,
    )
