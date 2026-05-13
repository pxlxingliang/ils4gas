import os
import subprocess
import signal
import threading

from backend.tools.registry import ToolInfo

DESCRIPTION = """\
Executes a given bash command in a persistent shell session with optional timeout, ensuring proper handling and security measures.

All commands run in the current working directory by default. Use the `workdir` parameter if you need to run a command in a different directory. AVOID using `cd <directory> && <command>` patterns - use `workdir` instead.

IMPORTANT: This tool is for terminal operations like git, npm, docker, etc. DO NOT use it for file operations (reading, writing, editing, searching, finding files) - use the specialized tools for this instead.

Before executing the command, please follow these steps:

1. Directory Verification:
   - If the command will create new directories or files, first use `ls` to verify the parent directory exists and is the correct location
   - For example, before running "mkdir foo/bar", first use `ls foo` to check that "foo" exists and is the intended parent directory

2. Command Execution:
   - Always quote file paths that contain spaces with double quotes (e.g., rm "path with spaces/file.txt")
   - Examples of proper quoting:
     - mkdir "/Users/name/My Documents" (correct)
     - mkdir /Users/name/My Documents (incorrect - will fail)
     - python "/path/with spaces/script.py" (correct)
     - python /path/with spaces/script.py (incorrect - will fail)
   - After ensuring proper quoting, execute the command.
   - Capture the output of the command.

Usage notes:
  - The command argument is required.
  - You can specify an optional timeout in milliseconds. If not specified, commands will time out after 120000ms (2 minutes).
  - It is very helpful if you write a clear, concise description of what this command does in 5-10 words.
  - If the output exceeds 2000 lines or 51200 bytes, it will be truncated and the full output will be written to a file. You can use Read with offset/limit to read specific sections or Grep to search the full content. Do NOT use `head`, `tail`, or other truncation commands to limit output; the full output will already be captured to a file for more precise searching.

  - Avoid using Bash with the `find`, `grep`, `cat`, `head`, `tail`, `sed`, `awk`, or `echo` commands, unless explicitly instructed or when these commands are truly necessary for the task. Instead, always prefer using the dedicated tools for these commands:
    - File search: Use Glob (NOT find or ls)
    - Content search: Use Grep (NOT grep or rg)
    - Read files: Use Read (NOT cat/head/tail)
    - Edit files: Use Edit (NOT sed/awk)
    - Write files: Use Write (NOT echo >/cat <<EOF)
    - Communication: Output text directly (NOT echo/printf)
  - When issuing multiple commands:
    - If the commands are independent and can run in parallel, make multiple bash tool calls in a single message.
    - If the commands depend on each other and must run sequentially, use a single Bash call with '&&' to chain them together.
    - Use ';' only when you need to run commands sequentially but don't care if earlier commands fail
    - DO NOT use newlines to separate commands (newlines are ok in quoted strings)
  - AVOID using `cd <directory> && <command>`. Use the `workdir` parameter to change directories instead."""

PARAMETERS = {
    "type": "object",
    "properties": {
        "command": {
            "type": "string",
            "description": "The command to execute",
        },
        "timeout": {
            "type": "integer",
            "description": "Optional timeout in milliseconds. Max 600000 (10 minutes). Default 120000 (2 minutes).",
        },
        "workdir": {
            "type": "string",
            "description": "The working directory to run the command in. Defaults to the current working directory. Use this instead of 'cd' commands.",
        },
        "description": {
            "type": "string",
            "description": "Clear, concise description of what this command does in 5-10 words.",
        },
    },
    "required": ["command", "description"],
}

DEFAULT_TIMEOUT = 120_000
MAX_TIMEOUT = 600_000
MAX_OUTPUT_LINES = 2000
MAX_OUTPUT_BYTES = 51200


def _execute(
    command: str,
    description: str,
    timeout: int = DEFAULT_TIMEOUT,
    workdir: str = "",
) -> str:
    cwd = workdir or os.getcwd()
    if not os.path.isdir(cwd):
        return f"Working directory not found: {cwd}"

    timeout = min(max(1, timeout), MAX_TIMEOUT)
    timeout_sec = timeout / 1000.0

    try:
        proc = subprocess.Popen(
            command,
            shell=True,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            preexec_fn=os.setsid if os.name != "nt" else None,
        )
    except Exception as e:
        return f"Error executing command: {e}"

    output_chunks: list[bytes] = []
    killed = False

    def read_output():
        try:
            for chunk in iter(lambda: proc.stdout.read(4096), b""):  # type: ignore
                output_chunks.append(chunk)
        except Exception:
            pass

    reader = threading.Thread(target=read_output)
    reader.start()

    try:
        proc.wait(timeout=timeout_sec)
        reader.join(timeout=5)
    except subprocess.TimeoutExpired:
        killed = True
        try:
            if os.name != "nt":
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            else:
                proc.kill()
            proc.wait(timeout=5)
        except Exception:
            pass
        reader.join(timeout=5)

    raw = b"".join(output_chunks)
    output = raw.decode("utf-8", errors="replace")

    lines = output.split("\n")
    byte_len = len(output.encode("utf-8"))

    truncated = len(lines) > MAX_OUTPUT_LINES or byte_len > MAX_OUTPUT_BYTES
    if truncated:
        output = "\n".join(lines[:MAX_OUTPUT_LINES])
        if len(output.encode("utf-8")) > MAX_OUTPUT_BYTES:
            output = output.encode("utf-8")[:MAX_OUTPUT_BYTES].decode("utf-8", errors="replace")

        output_file = os.path.join(cwd, ".opencode_bash_output.txt")
        try:
            with open(output_file, "w", encoding="utf-8") as f:
                full = raw.decode("utf-8", errors="replace")
                f.write(full)
            output += f"\n\n(Output truncated. Full output saved to {output_file})"
        except Exception:
            output += "\n\n(Output truncated)"

    if killed:
        output += f"\n\nCommand timed out after {timeout}ms"

    return (output or "(no output)").strip()


def create_bash_tool() -> ToolInfo:
    return ToolInfo(
        name="bash",
        description=DESCRIPTION,
        parameters=PARAMETERS,
        fn=_execute,
    )
