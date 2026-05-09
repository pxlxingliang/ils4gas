from textual.screen import Screen
from textual.containers import Vertical, Horizontal
from textual.widgets import Static, Input, Label, Header, ListView, ListItem
from textual.binding import Binding
from textual import on

from backend.tui.widgets.history import MessageHistory
from backend.tui.widgets.status_bar import StatusBar

from backend.agents.react_agent import ReActAgent
from backend.services.llm_service import LLMService
from backend.services.session_service import SessionService
from backend.tools.loader import load_tools_from_mcp_modules
from backend.core.context import WorkspaceContext
from backend.core.events import AgentEventType
from backend.skills.registry import SkillRegistry
from backend.services.title_service import generate_title

COMMAND_DESCRIPTIONS = {
    "/exit": "Exit ILS4GAS",
    "/quit": "Exit ILS4GAS",
    "/model": "Switch LLM model",
    "/help": "Show available commands",
    "/clear": "Start a new session",
    "/new": "Start a new session",
}


class ModelSelectScreen(Screen):
    def __init__(self, models: list[dict], callback):
        super().__init__()
        self._models = models
        self._callback = callback
        self._id_map: dict[str, str] = {}

    def compose(self):
        items = []
        for m in self._models:
            safe_id = m["id"].replace("/", "_").replace(".", "_")
            self._id_map[safe_id] = m["id"]
            items.append(ListItem(Label(f"{m['name']}  ({m['provider']})"), id=safe_id))
        self._list = ListView(*items)
        yield Header()
        yield Vertical(
            Static("Select a model (Enter to confirm, Esc to cancel):"),
            self._list,
        )

    def on_list_view_selected(self, event: ListView.Selected):
        if event.item and event.item.id:
            real_id = self._id_map.get(event.item.id, event.item.id)
            self._callback(real_id)
            self.app.pop_screen()

    def on_screen_resume(self):
        self._list.focus()


class CommandDropdown(ListView):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.can_focus = True
        self.display = False

    async def show_matches(self, prefix: str):
        prefix = prefix.strip().lower()
        await self.clear()
        if not prefix.startswith("/"):
            self.display = False
            return
        matches = [
            (cmd, desc)
            for cmd, desc in COMMAND_DESCRIPTIONS.items()
            if cmd.startswith(prefix)
        ]
        if matches:
            items = [
                ListItem(Label(f"{cmd}  ({desc})"), id="cmd_" + cmd.lstrip("/"))
                for cmd, desc in matches
            ]
            await self.extend(items)
            self.display = True
        else:
            self.display = False

    async def hide_dropdown(self):
        self.display = False
        await self.clear()


class ChatScreen(Screen):
    BINDINGS = [
        Binding("ctrl+n", "new_session", "New"),
        Binding("escape", "dismiss_dropdown", "Dismiss"),
    ]

    def __init__(self):
        super().__init__()
        self._llm = LLMService()
        self._sess = SessionService()
        self._tool_registry = load_tools_from_mcp_modules()
        self._llm.tool_registry = self._tool_registry
        self._session_id: str = ""

    def _make_agent(self):
        system_prompt = WorkspaceContext().build_system_prompt()
        registry = SkillRegistry()
        registry.load_all()
        skill_prompts = registry.get_active_prompts()
        if skill_prompts:
            system_prompt = system_prompt + "\n\n" + skill_prompts
        return ReActAgent(
            llm_service=self._llm,
            tool_registry=self._tool_registry,
            system_prompt=system_prompt,
        )

    def compose(self):
        sess = self._sess.create_session(
            model_provider=self._llm.get_current_model().get("provider", ""),
            model_name=self._llm.get_current_model().get("id", ""),
        )
        self._session_id = sess["id"]
        model_name = self._llm.get_current_model().get("name", "?")
        yield Header()
        yield MessageHistory(id="history")
        yield Vertical(
            Input(placeholder="Type a message... (Enter to send, / for commands)", id="chat-input"),
            CommandDropdown(id="command-dropdown"),
            Label("", id="thinking"),
            StatusBar(f"model: {model_name} | / for help | {len(self._tool_registry)} tools | Shift+Mouse to select | Ctrl+Q quit"),
            id="bottom-bar",
        )

    def on_mount(self):
        self.query_one("#chat-input", Input).focus()

    # ── command dropdown ──────────────────────────────────

    @on(Input.Changed)
    async def on_input_changed(self, event: Input.Changed):
        value = event.value
        dropdown = self.query_one("#command-dropdown", CommandDropdown)
        if value.startswith("/"):
            await dropdown.show_matches(value)
        else:
            await dropdown.hide_dropdown()

    async def action_dismiss_dropdown(self):
        dropdown = self.query_one("#command-dropdown", CommandDropdown)
        if dropdown.display:
            await dropdown.hide_dropdown()
        self.query_one("#chat-input", Input).focus()

    @on(ListView.Selected)
    async def on_command_selected(self, event: ListView.Selected):
        if event.control.id == "command-dropdown" and event.item:
            safe_id = event.item.id
            cmd = "/" + safe_id[4:]
            input_widget = self.query_one("#chat-input", Input)
            input_widget.value = cmd + " "
            input_widget.cursor_position = len(cmd) + 1
            await self.query_one("#command-dropdown", CommandDropdown).hide_dropdown()
            input_widget.focus()
            event.stop()

    # ── chat input ────────────────────────────────────────

    async def on_input_submitted(self, event: Input.Submitted):
        content = event.value.strip()
        if not content:
            return

        input_widget = self.query_one("#chat-input", Input)
        input_widget.clear()
        await self.query_one("#command-dropdown", CommandDropdown).hide_dropdown()

        if content.startswith("/"):
            self._handle_command(content)
            input_widget.focus()
            return

        history = self.query_one("#history", MessageHistory)
        history.add_message("user", content)
        self.run_worker(self._stream_response(content), exclusive=False)

    # ── commands ──────────────────────────────────────────

    def _handle_command(self, cmd: str):
        history = self.query_one("#history", MessageHistory)
        parts = cmd.split(maxsplit=1)
        cmd = parts[0].lower()

        if cmd in ("/exit", "/quit"):
            history.add_message("system", "Goodbye!")
            self.app.exit()

        elif cmd == "/help":
            help_text = (
                "Commands:\n"
                "  /exit, /quit    Exit ILS4GAS\n"
                "  /model          Switch LLM model\n"
                "  /clear, /new    Start a new session\n"
                "  /help           Show this help\n"
                "\n"
                "Tips:\n"
                "  Shift+Mouse drag to select and copy text"
            )
            history.add_message("system", help_text)

        elif cmd in ("/clear", "/new"):
            self.action_new_session()

        elif cmd == "/model":
            self._show_model_selector()

        else:
            history.add_message("system", f"Unknown command: {cmd}.  Type /help for commands.")

    def _show_model_selector(self):
        models = self._llm.list_models()

        def on_model_selected(model_id: str):
            self._llm.switch_model(model_id)
            info = self._llm.get_current_model()
            history = self.query_one("#history", MessageHistory)
            history.add_message("system", f"Switched to {info['name']} ({info['id']})")
            self._update_status()

        self.app.push_screen(ModelSelectScreen(models, on_model_selected))

    def _update_status(self):
        model_name = self._llm.get_current_model().get("name", "?")
        status = self.query_one(StatusBar)
        status.update(
            f"model: {model_name} | / for commands | {len(self._tool_registry)} tools"
            f" | Shift+Mouse to select text | Ctrl+Q to quit"
        )

    async def _auto_title(self, user_content: str):
        session = self._sess.get_session(self._session_id)
        if not session or session.get("title") != "New Chat":
            return
        try:
            title = await generate_title(user_content, self._llm)
            if title:
                self._sess.update_session(self._session_id, title=title)
        except Exception:
            pass

    # ── streaming via ReActAgent ──────────────────────────

    async def _stream_response(self, user_content: str):
        history = self.query_one("#history", MessageHistory)
        thinking = self.query_one("#thinking", Label)

        thinking.update("thinking...")
        self._sess.add_message(self._session_id, "user", user_content)

        # Build history, skipping tool-only messages (API rejects incomplete sequences)
        messages = []
        for msg in self._sess.get_messages(self._session_id)[:-1]:
            role = msg["role"]
            if role == "tool":
                continue
            if role == "assistant" and msg.get("tool_calls") and not msg.get("content"):
                continue
            entry = {"role": role, "content": msg["content"]}
            if msg.get("tool_calls"):
                entry["tool_calls"] = msg["tool_calls"]
            if msg.get("tool_call_id"):
                entry["tool_call_id"] = msg["tool_call_id"]
            messages.append(entry)
        messages.append({"role": "user", "content": user_content})

        try:
            agent = self._make_agent()
            accumulated = ""
            tool_calls: list[dict] = []

            async for event in agent.stream_run(messages):
                if event.type == AgentEventType.CONTENT_CHUNK:
                    accumulated += event.data["text"]
                    history.update_last(accumulated)

                elif event.type == AgentEventType.TOOL_CALL_START:
                    tool_name = event.data.get("tool_name", "?")
                    tc_id = event.data.get("tool_call_id", "")
                    tool_calls.append({
                        "id": tc_id,
                        "type": "function",
                        "function": {
                            "name": tool_name,
                            "arguments": event.data.get("args", ""),
                        },
                    })
                    thinking.update(f"calling {tool_name}...")

                elif event.type == AgentEventType.TOOL_CALL_END:
                    if tool_calls:
                        tool_calls[-1]["_result"] = event.data.get("result", "")
                    thinking.update("thinking...")

                elif event.type == AgentEventType.ERROR:
                    history.add_message("system", f"Error: {event.data.get('message', 'unknown')}")

            if accumulated:
                # Save assistant message with proper OpenAI tool_calls format
                saved_calls = [
                    {k: v for k, v in tc.items() if k != "_result"}
                    for tc in tool_calls
                ] if tool_calls else None
                self._sess.add_message(
                    self._session_id, "assistant", accumulated,
                    tool_calls=saved_calls,
                )
                # Save individual tool result messages for complete history
                for tc in tool_calls:
                    if "_result" in tc:
                        self._sess.add_message(
                            self._session_id, "tool", tc["_result"],
                            tool_call_id=tc["id"],
                        )
                # Auto-generate title on first user message
                self.run_worker(
                    self._auto_title(user_content),
                    exclusive=False,
                )

        except Exception as e:
            history.add_message("system", f"Error: {e}")
        finally:
            thinking.update("")

    def action_new_session(self):
        history = self.query_one("#history", MessageHistory)
        history.add_message("system", "Starting new session...")
        sess = self._sess.create_session(
            model_provider=self._llm.get_current_model().get("provider", ""),
            model_name=self._llm.get_current_model().get("id", ""),
        )
        self._session_id = sess["id"]
        history.clear()
        self._update_status()
