from rich.text import Text
from textual.widgets import Static
from textual.containers import VerticalScroll


class MessageHistory(VerticalScroll):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._msgs: list[tuple[str, str]] = []

    def add_message(self, role: str, content: str):
        self._msgs.append((role, content))
        self._render_all()

    def update_last(self, content: str):
        if self._msgs and self._msgs[-1][0] == "assistant":
            self._msgs[-1] = ("assistant", content)
        else:
            self._msgs.append(("assistant", content))
        self._render_all()

    def clear(self):
        self._msgs.clear()
        self.remove_children()

    def _render_all(self):
        self.remove_children()
        for role, content in self._msgs:
            if role == "user":
                prefix = "[bold blue]You:[/] "
            elif role == "system":
                prefix = "[bold grey50]●[/] "
            else:
                prefix = "[bold purple]AI:[/] "
            text = Text.from_markup(prefix + content)
            self.mount(Static(text))
        self.scroll_end(animate=False)
