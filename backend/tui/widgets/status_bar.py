from textual.app import ComposeResult
from textual.widgets import Static


class StatusBar(Static):
    def __init__(self, text: str = ""):
        super().__init__(text)
