from textual.app import App
from backend.tui.screen import ChatScreen


class ILS4GASApp(App):
    CSS = """
    Screen {
        background: #ffffff;
        color: #1a1a2e;
    }
    Header {
        background: #f0f0f5;
        color: #1a1a2e;
    }
    #history {
        height: 1fr;
    }
    #bottom-bar {
        height: auto;
    }
    Input {
        background: #f7f7f8;
        color: #1a1a2e;
        border: solid #e5e7eb;
    }
    Input:focus {
        border: solid #6366f1;
    }
    StatusBar {
        background: #f0f0f5;
        color: #6b7280;
    }
    #thinking {
        color: #6366f1;
    }
    ListView {
        background: #ffffff;
        color: #1a1a2e;
    }
    ListView > ListItem {
        background: #ffffff;
        color: #1a1a2e;
        padding: 1;
    }
    ListView > ListItem:hover {
        background: #eef2ff;
    }
    ListView > ListItem.--highlight {
        background: #e0e7ff;
    }
    #command-dropdown {
        background: #ffffff;
        border: solid #6366f1;
        height: auto;
        max-height: 8;
    }
    #command-dropdown > ListItem {
        background: #ffffff;
        color: #1a1a2e;
        padding: 0 1;
    }
    #command-dropdown > ListItem:hover {
        background: #eef2ff;
    }
    #command-dropdown > ListItem.--highlight {
        background: #e0e7ff;
        color: #6366f1;
    }
    """

    BINDINGS = [
        ("ctrl+q", "quit", "Quit"),
    ]

    def on_mount(self):
        self.push_screen(ChatScreen())
