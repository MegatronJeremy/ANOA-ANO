"""
NERV / MAGI-themed Textual TUI for the pipeline menu.

Convenience layer ONLY -- the canonical, graded entry point stays the headless
`run_pipeline.py --stage/--all` CLI. This app is pure presentation: it renders
stage status and, when you launch a stage, it *suspends* (hands the terminal
back) and calls the exact same run callback the CLI uses, so all the existing
rich logging / progress output is unchanged. Nothing about the analysis lives
here.

Theme: Evangelion NERV console -- amber on void black, alert red, terminal
green for "ready", with a MAGI header and a DNA motif.
"""
from __future__ import annotations

from typing import Callable

from rich.text import Text
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal
from textual.widgets import Button, Checkbox, DataTable, Footer, Rule, Static

# --- NERV palette -----------------------------------------------------------
AMBER = "#F5A100"
ALERT = "#E63329"
GREEN = "#7CFC00"
VOID = "#0A0A0A"
PANEL = "#1A1206"

_BANNER = (
    f"[b {AMBER}]N E R V[/]   [{ALERT}]//[/]   [b {AMBER}]MAGI SYSTEM[/]\n"
    f"[dim {AMBER}]scRNA-seq · immune response to nanoplastics · particle-size study[/]\n"
    f"[{GREEN}]⌁[/] [dim]DNA  A–T  G≡C  C≡G  T–A  G≡C  A–T  C≡G[/] [{GREEN}]⌁[/]"
)


class PipelineTUI(App):
    """Full-screen menu; delegates all real work to `run_selection`."""

    CSS = f"""
    Screen {{
        background: {VOID};
        color: {AMBER};
    }}
    #banner {{
        padding: 1 2 0 2;
        color: {AMBER};
        text-align: center;
    }}
    #hazard {{
        color: {ALERT};
        margin: 0 2;
    }}
    DataTable {{
        background: {VOID};
        color: {AMBER};
        border: heavy {AMBER};
        margin: 1 2;
        height: auto;
    }}
    DataTable > .datatable--header {{
        background: {PANEL};
        color: {AMBER};
        text-style: bold;
    }}
    DataTable > .datatable--cursor {{
        background: {AMBER};
        color: {VOID};
        text-style: bold;
    }}
    #controls {{
        height: auto;
        padding: 0 2 1 2;
        align: left middle;
    }}
    Checkbox {{
        background: {VOID};
        color: {AMBER};
        border: none;
        margin: 0 3 0 0;
        width: auto;
    }}
    Button {{ margin: 0 1 0 0; width: auto; }}
    Button#run    {{ background: {AMBER}; color: {VOID}; text-style: bold; }}
    Button#runall {{ background: {ALERT}; color: {VOID}; text-style: bold; }}
    Footer {{ background: {PANEL}; color: {AMBER}; }}
    """

    TITLE = "NERV // MAGI"
    SUB_TITLE = "scRNA-seq nanoplastic pipeline"

    BINDINGS = [
        Binding("a", "run_all", "Run ALL"),
        Binding("s", "toggle_smoke", "Smoke"),
        Binding("d", "toggle_debug", "Debug"),
        Binding("r", "refresh", "Refresh"),
        Binding("q", "quit", "Quit"),
    ]

    def __init__(self, registry: dict, status_fn: Callable, run_selection: Callable):
        super().__init__()
        self.registry = registry
        self.status_fn = status_fn
        self.run_selection = run_selection
        self.stage_keys = list(registry.keys())

    def compose(self) -> ComposeResult:
        yield Static(_BANNER, id="banner")
        yield Rule(id="hazard")
        yield DataTable(id="stages", cursor_type="row", zebra_stripes=False)
        with Horizontal(id="controls"):
            yield Checkbox("SMOKE", id="smoke")
            yield Checkbox("DEBUG", id="debug")
            yield Button("RUN ▶", id="run")
            yield Button("RUN ALL ⏩", id="runall")
            yield Button("QUIT ✖", id="quit")
        yield Footer()

    def on_mount(self) -> None:
        table = self.query_one("#stages", DataTable)
        # Fixed widths so the status columns (the point of the grid) are never
        # squeezed off-screen; FUNCTION is truncated rather than allowed to grow.
        table.add_column("STAGE", width=30)
        table.add_column("FUNCTION", width=40)
        table.add_column("INPUT", width=16)
        table.add_column("FULL", width=6)
        table.add_column("SMOKE", width=7)
        self._populate()
        table.focus()

    # -- rendering -----------------------------------------------------------
    def _populate(self) -> None:
        table = self.query_one("#stages", DataTable)
        table.clear()
        done = f"[{GREEN}]●[/]"
        todo = "[grey37]—[/]"
        for spec in self.registry.values():
            st = self.status_fn(spec)
            if st["input_ready"]:
                inp = f"[{GREEN}]● ready[/]"
            else:
                need = st["input_desc"].replace(".h5ad", "")
                inp = f"[{ALERT}]○[/] [dim]{need}[/]"
            table.add_row(
                Text.from_markup(f"[b]{spec.label}[/b]"),
                Text.from_markup(f"[dim]{spec.description}[/dim]"),
                Text.from_markup(inp),
                Text.from_markup(done if st["full_done"] else todo),
                Text.from_markup(done if st["smoke_done"] else todo),
                key=spec.key,
            )

    # -- actions -------------------------------------------------------------
    def _run(self, selected: str) -> None:
        smoke = self.query_one("#smoke", Checkbox).value
        debug = self.query_one("#debug", Checkbox).value
        with self.suspend():
            self.run_selection(selected, smoke, debug, None)
            try:
                input("\n[ press Enter to return to NERV ] ")
            except EOFError:
                pass
        self._populate()

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        # Enter on a stage row runs that stage.
        self._run(event.row_key.value)

    def action_run_all(self) -> None:
        self._run("__all__")

    def action_toggle_smoke(self) -> None:
        cb = self.query_one("#smoke", Checkbox)
        cb.value = not cb.value

    def action_toggle_debug(self) -> None:
        cb = self.query_one("#debug", Checkbox)
        cb.value = not cb.value

    def action_refresh(self) -> None:
        self._populate()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "run":
            table = self.query_one("#stages", DataTable)
            row = table.cursor_row
            if row is not None and 0 <= row < len(self.stage_keys):
                self._run(self.stage_keys[row])
        elif event.button.id == "runall":
            self.action_run_all()
        elif event.button.id == "quit":
            self.exit()
