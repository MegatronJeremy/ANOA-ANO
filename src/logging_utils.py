"""
Logging + debug-mode display helpers shared by every stage.

This module is a PRESENTATION layer only -- nothing in here computes or
changes any pipeline numbers. It controls how progress, structured numeric
dumps, and verification plots are rendered:

  - Plain log lines go through the standard `logging` module, rendered by
    `rich.logging.RichHandler`. Quiet runs stay at INFO, --debug switches to
    DEBUG -- same semantics as before, just nicer rendering.
  - Stage start/end banners are `rich` Panels, always visible (both modes),
    printed directly via the shared Console rather than through a log record.
  - Structured numeric dumps (filter summaries, percentile tables, per-sample
    counts) go through `log_table()` -- one reusable helper so stages 2-6
    get the same look for free.
  - `term_hist()` / `term_bar()` render compact terminal plots with
    `plotext`, gated behind --debug, used for at-a-glance verification only.
    They do not replace any matplotlib figures saved to disk.
  - `progress_spinner()` wraps slow steps (e.g. full, non-smoke-test runs)
    with a `rich` spinner, sharing the same Console as the logger so output
    doesn't interleave/garble.

Everything here degrades gracefully when stdout isn't a real terminal
(piped to a file, CI, etc.): `rich` automatically drops ANSI styling when
`Console.is_terminal` is False, and the `plotext` plots are skipped entirely
(replaced by a one-line debug note) rather than emitting unrenderable
escape-code soup into a file.
"""
import logging
import time
from contextlib import contextmanager

import numpy as np
from rich.console import Console
from rich.logging import RichHandler
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

_LOGGER_NAME = "nanoplastic"
_console: Console | None = None

# Compact, fixed plot size so terminal plots can't flood the screen
# regardless of how wide the actual terminal is.
_PLOT_WIDTH = 70
_PLOT_HEIGHT = 16


def get_console() -> Console:
    """Single shared Console instance -- used by the logger, progress bars,
    and table/panel output, so none of them interleave into garbled output."""
    global _console
    if _console is None:
        _console = Console()
    return _console


def setup_logging(debug: bool = False) -> logging.Logger:
    """Configure and return the project logger. Call once, at startup."""
    logger = logging.getLogger(_LOGGER_NAME)
    logger.handlers.clear()

    console = get_console()
    level = logging.DEBUG if debug else logging.INFO
    handler = RichHandler(
        console=console,
        show_time=debug,
        show_level=True,
        show_path=False,
        markup=False,
        rich_tracebacks=False,
    )
    logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False
    return logger


def get_logger() -> logging.Logger:
    return logging.getLogger(_LOGGER_NAME)


def _debug_enabled() -> bool:
    return get_logger().isEnabledFor(logging.DEBUG)


def silence_third_party_logger(name: str, debug: bool, debug_level=logging.INFO):
    """
    Some dependencies configure their own logging.Logger + handler at import
    time, completely independent of our setup_logging() (notably harmonypy,
    which attaches its own StreamHandler at DEBUG unconditionally). Without
    this, their internal messages would print in quiet runs too. Called once
    per stage that pulls in a chatty dependency.
    """
    logging.getLogger(name).setLevel(debug_level if debug else logging.WARNING)


@contextmanager
def log_stage_banner(name: str):
    """
    Wrap a pipeline stage: prints a start/end rich Panel with wall-clock
    time. Always visible (quiet and --debug runs alike) -- printed directly
    to the shared Console rather than gated by the logger's level, matching
    the previous INFO-level banner behaviour.
    """
    console = get_console()
    console.print(Panel(f"[bold]{name}[/bold]", title="STAGE START", border_style="cyan"))
    t0 = time.time()
    try:
        yield
    finally:
        dt = time.time() - t0
        console.print(
            Panel(f"[bold]{name}[/bold]\nelapsed: {dt:.1f}s", title="STAGE END", border_style="green")
        )


# Backwards-compatible alias (previous name used by run_pipeline.py).
stage_banner = log_stage_banner


@contextmanager
def progress_spinner(description: str):
    """
    Wrap a slow step with a rich spinner + elapsed time, sharing the same
    Console as the logger. Intended for non-smoke-test (full) runs where a
    step can take ~minutes; reusable for any later stage.

    Falls back to a single before/after log line when stdout isn't a real
    terminal, so piped/CI output stays plain text instead of getting a
    spinner's redraw escape codes.
    """
    console = get_console()
    log = get_logger()

    if not console.is_terminal:
        log.info(f"{description} ...")
        t0 = time.time()
        yield
        log.info(f"{description} done ({time.time() - t0:.1f}s)")
        return

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    )
    with progress:
        task = progress.add_task(description, total=None)
        yield
        progress.update(task, completed=1)


def require(condition: bool, message: str):
    """
    Loud assertion for stage boundaries. Raises with a clear, actionable
    message instead of letting a downstream stage fail cryptically.
    """
    if not condition:
        log = get_logger()
        log.error(f"GUARD FAILED: {message}")
        raise AssertionError(message)


def log_table(title: str, columns: list, rows: list, debug_only: bool = True):
    """
    Render a structured numeric dump (filter summaries, percentile tables,
    per-sample counts, ...) as a rich Table instead of a printed dict/line.
    One reusable entry point so stages 2-6 inherit the same look.

    `columns` is a list of header strings; `rows` is a list of row
    tuples/lists, each matching len(columns).
    """
    if debug_only and not _debug_enabled():
        return
    console = get_console()
    table = Table(title=title, show_header=True, header_style="bold cyan")
    for col in columns:
        table.add_column(str(col))
    for row in rows:
        table.add_row(*[str(c) for c in row])
    console.print(table)


def describe_adata(adata, label: str, debug_only: bool = True):
    """
    Log shape + structure of an AnnData object: n_obs/n_vars, obs columns,
    obsm keys, layers, HVG count if present, and per-sample cell counts (as
    a table). Used at every stage boundary.
    """
    log = get_logger()
    emit = log.debug if debug_only else log.info

    emit(f"[{label}] shape: {adata.n_obs} cells x {adata.n_vars} genes")
    emit(f"[{label}] obs columns: {list(adata.obs.columns)}")
    emit(f"[{label}] obsm keys: {list(adata.obsm.keys())}")
    emit(f"[{label}] layers: {list(adata.layers.keys())}")
    if "highly_variable" in adata.var.columns:
        n_hvg = int(adata.var["highly_variable"].sum())
        emit(f"[{label}] highly_variable genes flagged: {n_hvg}")
    if adata.obs.get("sample") is not None:
        try:
            counts = adata.obs["sample"].value_counts().to_dict()
            log_table(
                f"[{label}] per-sample cell counts",
                ["sample", "n_cells"],
                sorted(counts.items()),
                debug_only=debug_only,
            )
        except Exception:
            pass


def log_percentiles(values, name: str, percentiles=(1, 5, 25, 50, 75, 95, 99)):
    """Compute a percentile summary of a 1-D array-like, for QC distributions.
    Returns the dict so callers can also feed it into log_table / term_hist."""
    arr = np.asarray(values)
    pct = np.percentile(arr, percentiles)
    return dict(zip(percentiles, pct)), float(arr.min()), float(arr.max())


def log_percentile_table(metrics: dict, percentiles=(1, 5, 25, 50, 75, 95, 99), debug_only: bool = True):
    """
    Render percentile summaries for several QC metrics (e.g. n_genes_by_counts,
    total_counts, pct_counts_mt) as one combined rich Table.
    `metrics` maps metric name -> 1-D array-like of values.
    """
    columns = ["metric"] + [f"p{p}" for p in percentiles] + ["min", "max"]
    rows = []
    for name, values in metrics.items():
        pct, vmin, vmax = log_percentiles(values, name, percentiles)
        rows.append([name] + [f"{pct[p]:.1f}" for p in percentiles] + [f"{vmin:.1f}", f"{vmax:.1f}"])
    log_table("QC percentile summary", columns, rows, debug_only=debug_only)


# ---------------------------------------------------------------------------
# Terminal plots (plotext). Verification-only: never replace saved
# matplotlib figures, only give an at-a-glance look in the terminal itself.
# Gated behind --debug by callers; additionally guarded here against
# non-terminal stdout so piped/CI output never gets raw escape codes.
# ---------------------------------------------------------------------------

def _plot_allowed() -> bool:
    return _debug_enabled() and get_console().is_terminal


def term_hist(values, title: str, bins: int = 30, vlines=None, xlabel: str = ""):
    """
    Render a compact terminal histogram via plotext, with optional vertical
    lines marking thresholds (e.g. QC cutoffs) so they can be checked
    against the real distribution at a glance.

    `vlines` is an optional list of (x_value, color) tuples.
    Skips (with a one-line debug note) when not in --debug or not a real
    terminal, instead of emitting broken escape codes into a file/CI log.
    """
    log = get_logger()
    if not _plot_allowed():
        log.debug(f"[plot skipped: not a TTY or not --debug] histogram '{title}'")
        return

    import plotext as tplt

    tplt.clear_figure()
    tplt.theme("clear")
    tplt.plotsize(_PLOT_WIDTH, _PLOT_HEIGHT)
    tplt.hist(np.asarray(values), bins=bins)
    for x, color in (vlines or []):
        tplt.vline(x, color=color)
    tplt.title(title)
    if xlabel:
        tplt.xlabel(xlabel)
    get_console().print(tplt.build())
    tplt.clear_figure()


def term_bar(categories, series: dict, title: str, ylabel: str = ""):
    """
    Render a compact terminal bar chart via plotext. `series` maps a series
    label (e.g. "before", "after") to a list of values aligned with
    `categories`. A single-series chart is just `{"value": [...]}`.

    Skips (with a one-line debug note) when not in --debug or not a real
    terminal, same as term_hist.
    """
    log = get_logger()
    if not _plot_allowed():
        log.debug(f"[plot skipped: not a TTY or not --debug] bar chart '{title}'")
        return

    import plotext as tplt

    tplt.clear_figure()
    tplt.theme("clear")
    tplt.plotsize(_PLOT_WIDTH, _PLOT_HEIGHT)
    labels = list(categories)
    if len(series) == 1:
        ((_, values),) = series.items()
        tplt.bar(labels, list(values))
    else:
        tplt.multiple_bar(labels, [list(v) for v in series.values()], labels=list(series.keys()))
    tplt.title(title)
    if ylabel:
        tplt.ylabel(ylabel)
    get_console().print(tplt.build())
    tplt.clear_figure()


def term_scatter(x, y, groups=None, title: str = "", xlabel: str = "", ylabel: str = "",
                  max_points: int = 2000, seed: int = 0):
    """
    Render a compact terminal scatter plot via plotext, optionally split into
    coloured series by `groups` (e.g. sample/batch, for an at-a-glance look
    at batch separation before/after integration, or cluster identity).

    Subsamples to `max_points` total points (reproducibly, via `seed`) purely
    for terminal-plot legibility -- this never touches the actual data the
    pipeline operates on, only what gets drawn in this one plot.

    Skips (with a one-line debug note) when not in --debug or not a real
    terminal, same as term_hist/term_bar.
    """
    log = get_logger()
    if not _plot_allowed():
        log.debug(f"[plot skipped: not a TTY or not --debug] scatter '{title}'")
        return

    import plotext as tplt

    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)
    if n > max_points:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n, size=max_points, replace=False)
        x, y = x[idx], y[idx]
        groups = np.asarray(groups)[idx] if groups is not None else None

    tplt.clear_figure()
    tplt.theme("clear")
    tplt.plotsize(_PLOT_WIDTH, _PLOT_HEIGHT)
    if groups is None:
        tplt.scatter(x, y, marker="dot")
    else:
        groups = np.asarray(groups)
        for g in sorted(set(groups)):
            mask = groups == g
            tplt.scatter(x[mask], y[mask], label=str(g), marker="dot")
    tplt.title(title)
    if xlabel:
        tplt.xlabel(xlabel)
    if ylabel:
        tplt.ylabel(ylabel)
    get_console().print(tplt.build())
    tplt.clear_figure()
