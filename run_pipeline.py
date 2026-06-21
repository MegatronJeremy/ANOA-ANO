#!/usr/bin/env python
"""
Driver for the nanoplastic scRNA-seq pipeline.

Headless (canonical, scriptable/CI-safe -- this is what reproducibility and
grading depend on):
    python run_pipeline.py --stage qc                  # run Stage 1 only
    python run_pipeline.py --stage qc --debug           # ... with rich debug logging
    python run_pipeline.py --stage qc --smoke-test       # ... on a small subsample
    python run_pipeline.py --stage qc --smoke-test --debug
    python run_pipeline.py --stage qc --subsample 200 --debug

Interactive (convenience only):
    python run_pipeline.py                              # launches a menu to
                                                          # pick a stage + toggle
                                                          # options, if stdin/stdout
                                                          # is a real terminal.
                                                          # Otherwise prints help
                                                          # and exits -- never blocks
                                                          # a non-interactive run.

Both paths call the exact same per-stage run function (see STAGE_REGISTRY) --
the menu is sugar on top of the CLI, not a separate code path.

Stages are added incrementally; only "qc" exists today. Each stage reads the
previous stage's checkpoint from data/processed/ (or, for "qc", loads the
raw samples) and writes its own checkpoint for the next stage to consume.
"""
import argparse
import random
import sys
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from src import config as cfg
from src import io as pio
from src import qc
from src import integration
from src import annotation
from src.logging_utils import setup_logging, log_stage_banner, progress_spinner, get_console, get_logger


def _fix_console_encoding():
    """
    On this Windows box, sys.stdout/stderr can report a legacy codepage
    (e.g. cp1252) even when the terminal itself is UTF-8, which crashes the
    moment rich/plotext print a box-drawing or other non-ASCII character.
    Force UTF-8 with a safe fallback so output never crashes on encoding,
    regardless of how the process was launched (terminal, piped, redirected).
    """
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")


def set_seeds(seed: int):
    random.seed(seed)
    np.random.seed(seed)


def run_qc_stage(args, log):
    is_smoke = bool(args.smoke_test or args.subsample)
    with log_stage_banner("Stage 1: QC & preprocessing"):
        adata = pio.load_all_samples()
        if is_smoke:
            n = args.subsample or cfg.SMOKE_TEST_CELLS_PER_SAMPLE
            adata = pio.subsample_for_smoke_test(adata, n_per_sample=n, seed=cfg.RANDOM_SEED)
            adata = qc.run(adata, debug=args.debug, smoke=True)
        else:
            # Full runs take ~minutes; smoke tests are seconds and don't need a spinner.
            with progress_spinner("Running QC filtering, normalization & HVG selection"):
                adata = qc.run(adata, debug=args.debug, smoke=False)
        checkpoint_name = "01_qc_done_smoke" if is_smoke else "01_qc_done"
        pio.save_checkpoint(adata, checkpoint_name)
    return adata


def run_integration_stage(args, log):
    is_smoke = bool(args.smoke_test or args.subsample)
    input_name = "01_qc_done_smoke" if is_smoke else "01_qc_done"
    with log_stage_banner("Stage 2: Integration & clustering"):
        adata = pio.load_checkpoint(input_name)
        if is_smoke:
            adata = integration.run(adata, debug=args.debug, smoke=True)
        else:
            # PCA + Harmony + neighbors/UMAP (x2) + Leiden take a few minutes on the full dataset.
            with progress_spinner("Running PCA, Harmony integration, UMAP & Leiden clustering"):
                adata = integration.run(adata, debug=args.debug, smoke=False)
        checkpoint_name = "02_clustered_smoke" if is_smoke else "02_clustered"
        pio.save_checkpoint(adata, checkpoint_name)
    return adata


def run_annotation_stage(args, log):
    is_smoke = bool(args.smoke_test or args.subsample)
    input_name = "02_clustered_smoke" if is_smoke else "02_clustered"
    with log_stage_banner("Stage 3: Cell-type annotation"):
        adata = pio.load_checkpoint(input_name)
        if is_smoke:
            adata = annotation.run(adata, debug=args.debug, smoke=True)
        else:
            with progress_spinner("Running celltypist annotation & cross-check"):
                adata = annotation.run(adata, debug=args.debug, smoke=False)
        checkpoint_name = "03_annotated_smoke" if is_smoke else "03_annotated"
        pio.save_checkpoint(adata, checkpoint_name)
    return adata


# ---------------------------------------------------------------------------
# Stage registry -- the single source of truth for both the --stage CLI
# argument and the interactive menu. Adding a future stage (2-6) is a
# one-line entry here; nothing else needs to change for it to show up in
# both interfaces.
# ---------------------------------------------------------------------------

@dataclass
class StageSpec:
    key: str                          # value passed to --stage
    label: str                        # short name shown in the menu
    description: str                  # one-line description shown in the menu
    input_checkpoint: Optional[str]   # checkpoint this stage reads, or None if it loads raw data
    output_checkpoint: str            # checkpoint name this stage writes (full-run name)
    run_fn: Callable                  # run_fn(args, log) -> AnnData; same function for CLI and menu


STAGE_REGISTRY = {
    "qc": StageSpec(
        key="qc",
        label="Stage 1: QC & preprocessing",
        description="filter cells/genes, normalize, log1p, select HVGs",
        input_checkpoint=None,   # loads the 4 raw samples directly, not a checkpoint
        output_checkpoint="01_qc_done",
        run_fn=run_qc_stage,
    ),
    "integration": StageSpec(
        key="integration",
        label="Stage 2: Integration & clustering",
        description="Harmony batch correction, UMAP, Leiden clusters",
        input_checkpoint="01_qc_done",
        output_checkpoint="02_clustered",
        run_fn=run_integration_stage,
    ),
    "annotation": StageSpec(
        key="annotation",
        label="Stage 3: Cell-type annotation",
        description="celltypist + marker genes -> lineage labels, cross-checked vs Azimuth & CoDi",
        input_checkpoint="02_clustered",
        output_checkpoint="03_annotated",
        run_fn=run_annotation_stage,
    ),
    # Stage 4-6 register here, one line each.
}


def stage_status(spec: StageSpec) -> dict:
    """
    Checkpoint-based status for the menu: is this stage's input ready, and
    has it already produced a (full or smoke) output checkpoint?
    """
    if spec.input_checkpoint is None:
        input_ready = all((cfg.RAW_DIR / fname).exists() for fname in cfg.SAMPLES.values())
        input_desc = "raw data files"
    else:
        input_ready = (cfg.PROCESSED_DIR / f"{spec.input_checkpoint}.h5ad").exists()
        input_desc = f"{spec.input_checkpoint}.h5ad"

    full_done = (cfg.PROCESSED_DIR / f"{spec.output_checkpoint}.h5ad").exists()
    smoke_done = (cfg.PROCESSED_DIR / f"{spec.output_checkpoint}_smoke.h5ad").exists()
    return {
        "input_ready": input_ready,
        "input_desc": input_desc,
        "full_done": full_done,
        "smoke_done": smoke_done,
    }


def run_stage(spec: StageSpec, args, log):
    """The one code path both the CLI and the menu call to actually run a stage."""
    try:
        spec.run_fn(args, log)
        return True
    except (FileNotFoundError, AssertionError) as exc:
        log.error(f"Pipeline stopped: {exc}")
        return False


# ---------------------------------------------------------------------------
# Environment doctor (`--check`): a quick, side-effect-free readiness report.
# Verifies the things a run depends on -- Python version, raw data files,
# existing checkpoints, and importability of the heavy libraries -- without
# loading any data or running an analysis. Mirrors the PSZ `debug` command.
# ---------------------------------------------------------------------------

def run_doctor() -> bool:
    from rich.table import Table
    from rich.panel import Panel
    console = get_console()
    console.print(Panel("[bold]Environment check[/bold]", title="DOCTOR", border_style="cyan"))

    ok = True

    # Python
    py = ".".join(str(v) for v in sys.version_info[:3])
    py_ok = sys.version_info[:2] >= (3, 10)
    ok = ok and py_ok

    # Raw data files
    data_tbl = Table(title="Raw data files (data/raw/)", header_style="bold cyan")
    data_tbl.add_column("sample"); data_tbl.add_column("file"); data_tbl.add_column("present")
    for label, fname in cfg.SAMPLES.items():
        path = cfg.RAW_DIR / fname
        present = path.exists()
        ok = ok and present
        size = f"{path.stat().st_size / 1e6:.0f} MB" if present else "-"
        data_tbl.add_row(label, fname, f"[green]yes[/green] ({size})" if present else "[red]MISSING[/red]")

    # Checkpoints
    ckpts = sorted(p.name for p in cfg.PROCESSED_DIR.glob("*.h5ad"))
    ckpt_str = ", ".join(ckpts) if ckpts else "[yellow]none yet[/yellow]"

    # Key libraries
    libs = ["scanpy", "anndata", "harmonypy", "leidenalg", "celltypist",
            "gseapy", "skmisc", "rich", "plotext", "questionary", "pytest"]
    lib_tbl = Table(title="Key libraries", header_style="bold cyan")
    lib_tbl.add_column("library"); lib_tbl.add_column("import")
    import importlib.util
    for name in libs:
        found = importlib.util.find_spec(name) is not None
        ok = ok and found
        lib_tbl.add_row(name, "[green]ok[/green]" if found else "[red]MISSING[/red]")

    console.print(f"Python: {py}  " + ("[green](>= 3.10, ok)[/green]" if py_ok else "[red](need >= 3.10)[/red]"))
    console.print(data_tbl)
    console.print(f"Checkpoints in data/processed/: {ckpt_str}")
    console.print(lib_tbl)
    console.print(Panel("[bold green]All good.[/bold green]" if ok else "[bold red]Problems found -- see above.[/bold red]",
                        border_style="green" if ok else "red"))
    return ok


# ---------------------------------------------------------------------------
# Interactive menu (convenience wrapper only -- see module docstring)
# ---------------------------------------------------------------------------

def _is_interactive() -> bool:
    return sys.stdin.isatty() and sys.stdout.isatty()


def _print_status_table():
    from rich.table import Table
    console = get_console()
    table = Table(title="Pipeline stages", show_header=True, header_style="bold cyan")
    table.add_column("stage")
    table.add_column("description")
    table.add_column("input ready?")
    table.add_column("full run done?")
    table.add_column("smoke run done?")
    for spec in STAGE_REGISTRY.values():
        st = stage_status(spec)
        table.add_row(
            spec.label,
            spec.description,
            f"[green]yes[/green] ({st['input_desc']})" if st["input_ready"] else f"[red]no[/red] ({st['input_desc']})",
            "[green]yes[/green]" if st["full_done"] else "[yellow]no[/yellow]",
            "[green]yes[/green]" if st["smoke_done"] else "[yellow]no[/yellow]",
        )
    console.print(table)


def run_menu():
    """
    Interactive launcher: pick a stage, toggle smoke-test/debug/subsample,
    run it through the exact same run_stage()/run_fn used by the headless
    CLI, then loop back to the menu. Never reached for non-TTY invocations.
    """
    import questionary

    console = get_console()
    log = get_logger()

    while True:
        console.print()
        _print_status_table()

        choices = [
            questionary.Choice(title=f"{spec.label} -- {spec.description}", value=spec.key)
            for spec in STAGE_REGISTRY.values()
        ]
        choices.append(questionary.Choice(title="Quit", value=None))

        selected = questionary.select("Select a stage to run:", choices=choices).ask()
        if selected is None:
            console.print("[bold]Bye.[/bold]")
            return

        spec = STAGE_REGISTRY[selected]

        smoke_test = questionary.confirm(
            "Run as a smoke test (small random subsample)?", default=False
        ).ask()
        if smoke_test is None:  # Ctrl+C
            continue

        subsample = None
        if smoke_test:
            raw = questionary.text(
                f"Cells per sample (default {cfg.SMOKE_TEST_CELLS_PER_SAMPLE}):",
                default=str(cfg.SMOKE_TEST_CELLS_PER_SAMPLE),
            ).ask()
            if raw is None:
                continue
            n = int(raw) if raw.strip() else cfg.SMOKE_TEST_CELLS_PER_SAMPLE
            if n != cfg.SMOKE_TEST_CELLS_PER_SAMPLE:
                # mirrors --subsample N on the CLI: custom size, smoke_test flag itself off
                smoke_test, subsample = False, n

        debug = questionary.confirm("Enable --debug output?", default=False).ask()
        if debug is None:
            continue

        args = argparse.Namespace(smoke_test=smoke_test, debug=debug, subsample=subsample)

        log.handlers.clear()
        setup_logging(debug=debug)
        log = get_logger()
        set_seeds(cfg.RANDOM_SEED)

        run_stage(spec, args, log)

        questionary.text("Press Enter to return to the menu...").ask()


# ---------------------------------------------------------------------------
# Headless CLI entry point
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stage", choices=list(STAGE_REGISTRY.keys()), default=None,
                         help="Which stage to run. Omit (with no other flags) to launch the "
                              "interactive menu in a real terminal.")
    parser.add_argument("--debug", action="store_true",
                         help="Verbose per-stage debug logging.")
    parser.add_argument("--smoke-test", action="store_true",
                         help=f"Run on a random subsample (~{cfg.SMOKE_TEST_CELLS_PER_SAMPLE} "
                              f"cells/sample) to verify the pipeline executes end-to-end quickly.")
    parser.add_argument("--subsample", type=int, default=None, metavar="N",
                         help="Like --smoke-test but with a custom N cells/sample.")
    parser.add_argument("--check", action="store_true",
                         help="Environment doctor: verify Python, raw data, checkpoints and "
                              "libraries without running any analysis, then exit.")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    _fix_console_encoding()

    if args.check:
        setup_logging(debug=False)
        sys.exit(0 if run_doctor() else 1)

    if args.stage is None:
        # No --stage given. Only launch the menu in a real interactive terminal;
        # never block a piped/redirected/CI invocation.
        if _is_interactive():
            setup_logging(debug=False)
            run_menu()
            return
        else:
            parser.print_help()
            sys.exit(0)

    # Headless path -- unchanged behaviour, no menu/prompt involved.
    log = setup_logging(debug=args.debug)
    set_seeds(cfg.RANDOM_SEED)

    log.info(f"random seed = {cfg.RANDOM_SEED}")
    if args.smoke_test or args.subsample:
        log.info(f"SMOKE TEST MODE: subsampling to {args.subsample or cfg.SMOKE_TEST_CELLS_PER_SAMPLE} cells/sample")

    ok = run_stage(STAGE_REGISTRY[args.stage], args, log)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
