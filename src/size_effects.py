"""
Stage 6: Size-specific effects -- the main biological conclusion.

Takes the per-(lineage, condition) DE results from Stage 5 and, for each
lineage, splits the responding genes into four set-membership categories:

  - unique_40nm        : significant in 40 nm vs control, NOT in 200 nm
  - unique_200nm       : significant in 200 nm vs control, NOT in 40 nm
  - shared             : significant in BOTH single sizes
  - mixture_emergent   : significant in the mixture, but in NEITHER single size
                         (a response that only appears when both sizes are present)

The set logic is the load-bearing, easy-to-get-subtly-wrong part, so it lives
in one pure function (`categorize`) covered by tests for all four categories.
This stage reads the 05_DE_*.csv tables (not an AnnData checkpoint) and writes
tables + UpSet figures.
"""
from . import config as cfg
from .logging_utils import get_logger, log_table
from .differential_expression import significant_genes, _safe


def categorize(genes_40, genes_200, genes_mix):
    """Pure set logic. Returns dict category -> sorted gene list."""
    s40, s200, smix = set(genes_40), set(genes_200), set(genes_mix)
    return {
        "unique_40nm": sorted(s40 - s200),
        "unique_200nm": sorted(s200 - s40),
        "shared": sorted(s40 & s200),
        "mixture_emergent": sorted(smix - s40 - s200),
    }


def load_de_results(smoke: bool = False):
    """Read Stage 5's per-(lineage, condition) DE tables back in as
    {(lineage, condition): df}."""
    import pandas as pd
    sfx = "_smoke" if smoke else ""
    conds = [s for s in cfg.SAMPLES if s != cfg.CONTROL_LABEL]
    out = {}
    for f in cfg.TABLE_DIR.glob(f"05_DE_*_vs_control{sfx}.csv"):
        stem = f.name[len("05_DE_"):].rsplit("_vs_control", 1)[0]  # "{safelin}_{cond}"
        for cond in conds:
            if stem.endswith("_" + cond):
                lineage = stem[: -(len(cond) + 1)]
                out[(lineage, cond)] = pd.read_csv(f)
                break
    return out


def save_category_figure(summary, suffix=""):
    """Robust (matplotlib-only) grouped bar chart: the four size-specific
    categories per lineage. This is the reliable deliverable figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    log = get_logger()
    cats = ["unique_40nm", "unique_200nm", "shared", "mixture_emergent"]
    ax = summary.set_index("lineage")[cats].plot(kind="bar", figsize=(9, 5), width=0.8)
    ax.set_ylabel("number of DE genes"); ax.set_title("Size-specific DE genes per lineage")
    ax.tick_params(axis="x", rotation=20)
    ax.legend(fontsize=8, frameon=False)
    fig = ax.get_figure(); fig.tight_layout()
    fig.savefig(cfg.FIG_DIR / f"06_size_categories{suffix}.png", dpi=120)
    plt.close(fig)
    log.info(f"saved size-category figure -> 06_size_categories{suffix}.png")


def save_upset(sets_by_size, lineage, suffix=""):
    """Best-effort UpSet plot per lineage. upsetplot can be incompatible with
    some matplotlib/numpy combos, so failures are caught -- the bar chart and
    the CSVs are the guaranteed deliverables, this is a bonus."""
    log = get_logger()
    contents = {k: set(v) for k, v in sets_by_size.items() if len(v)}
    if len(contents) < 2:
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from upsetplot import from_contents, UpSet
        data = from_contents(contents)
        UpSet(data, subset_size="count", show_counts=True).plot()
        fig = plt.gcf()
        fig.suptitle(f"Size-specific DE genes -- {lineage}")
        fig.savefig(cfg.FIG_DIR / f"06_upset_{_safe(lineage)}{suffix}.png", dpi=110, bbox_inches="tight")
        plt.close(fig)
    except Exception as exc:
        log.warning(f"UpSet plot skipped for {lineage} ({type(exc).__name__}); "
                    f"see 06_size_categories figure instead.")


def run(adata=None, debug: bool = False, smoke: bool = False):
    """Run Stage 6. `adata` is accepted for a uniform stage signature but unused
    -- the inputs are Stage 5's DE CSVs."""
    import pandas as pd
    log = get_logger()
    sfx = "_smoke" if smoke else ""
    de = load_de_results(smoke=smoke)
    if not de:
        log.error("No Stage 5 DE tables found -- run `--stage de` first.")
        raise FileNotFoundError("missing 05_DE_*_vs_control CSVs")

    cond40 = next(c for c in cfg.SAMPLES if c.endswith("40nm"))
    cond200 = next(c for c in cfg.SAMPLES if c.endswith("200nm"))
    condmix = next(c for c in cfg.SAMPLES if "mix" in c.lower())

    lineages = sorted({lin for (lin, _) in de})
    long_rows, summary_rows, sets_per_lineage = [], [], {}
    for lin in lineages:
        def sig(cond):
            df = de.get((lin, cond))
            return significant_genes(df, direction="any") if df is not None else []

        g40, g200, gmix = sig(cond40), sig(cond200), sig(condmix)
        cats = categorize(g40, g200, gmix)
        for cat, genes in cats.items():
            for g in genes:
                long_rows.append([lin, cat, g])
        summary_rows.append([lin, len(cats["unique_40nm"]), len(cats["unique_200nm"]),
                             len(cats["shared"]), len(cats["mixture_emergent"])])
        sets_per_lineage[lin] = {"40nm": g40, "200nm": g200, "mixture": gmix}

    # Write the analysis (the deliverable) BEFORE any plotting, so a figure-
    # library hiccup can never lose the results.
    pd.DataFrame(long_rows, columns=["lineage", "category", "gene"]).to_csv(
        cfg.TABLE_DIR / f"06_size_specific_genes{sfx}.csv", index=False)
    summary = pd.DataFrame(summary_rows,
                           columns=["lineage", "unique_40nm", "unique_200nm", "shared", "mixture_emergent"])
    summary.to_csv(cfg.TABLE_DIR / f"06_size_specific_summary{sfx}.csv", index=False)

    log_table("Size-specific DE gene counts per lineage",
              list(summary.columns), summary.values.tolist(), debug_only=False)

    save_category_figure(summary, suffix=sfx)
    for lin, sets in sets_per_lineage.items():
        save_upset(sets, lin, suffix=sfx)

    # Headline biological readout: where does the mixture do something neither size does?
    emergent = summary.sort_values("mixture_emergent", ascending=False)
    top = emergent.iloc[0]
    log.info(f"BIOLOGICAL READOUT: mixture-emergent response is largest in '{top['lineage']}' "
             f"({int(top['mixture_emergent'])} genes significant only in the mixture). "
             f"200nm-unique vs 40nm-unique totals: "
             f"{int(summary['unique_200nm'].sum())} vs {int(summary['unique_40nm'].sum())}.")
    return adata
