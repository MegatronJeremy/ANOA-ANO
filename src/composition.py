"""
Stage 4: Composition analysis.

Ported from legacy/notebooks/04_composition.ipynb. Per-sample cell-type
proportions and how they shift relative to the control.

DESIGN NOTE (one donor, one sample per condition -> NO biological replicates):
the PRIMARY result is descriptive -- proportions and log2 fold-change vs
control. Replicate-based tests (propeller / ANOVA over conditions) are NOT
valid here. A per-lineage chi-square is offered only as an EXPLORATORY signal
and is explicitly flagged as pseudoreplication (cells are treated as
independent, which they are not). It is not evidence on its own.
"""
import numpy as np

from . import config as cfg
from .logging_utils import get_logger, describe_adata, log_table, require


def composition_tables(adata):
    """Return (counts, proportions, log2fc_vs_control) DataFrames.
    Rows = samples, columns = lineages. Pure / no I/O, so it is unit-tested."""
    import pandas as pd
    require("lineage" in adata.obs.columns,
            "No 'lineage' column -- run Stage 3 (annotation) first.")
    require(cfg.CONTROL_LABEL in set(adata.obs[cfg.BATCH_KEY]),
            f"Control sample '{cfg.CONTROL_LABEL}' not present -- cannot compute "
            f"relative-to-control composition.")

    counts = pd.crosstab(adata.obs[cfg.BATCH_KEY], adata.obs["lineage"])
    proportions = counts.div(counts.sum(axis=1), axis=0)

    ctrl = proportions.loc[cfg.CONTROL_LABEL]
    with np.errstate(divide="ignore", invalid="ignore"):
        log2fc = np.log2(proportions.div(ctrl, axis=1))
    log2fc = log2fc.replace([np.inf, -np.inf], np.nan)  # 0-proportion -> undefined, not +/-inf
    return counts, proportions, log2fc


def exploratory_shift_tests(counts):
    """EXPLORATORY only (see module docstring). Per exposed sample and lineage,
    a 2x2 chi-square of (this lineage vs rest) x (exposed vs control), BH-FDR
    across lineages within each sample. Returns a tidy DataFrame with a caveat
    column so nobody mistakes it for replicate-backed significance."""
    import pandas as pd
    from scipy.stats import chi2_contingency, false_discovery_control

    ctrl_row = counts.loc[cfg.CONTROL_LABEL]
    ctrl_total = ctrl_row.sum()
    rows = []
    for sample in counts.index:
        if sample == cfg.CONTROL_LABEL:
            continue
        exp_row = counts.loc[sample]
        exp_total = exp_row.sum()
        pvals, lineages = [], []
        for lineage in counts.columns:
            a, c = exp_row[lineage], ctrl_row[lineage]
            table = [[a, exp_total - a], [c, ctrl_total - c]]
            try:
                _, p, _, _ = chi2_contingency(table)
            except ValueError:
                p = np.nan
            pvals.append(p); lineages.append(lineage)
        valid = ~np.isnan(pvals)
        fdr = np.full(len(pvals), np.nan)
        if valid.any():
            fdr[valid] = false_discovery_control(np.array(pvals)[valid], method="bh")
        for lineage, p, q in zip(lineages, pvals, fdr):
            rows.append([sample, lineage, p, q, "exploratory: pseudoreplication, no replicates"])
    return pd.DataFrame(rows, columns=["sample", "lineage", "p_chi2", "fdr_bh", "caveat"])


def save_composition_figures(proportions, suffix: str = ""):
    """Stacked bar (composition per sample) + grouped bar (per-lineage across
    samples), to results/figures/."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    log = get_logger()

    # order samples so control is first, for at-a-glance comparison
    order = [cfg.CONTROL_LABEL] + [s for s in proportions.index if s != cfg.CONTROL_LABEL]
    prop = proportions.loc[order]

    fig, ax = plt.subplots(figsize=(8, 5))
    prop.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", width=0.8)
    ax.set_ylabel("proportion of cells"); ax.set_title("Cell-type composition per sample")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7, frameon=False)
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(cfg.FIG_DIR / f"04_composition_stacked{suffix}.png", dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    prop.T.plot(kind="bar", ax=ax, width=0.8)
    ax.set_ylabel("proportion of cells"); ax.set_title("Per-lineage proportion across samples")
    ax.legend(fontsize=7, frameon=False)
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(cfg.FIG_DIR / f"04_composition_grouped{suffix}.png", dpi=120)
    plt.close(fig)
    log.info(f"saved composition figures -> 04_composition_stacked{suffix}.png, "
             f"04_composition_grouped{suffix}.png")


def run(adata, debug: bool = False, smoke: bool = False):
    """Run Stage 4 on a Stage-3 (annotated) checkpoint AnnData. Writes tables +
    figures only -- no new AnnData checkpoint."""
    log = get_logger()
    describe_adata(adata, "composition:input")
    sfx = "_smoke" if smoke else ""

    counts, proportions, log2fc = composition_tables(adata)
    counts.to_csv(cfg.TABLE_DIR / f"04_composition_counts{sfx}.csv")
    proportions.to_csv(cfg.TABLE_DIR / f"04_composition_proportions{sfx}.csv")
    log2fc.to_csv(cfg.TABLE_DIR / f"04_composition_relative_to_control{sfx}.csv")

    tests = exploratory_shift_tests(counts)
    tests.to_csv(cfg.TABLE_DIR / f"04_composition_exploratory_tests{sfx}.csv", index=False)

    log_table("Cell-type proportions per sample (%)",
              ["sample"] + list(proportions.columns),
              [[idx] + [f"{proportions.loc[idx, c]:.1%}" for c in proportions.columns]
               for idx in proportions.index], debug_only=False)
    log_table("log2 fold-change vs control",
              ["sample"] + list(log2fc.columns),
              [[idx] + [("." if np.isnan(log2fc.loc[idx, c]) else f"{log2fc.loc[idx, c]:+.2f}")
                        for c in log2fc.columns]
               for idx in log2fc.index if idx != cfg.CONTROL_LABEL], debug_only=False)

    save_composition_figures(proportions, suffix=sfx)
    return adata
