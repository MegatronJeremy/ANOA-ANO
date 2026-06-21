"""
Stage 1: QC & preprocessing.

filter cells/genes -> normalize -> log1p -> highly-variable genes.

Design notes (see config.py for the longer version):
- `.X` is reset to raw integer counts on load (io.load_all_samples), so this
  stage's normalize_total/log1p are the *only* normalization applied -- we
  are not stacking scanpy normalization on top of the Seurat normalization
  already baked into the source files.
- HVG selection uses flavor="seurat_v3", which scanpy requires raw counts
  for (it fits a mean-variance trend on counts, not log data), and supports
  batch_key so the HVG set isn't dominated by one sample.
"""
import numpy as np
import scanpy as sc

from . import config as cfg
from .logging_utils import (
    get_logger,
    describe_adata,
    log_percentile_table,
    log_table,
    require,
    term_hist,
    term_bar,
)


def compute_qc_metrics(adata):
    """Flag mitochondrial genes and compute per-cell QC metrics in place."""
    log = get_logger()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_mt = int(adata.var["mt"].sum())
    log.debug(f"flagged {n_mt} mitochondrial genes (var_names starting with 'MT-')")
    require(n_mt > 0, "No mitochondrial genes found (var_names.startswith('MT-')) -- "
                       "pct_counts_mt will be meaningless. Check gene naming convention.")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    return adata


def log_qc_distributions(adata):
    """
    Debug-mode: percentile summary of the 3 QC metrics, pre-filtering, as one
    table, plus a terminal histogram per metric with the configured
    threshold(s) marked as vertical lines -- the main way to actually SEE
    whether a threshold sits sanely on the real distribution rather than
    trusting a printed number.
    """
    log_percentile_table({
        "n_genes_by_counts": adata.obs["n_genes_by_counts"],
        "total_counts": adata.obs["total_counts"],
        "pct_counts_mt": adata.obs["pct_counts_mt"],
    })

    term_hist(
        adata.obs["n_genes_by_counts"],
        title=f"genes/cell (min={cfg.QC_MIN_GENES_PER_CELL}, max={cfg.QC_MAX_GENES_PER_CELL})",
        vlines=[(cfg.QC_MIN_GENES_PER_CELL, "red"), (cfg.QC_MAX_GENES_PER_CELL, "red")],
        xlabel="n_genes_by_counts",
    )
    term_hist(
        adata.obs["total_counts"],
        title="counts/cell (no filter applied on this metric)",
        xlabel="total_counts",
    )
    term_hist(
        adata.obs["pct_counts_mt"],
        title=f"%mito (max={cfg.QC_MAX_PCT_MITO})",
        vlines=[(cfg.QC_MAX_PCT_MITO, "red")],
        xlabel="pct_counts_mt",
    )


def apply_filters(adata):
    """
    Apply the cell/gene QC filters from config.py. Records how many
    cells/genes each individual filter removes (not just the combined
    total), rendered as one summary table, so a threshold that's too
    aggressive is visible immediately rather than buried in one final
    number.
    """
    log = get_logger()
    n0 = adata.n_obs
    sample_counts_before = adata.obs["sample"].value_counts().to_dict()
    rows = []

    sc.pp.filter_cells(adata, min_genes=cfg.QC_MIN_GENES_PER_CELL)
    n1 = adata.n_obs
    rows.append(["min_genes_per_cell", n0 - n1, n1, f">= {cfg.QC_MIN_GENES_PER_CELL}"])

    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=cfg.QC_MIN_CELLS_PER_GENE)
    rows.append(["min_cells_per_gene", n_genes_before - adata.n_vars, adata.n_vars,
                 f">= {cfg.QC_MIN_CELLS_PER_GENE} (genes)"])

    n2 = adata.n_obs
    adata = adata[adata.obs["n_genes_by_counts"] < cfg.QC_MAX_GENES_PER_CELL].copy()
    n3 = adata.n_obs
    rows.append(["max_genes_per_cell", n2 - n3, n3, f"< {cfg.QC_MAX_GENES_PER_CELL}"])

    n4 = adata.n_obs
    adata = adata[adata.obs["pct_counts_mt"] < cfg.QC_MAX_PCT_MITO].copy()
    n5 = adata.n_obs
    rows.append(["max_pct_mito", n4 - n5, n5, f"< {cfg.QC_MAX_PCT_MITO}"])

    log_table("QC filter summary", ["filter", "removed", "remaining", "threshold"], rows)
    log.info(f"QC filtering: {n0} -> {n5} cells ({n0 - n5} removed, {100*(n0-n5)/n0:.1f}%)")

    sample_counts_after = adata.obs["sample"].value_counts().to_dict()
    samples = sorted(set(sample_counts_before) | set(sample_counts_after))
    term_bar(
        samples,
        {
            "before": [sample_counts_before.get(s, 0) for s in samples],
            "after": [sample_counts_after.get(s, 0) for s in samples],
        },
        title="per-sample cell counts: before vs after QC",
        ylabel="n_cells",
    )

    require(adata.n_obs > 0, "AnnData is empty after QC filtering -- thresholds in "
                              "config.py are too aggressive for this dataset.")
    require(adata.n_vars > 0, "AnnData has zero genes after QC filtering.")

    return adata


def normalize_and_log(adata):
    """normalize_total (CP10K) -> log1p. Stash a log-normalized copy in .raw
    for later DE / plotting, which conventionally expect that, not scaled data."""
    log = get_logger()
    sc.pp.normalize_total(adata, target_sum=cfg.NORM_TARGET_SUM)
    sc.pp.log1p(adata)
    adata.raw = adata
    log.debug(f"normalized to {cfg.NORM_TARGET_SUM:.0f} counts/cell and log1p-transformed; "
              f".raw set to this log-normalized matrix")
    return adata


def select_hvgs(adata):
    """
    NOTE: flavor="seurat_v3" fits a mean-variance trend on raw counts and
    requires non-negative integers, not log data. By this point in the
    pipeline .X already holds log-normalized values (normalize_and_log runs
    first so .raw is captured before HVG selection touches anything), so we
    point highly_variable_genes at the untouched `counts` layer explicitly
    rather than relying on .X. Without `layer="counts"` scanpy silently
    proceeds anyway but warns ("expects raw count data, but non-integers
    were found") and the variance-stabilizing fit is wrong -- caught by
    running this stage with --debug on real data (2026-06-16).
    """
    log = get_logger()
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=cfg.N_TOP_GENES,
        flavor=cfg.HVG_FLAVOR,
        batch_key=cfg.BATCH_KEY,
        layer="counts",
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    log.info(f"selected {n_hvg} highly variable genes (flavor={cfg.HVG_FLAVOR}, "
             f"target={cfg.N_TOP_GENES}, batch_key={cfg.BATCH_KEY})")
    require(n_hvg > 0, "No highly variable genes were selected.")
    return adata


_QC_METRICS = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]


def save_qc_figures(obs_before, obs_after, suffix: str = ""):
    """
    Save matplotlib QC figures to results/figures/ (PNG): per-sample violins of
    the 3 QC metrics before and after filtering, and a total_counts vs %mito
    scatter (pre-filter) with the thresholds drawn as lines. These are the
    persisted deliverable figures; the --debug terminal plots are only a quick
    in-session look. `suffix` is "_smoke" on smoke runs so they never overwrite
    the full-run figures.
    """
    import matplotlib
    matplotlib.use("Agg")  # headless: never needs a display, safe under CI/pipe
    import matplotlib.pyplot as plt
    import seaborn as sns

    log = get_logger()

    def _violin(obs, when, fname):
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        for ax, metric in zip(axes, _QC_METRICS):
            sns.violinplot(data=obs, x="sample", y=metric, ax=ax, cut=0, density_norm="width")
            ax.set_title(f"{metric} ({when} QC)")
            ax.tick_params(axis="x", rotation=30)
        fig.tight_layout()
        out = cfg.FIG_DIR / fname
        fig.savefig(out, dpi=120)
        plt.close(fig)
        return out

    f1 = _violin(obs_before, "before", f"01_qc_violin_before{suffix}.png")
    f2 = _violin(obs_after, "after", f"01_qc_violin_after{suffix}.png")

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.scatterplot(data=obs_before, x="total_counts", y="pct_counts_mt",
                    hue="sample", s=6, alpha=0.4, edgecolor=None, ax=ax)
    ax.axhline(cfg.QC_MAX_PCT_MITO, color="red", ls="--", lw=1,
               label=f"max %mito = {cfg.QC_MAX_PCT_MITO}")
    ax.set_title("total_counts vs %mito (pre-filter)")
    ax.legend(markerscale=2, fontsize=8)
    fig.tight_layout()
    f3 = cfg.FIG_DIR / f"01_qc_scatter_counts_mito{suffix}.png"
    fig.savefig(f3, dpi=120)
    plt.close(fig)

    log.info(f"saved QC figures -> {f1.name}, {f2.name}, {f3.name}")


def run(adata, debug: bool = False, smoke: bool = False):
    """Run the full Stage 1 pipeline on an already-loaded, raw-counts AnnData."""
    log = get_logger()
    describe_adata(adata, "qc:input")

    adata = compute_qc_metrics(adata)
    if debug:
        log_qc_distributions(adata)

    # Snapshot the pre-filter QC metrics (for the "before" violin) before
    # apply_filters subsets the cells.
    obs_before = adata.obs[_QC_METRICS + ["sample"]].copy()

    adata = apply_filters(adata)
    save_qc_figures(obs_before, adata.obs[_QC_METRICS + ["sample"]].copy(),
                    suffix="_smoke" if smoke else "")

    adata = normalize_and_log(adata)
    adata = select_hvgs(adata)

    describe_adata(adata, "qc:output")
    return adata
