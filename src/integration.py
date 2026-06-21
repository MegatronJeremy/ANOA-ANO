"""
Stage 2: Integration & clustering.

PCA -> neighbors+UMAP on raw PCA ("before", for comparison) -> Harmony batch
correction -> neighbors+UMAP on the Harmony embedding ("after", the one
downstream stages use) -> Leiden clustering.

Design notes:
- PCA is computed on the HVG subset scanpy already flagged in Stage 1
  (`sc.pp.pca` automatically restricts to `.var['highly_variable']` when that
  column is present -- no extra step needed here).
- The "before" UMAP exists purely so the before/after batch-mixing comparison
  in --debug is meaningful: scanpy always writes neighbors/UMAP results to
  the same default keys, so we copy the pre-Harmony result out to
  `obsm['X_umap_pre_harmony']` before recomputing on the corrected embedding.
  `obsm['X_umap']` and `obsm['X_pca_harmony']` end up holding the
  post-integration state that Stage 3+ should use.
- n_comps/n_neighbors adapt downward (with a warning) if a run has fewer
  cells than the configured value -- relevant for small --subsample smoke
  runs, not the full dataset.
"""
import warnings

import numpy as np
import scanpy as sc

from . import config as cfg
from .logging_utils import (
    get_logger,
    describe_adata,
    log_table,
    require,
    silence_third_party_logger,
    term_bar,
    term_scatter,
)


def run_pca(adata):
    log = get_logger()
    require("highly_variable" in adata.var.columns,
             "No 'highly_variable' column found -- run Stage 1 (qc) first.")

    n_comps = min(cfg.N_PCS, adata.n_obs - 1, int(adata.var["highly_variable"].sum()) - 1)
    if n_comps < cfg.N_PCS:
        log.warning(f"n_comps reduced from {cfg.N_PCS} to {n_comps} "
                    f"(too few cells/HVGs for the configured value -- expected on small subsamples)")

    sc.pp.pca(adata, n_comps=n_comps, random_state=cfg.RANDOM_SEED)
    log.debug(f"PCA: {n_comps} components on the HVG subset")
    return adata


def _neighbors_and_umap(adata, use_rep=None, n_pcs=None):
    n_neighbors = min(cfg.N_NEIGHBORS, adata.n_obs - 1)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep,
                     random_state=cfg.RANDOM_SEED)
    sc.tl.umap(adata, random_state=cfg.RANDOM_SEED)
    return adata


def compute_pre_integration_umap(adata, debug: bool = False):
    """UMAP on the raw (uncorrected) PCA -- the 'before' picture, kept only
    for the --debug batch-mixing comparison plot."""
    log = get_logger()
    adata = _neighbors_and_umap(adata, use_rep="X_pca")
    adata.obsm["X_umap_pre_harmony"] = adata.obsm["X_umap"].copy()
    log.debug("computed pre-integration UMAP -> obsm['X_umap_pre_harmony']")

    if debug:
        term_scatter(
            adata.obsm["X_umap_pre_harmony"][:, 0],
            adata.obsm["X_umap_pre_harmony"][:, 1],
            groups=adata.obs[cfg.BATCH_KEY].astype(str).values,
            title="UMAP before integration (coloured by sample)",
            xlabel="UMAP1", ylabel="UMAP2",
        )
    return adata


def integrate_harmony(adata, debug: bool = False):
    """
    Harmony batch correction on the BATCH_KEY column, writing
    `obsm['X_pca_harmony']`. Guards that the batch key exists, is present in
    every cell, and actually has >1 category -- Harmony on a single batch is
    a no-op that would silently hide a broken batch_key upstream.
    """
    log = get_logger()
    require(cfg.BATCH_KEY in adata.obs.columns,
             f"Batch key '{cfg.BATCH_KEY}' not found in .obs -- check Stage 1 output.")
    n_batches = adata.obs[cfg.BATCH_KEY].nunique()
    require(n_batches > 1,
             f"Batch key '{cfg.BATCH_KEY}' has only {n_batches} unique value(s) -- "
             f"nothing to integrate. Check sample loading in io.load_all_samples().")
    require(adata.obs[cfg.BATCH_KEY].isna().sum() == 0,
             f"Batch key '{cfg.BATCH_KEY}' has missing values.")

    # harmonypy sets its own logger to DEBUG + attaches its own StreamHandler
    # at *import* time, independent of our logging setup. sc.external.pp.
    # harmony_integrate() imports it lazily on first call, which would
    # silently re-clobber a level we set beforehand -- so import explicitly
    # first, then silence, so our setting is the one that sticks.
    import harmonypy  # noqa: F401
    silence_third_party_logger("harmonypy", debug=debug)
    sc.external.pp.harmony_integrate(adata, key=cfg.BATCH_KEY, random_state=cfg.RANDOM_SEED)

    require("X_pca_harmony" in adata.obsm,
             "Harmony did not produce 'X_pca_harmony' -- integration failed silently.")
    log.info(f"Harmony integration complete ({n_batches} batches: {cfg.BATCH_KEY})")
    return adata


def compute_post_integration_umap(adata, debug: bool = False):
    """UMAP on the Harmony-corrected embedding -- the embedding Stage 3+ use."""
    log = get_logger()
    adata = _neighbors_and_umap(adata, use_rep="X_pca_harmony")
    log.debug("computed post-integration UMAP -> obsm['X_umap'] (neighbors graph built on X_pca_harmony)")

    if debug:
        term_scatter(
            adata.obsm["X_umap"][:, 0],
            adata.obsm["X_umap"][:, 1],
            groups=adata.obs[cfg.BATCH_KEY].astype(str).values,
            title="UMAP after integration (coloured by sample)",
            xlabel="UMAP1", ylabel="UMAP2",
        )
    return adata


def cluster_leiden(adata, debug: bool = False):
    """
    Leiden clustering on the (Harmony-corrected) neighbor graph. flavor is
    pinned explicitly to "leidenalg" -- scanpy's current default -- since
    scanpy unconditionally raises a FutureWarning about the future default
    changing to "igraph" regardless of what's passed; we deliberately keep
    today's default rather than switching algorithms as a side effect of
    silencing a warning, so the warning is suppressed here instead.
    """
    log = get_logger()
    with warnings.catch_warnings():
        # NOTE: scanpy reports this warning's origin at the *caller's*
        # stacklevel (i.e. this module, not scanpy's), so a module= filter
        # here would not match -- confirmed by testing. Scoped to this one
        # call only, so it can't hide an unrelated FutureWarning elsewhere.
        warnings.filterwarnings("ignore", category=FutureWarning,
                                 message=".*default backend for leiden.*")
        sc.tl.leiden(adata, resolution=cfg.LEIDEN_RES, random_state=cfg.RANDOM_SEED,
                     flavor="leidenalg")

    sizes = adata.obs["leiden"].value_counts().sort_index()
    n_clusters = len(sizes)
    n_total = int(sizes.sum())
    require(n_clusters > 1, "Leiden produced only 1 cluster -- check the neighbor graph / "
                             "resolution; LEIDEN_RES in config.py may need adjusting.")

    log.info(f"Leiden clustering: {n_clusters} clusters (resolution={cfg.LEIDEN_RES})")
    log_table("Cluster sizes", ["cluster", "n_cells", "fraction"],
              [[c, n, f"{n / n_total:.1%}"] for c, n in sizes.items()])

    tiny = sizes[sizes / n_total < cfg.CLUSTER_TINY_FRACTION]
    if len(tiny):
        log.warning(f"{len(tiny)} cluster(s) have < {cfg.CLUSTER_TINY_FRACTION:.0%} of all cells: "
                    f"{dict(tiny)} -- could be a rare population or an artifact, worth a manual look")
    dominant = sizes[sizes / n_total > cfg.CLUSTER_DOMINANT_FRACTION]
    if len(dominant):
        log.warning(f"{len(dominant)} cluster(s) have > {cfg.CLUSTER_DOMINANT_FRACTION:.0%} of all cells: "
                    f"{dict(dominant)} -- one cluster may be swallowing distinct cell types; "
                    f"consider raising LEIDEN_RES")

    if debug:
        term_bar(list(sizes.index), {"n_cells": list(sizes.values)}, title="Leiden cluster sizes")

    return adata


def cluster_sample_composition(adata):
    """
    Quantify batch mixing: per Leiden cluster, the fraction coming from each
    sample. Good integration means most clusters are a blend of all four
    samples; a cluster that is almost entirely one sample is flagged (either a
    genuinely condition-specific population or a sign integration failed).
    Returns the fraction DataFrame (clusters x samples) and also saves it as a
    CSV deliverable.
    """
    import pandas as pd
    log = get_logger()
    counts = pd.crosstab(adata.obs["leiden"], adata.obs[cfg.BATCH_KEY])
    frac = counts.div(counts.sum(axis=1), axis=0)

    log_table(
        "Per-cluster sample composition (fraction)",
        ["cluster"] + list(frac.columns),
        [[idx] + [f"{frac.loc[idx, c]:.2f}" for c in frac.columns] for idx in frac.index],
    )

    dominant = frac.max(axis=1)
    impure = dominant[dominant > cfg.CLUSTER_BATCH_PURITY_WARN]
    for cl, val in impure.items():
        top_sample = frac.loc[cl].idxmax()
        log.warning(f"cluster {cl} is {val:.0%} '{top_sample}' (> "
                    f"{cfg.CLUSTER_BATCH_PURITY_WARN:.0%}): check whether this is a real "
                    f"condition-specific population or a batch-correction failure")
    return frac


def save_integration_figures(adata, suffix: str = ""):
    """
    Save UMAP PNGs to results/figures/: before vs after Harmony coloured by
    sample (the batch-correction effect), and the post-integration UMAP
    coloured by Leiden cluster. `suffix` is "_smoke" on smoke runs.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    log = get_logger()

    def _umap(coords, labels, title, fname):
        fig, ax = plt.subplots(figsize=(7, 6))
        labels = np.asarray(labels).astype(str)
        for g in sorted(set(labels), key=lambda x: (len(x), x)):
            m = labels == g
            ax.scatter(coords[m, 0], coords[m, 1], s=4, alpha=0.5, label=g, edgecolor="none")
        ax.set_title(title); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        ax.legend(markerscale=3, fontsize=7, ncol=2, frameon=False)
        fig.tight_layout()
        fig.savefig(cfg.FIG_DIR / fname, dpi=120)
        plt.close(fig)
        return fname

    sample = adata.obs[cfg.BATCH_KEY].values
    f1 = _umap(adata.obsm["X_umap_pre_harmony"], sample,
               "UMAP before Harmony (by sample)", f"02_umap_pre_harmony_by_sample{suffix}.png")
    f2 = _umap(adata.obsm["X_umap"], sample,
               "UMAP after Harmony (by sample)", f"02_umap_by_sample{suffix}.png")
    f3 = _umap(adata.obsm["X_umap"], adata.obs["leiden"].values,
               "UMAP after Harmony (by cluster)", f"02_umap_by_cluster{suffix}.png")
    log.info(f"saved integration figures -> {f1}, {f2}, {f3}")


def run(adata, debug: bool = False, smoke: bool = False):
    """Run the full Stage 2 pipeline on a Stage-1 checkpoint AnnData."""
    describe_adata(adata, "integration:input")

    adata = run_pca(adata)
    adata = compute_pre_integration_umap(adata, debug=debug)
    adata = integrate_harmony(adata, debug=debug)
    adata = compute_post_integration_umap(adata, debug=debug)
    adata = cluster_leiden(adata, debug=debug)

    frac = cluster_sample_composition(adata)
    frac.to_csv(cfg.TABLE_DIR / f"02_cluster_composition{'_smoke' if smoke else ''}.csv")
    save_integration_figures(adata, suffix="_smoke" if smoke else "")

    describe_adata(adata, "integration:output")
    return adata
