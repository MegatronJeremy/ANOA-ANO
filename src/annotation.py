"""
Stage 3: Cell-type annotation.

Ported from legacy/notebooks/03_annotation.ipynb into a clean, registered
stage. Two signals, used together:

1. celltypist (Immune_All_Low PBMC model, majority-voting per Leiden cluster) --
   the practical Python stand-in for Azimuth, which ships as an R/RDS object.
2. canonical marker genes (config.MARKER_GENES) for a manual sanity dotplot.

Differential verification: we have TWO independent references already in the
data -- `obs['predicted.celltype']` (Azimuth, run in R) and the per-cell CoDi
labels in `data/raw/*_CoDi_KLD.csv`. All three vocabularies are harmonized to a
common coarse lineage and cross-tabulated, so disagreement is reported, not
hidden. There is no single ground truth; three-way agreement is the check.
"""
import re

import numpy as np

from . import config as cfg
from .logging_utils import get_logger, describe_adata, log_table, require


# ---------------------------------------------------------------------------
# Label harmonization. celltypist (fine immune labels), Azimuth
# (predicted.celltype, e.g. "CD4+ T cell") and CoDi (e.g. "Cytotoxic T cell")
# each use a different vocabulary. Collapse any of them to one coarse lineage
# so the three can be compared. Order matters: "plasmacytoid dendritic cell"
# must hit dendritic before "plasma" would send it to the B lineage.
# ---------------------------------------------------------------------------
# Ordered (lineage, regex) rules. Word boundaries (\b) matter: without them
# "mast cells" matches the "t cell" inside "mas-t cell-s" and is mislabelled a
# T cell. Order matters too: dendritic is checked before B so "plasmacytoid
# dendritic cell" lands in DC, not B (via "plasma").
_LINEAGE_RULES = [
    ("Dendritic cell", re.compile(r"dendritic|\bp?dc\b|\bdc[12]\b|asdc")),
    ("Platelet",       re.compile(r"platelet|megakaryo")),
    ("NK cell",        re.compile(r"\bnk\b|natural killer|killer cell")),
    ("Monocyte",       re.compile(r"monocyt|macrophage")),
    ("T cell",         re.compile(r"\bt[\s\-]?cells?\b|\bcd4|\bcd8|helper t|cytotoxic t|"
                                  r"regulatory t|\btreg|\bmait\b|thymocyte")),
    ("B cell",         re.compile(r"\bb[\s\-]?cells?\b|plasma")),
]


def to_lineage(label) -> str:
    """Map any cell-type label (from celltypist / Azimuth / CoDi) to a coarse
    immune lineage. Unknown / NaN -> 'Other'."""
    if label is None:
        return "Other"
    s = str(label).lower().strip()
    if s in ("nan", "", "unassigned", "unknown"):
        return "Other"
    for lineage, pattern in _LINEAGE_RULES:
        if pattern.search(s):
            return lineage
    return "Other"


def annotate_celltypist(adata, debug: bool = False):
    """Run celltypist (majority voting per Leiden cluster) and write
    `obs['cell_type']` (fine label) and `obs['lineage']` (coarse)."""
    import celltypist
    from celltypist import models

    log = get_logger()
    require("leiden" in adata.obs.columns,
            "No 'leiden' column -- run Stage 2 (integration) first.")

    # Model is cached after first download; needs network only the first time.
    models.download_models(model=["Immune_All_Low.pkl"], force_update=False)
    pred = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True)

    labels = pred.predicted_labels
    adata.obs["cell_type"] = labels["majority_voting"].astype(str).values
    adata.obs["lineage"] = [to_lineage(x) for x in adata.obs["cell_type"]]

    counts = adata.obs["lineage"].value_counts()
    log.info(f"celltypist annotation: {adata.obs['cell_type'].nunique()} fine labels -> "
             f"{counts.shape[0]} lineages")
    log_table("Lineage counts (our annotation)", ["lineage", "n_cells"],
              list(counts.items()), debug_only=False)
    return adata


def load_codi_labels(adata):
    """Align the per-sample CoDi annotations (data/raw/*_CoDi_KLD.csv) to the
    cells of this AnnData. obs_names are 'barcode-<sample>'; the CoDi files are
    keyed by the bare barcode, per sample -- so strip the sample suffix and
    look up the matching file."""
    import pandas as pd
    log = get_logger()
    out = pd.Series(index=adata.obs_names, dtype=object)
    samples = adata.obs[cfg.BATCH_KEY].astype(str).values
    for label, fname in cfg.SAMPLES.items():
        codi_path = cfg.RAW_DIR / (fname.replace(".h5ad", "") + "_CoDi_KLD.csv")
        if not codi_path.exists():
            log.warning(f"CoDi reference missing for {label}: {codi_path.name}")
            continue
        ref = pd.read_csv(codi_path).set_index("cell_id")["CoDi"]
        mask = samples == label
        suffix = "-" + label
        barcodes = [n[:-len(suffix)] if n.endswith(suffix) else n
                    for n in adata.obs_names[mask]]
        out[mask] = ref.reindex(barcodes).values
    return out


def three_way_agreement(ours, azimuth, codi):
    """Pure function: given three label arrays (already coarse lineages),
    return overall pairwise + 3-way agreement fractions. Used both for the
    --debug report and the saved cross-check table."""
    ours = np.asarray(ours, dtype=object)
    azimuth = np.asarray(azimuth, dtype=object)
    codi = np.asarray(codi, dtype=object)
    n = len(ours)
    return {
        "n_cells": n,
        "ours_vs_azimuth": float(np.mean(ours == azimuth)) if n else 0.0,
        "ours_vs_codi": float(np.mean(ours == codi)) if n else 0.0,
        "azimuth_vs_codi": float(np.mean(azimuth == codi)) if n else 0.0,
        "all_three": float(np.mean((ours == azimuth) & (ours == codi))) if n else 0.0,
    }


def cross_check(adata):
    """Compare our lineage against the two independent references (Azimuth's
    predicted.celltype and CoDi), harmonized to lineages. Logs an agreement
    summary and a per-lineage breakdown; returns a per-lineage DataFrame."""
    import pandas as pd
    log = get_logger()

    ours = adata.obs["lineage"].astype(str).values
    az_raw = adata.obs["predicted.celltype"] if "predicted.celltype" in adata.obs else None
    azimuth = np.array([to_lineage(x) for x in az_raw]) if az_raw is not None \
        else np.array(["Other"] * adata.n_obs)
    codi_raw = load_codi_labels(adata)
    codi = np.array([to_lineage(x) for x in codi_raw])

    adata.obs["lineage_azimuth"] = azimuth
    adata.obs["lineage_codi"] = codi

    summary = three_way_agreement(ours, azimuth, codi)
    log_table("Annotation agreement (lineage-level)",
              ["comparison", "agreement"],
              [["ours vs Azimuth", f"{summary['ours_vs_azimuth']:.1%}"],
               ["ours vs CoDi", f"{summary['ours_vs_codi']:.1%}"],
               ["Azimuth vs CoDi", f"{summary['azimuth_vs_codi']:.1%}"],
               ["all three agree", f"{summary['all_three']:.1%}"]],
              debug_only=False)

    # Per-lineage: how often the references back up each of our calls.
    rows = []
    for lin in sorted(set(ours)):
        m = ours == lin
        rows.append([lin, int(m.sum()),
                     f"{np.mean(azimuth[m] == lin):.1%}",
                     f"{np.mean(codi[m] == lin):.1%}"])
    df = pd.DataFrame(rows, columns=["our_lineage", "n_cells", "azimuth_backs", "codi_backs"])
    log_table("Per-lineage reference support", list(df.columns), df.values.tolist(),
              debug_only=False)
    return df


def save_annotation_figures(adata, suffix: str = ""):
    """UMAP coloured by lineage + marker-gene dotplot, to results/figures/."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc
    log = get_logger()

    # UMAP by lineage
    coords = adata.obsm["X_umap"]
    labels = adata.obs["lineage"].astype(str).values
    fig, ax = plt.subplots(figsize=(8, 6))
    for g in sorted(set(labels)):
        m = labels == g
        ax.scatter(coords[m, 0], coords[m, 1], s=4, alpha=0.5, label=g, edgecolor="none")
    ax.set_title("Cell-type lineage (our annotation)")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.legend(markerscale=3, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(cfg.FIG_DIR / f"03_umap_lineage{suffix}.png", dpi=120)
    plt.close(fig)

    # Marker dotplot (only genes actually present)
    markers = {k: [g for g in v if g in adata.var_names] for k, v in cfg.MARKER_GENES.items()}
    markers = {k: v for k, v in markers.items() if v}
    sc.pl.dotplot(adata, markers, groupby="leiden", show=False)
    fig = plt.gcf()
    fig.savefig(cfg.FIG_DIR / f"03_marker_dotplot{suffix}.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    log.info(f"saved annotation figures -> 03_umap_lineage{suffix}.png, "
             f"03_marker_dotplot{suffix}.png")


def run(adata, debug: bool = False, smoke: bool = False):
    """Run Stage 3 on a Stage-2 (clustered) checkpoint AnnData."""
    describe_adata(adata, "annotation:input")

    adata = annotate_celltypist(adata, debug=debug)
    df = cross_check(adata)
    df.to_csv(cfg.TABLE_DIR / f"03_annotation_agreement{'_smoke' if smoke else ''}.csv",
              index=False)
    save_annotation_figures(adata, suffix="_smoke" if smoke else "")

    describe_adata(adata, "annotation:output")
    return adata
