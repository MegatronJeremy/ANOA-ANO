"""
Central configuration for the pipeline.

Everything tunable lives here so thresholds have one documented home and
`run_pipeline.py` / the stage modules stay free of magic numbers.
"""
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths. PROJECT_ROOT is resolved relative to this file so it works no matter
# where the driver is launched from.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DIR       = PROJECT_ROOT / "data" / "raw"        # untouched downloaded .h5ad files
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"  # per-stage checkpoints
FIG_DIR       = PROJECT_ROOT / "results" / "figures"
TABLE_DIR     = PROJECT_ROOT / "results" / "tables"

for _d in (RAW_DIR, PROCESSED_DIR, FIG_DIR, TABLE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Samples. Friendly label -> raw filename, confirmed against the actual
# Zenodo record (10.5281/zenodo.15866724): the un-suffixed file is Sample 1,
# and the record's own description states 40 nm = Sample 1, 200 nm =
# Sample 2, mixture = Sample 3, control = Sample 4.
# ---------------------------------------------------------------------------
SAMPLES = {
    "PSNP_40nm":    "filtered_feature_bc_matrix.h5ad",          # Sample 1 - 40 nm
    "PSNP_200nm":   "filtered_feature_bc_matrix_Sample2.h5ad",  # Sample 2 - 200 nm
    "PSNP_mixture": "filtered_feature_bc_matrix_Sample3.h5ad",  # Sample 3 - 40+200 nm mix
    "control":      "filtered_feature_bc_matrix_Sample4.h5ad",  # Sample 4 - no exposure
}
CONTROL_LABEL = "control"   # every exposed sample is compared against this
BATCH_KEY     = "sample"    # .obs column distinguishing the four samples

# ---------------------------------------------------------------------------
# IMPORTANT data note (found while building Stage 1, 2026-06-16):
# These .h5ad files are NOT raw scanpy exports -- they are Seurat objects
# converted to AnnData. `.X` already holds Seurat's LogNormalize output with
# scale.factor=1000 (verified: sum(expm1(X[0])) == 1000, not scanpy's default
# 10000). `.layers['counts']` holds the untouched integer raw counts.
# `.var`/`.obs` also carry Seurat-side HVG flags (`vf_vst_counts_*`) and an
# Azimuth annotation already run in R (`obs['predicted.celltype']`), plus
# Seurat PCA/UMAP embeddings in `.obsm`.
#
# To keep this pipeline self-contained and avoid silently inheriting another
# tool's normalization/HVG/clustering choices, Stage 1 resets `.X` to
# `.layers['counts']` and recomputes normalization, log1p and HVGs from
# scratch with scanpy. The Seurat-derived `predicted.celltype` column is left
# in `.obs` untouched -- it becomes a useful independent cross-check against
# our own Stage 3 celltypist annotation, not an input to this pipeline.
#
# Per-sample gene sets also differ slightly (Seurat applied its own per-sample
# feature filtering before export): 22613 / 23206 / 21715 / 21961 genes for
# samples 1-4 respectively. ad.concat's default join="inner" handles this by
# keeping only the genes shared across all 4 -- intentional, not a bug.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# QC thresholds (Stage 1). Justified against the real distribution of all
# 34,078 cells (printed by `--stage qc --debug`; measured 2026-06-21):
#   n_genes_by_counts: min=203, p1=541, p50=2122, p95=4074, p99=6296, max=10956
#   total_counts:      min=499, p50=5468, p95=15409, p99=35629, max=215279
#   pct_counts_mt:     min=0.0, p50=4.0, p95=7.9, p99=17.3, max=89.6
# The data is clearly already pre-filtered upstream (no cell below 203 genes /
# 499 counts), so the lower gene floor is a safety net, not the active filter.
# ---------------------------------------------------------------------------
QC_MIN_GENES_PER_CELL = 200     # safety floor: min observed is 203, so this is a guard
                                 # against empty droplets, not the binding filter here
QC_MAX_GENES_PER_CELL = 6000    # ~p99 (6296): drops the top ~1% as putative doublets
QC_MAX_PCT_MITO       = 15.0    # between p95 (7.9) and p99 (17.3): drops ~2-3% dying/
                                 # stressed cells with high mito% without cutting the bulk
QC_MIN_CELLS_PER_GENE = 3       # standard: drop genes seen in < 3 cells (noise / near-absent)

# ---------------------------------------------------------------------------
# Normalization / HVG (Stage 1)
# ---------------------------------------------------------------------------
NORM_TARGET_SUM = 1e4          # scanpy convention (CP10K); independent of Seurat's CP1K in source X
N_TOP_GENES     = 2000         # highly variable genes to keep
HVG_FLAVOR      = "seurat_v3"  # operates on raw counts directly, robust across batches

# ---------------------------------------------------------------------------
# Integration / clustering (Stage 2)
# ---------------------------------------------------------------------------
N_PCS       = 30    # PCs feeding the neighbour graph; standard scanpy starting point for PBMC-scale data
N_NEIGHBORS = 15     # kNN graph size for UMAP / Leiden; scanpy default
LEIDEN_RES  = 1.0    # clustering resolution; chosen to land roughly in the ~8-12 PBMC cell type range

# Cluster sanity-check thresholds (--debug warnings, not hard failures -- a
# tiny or dominant cluster might be real biology, but it's worth flagging
# for a manual look rather than letting it pass silently).
CLUSTER_TINY_FRACTION     = 0.01   # warn if a cluster has < 1% of all cells
CLUSTER_DOMINANT_FRACTION = 0.60   # warn if a cluster has > 60% of all cells
CLUSTER_BATCH_PURITY_WARN = 0.90   # warn if a cluster is > 90% one sample (integration
                                    # likely failed -- clusters should mix samples, not be
                                    # one-sample blobs, unless a population is truly condition-specific)

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
RANDOM_SEED = 0   # set everywhere: numpy, scanpy, harmony, leiden, UMAP

# ---------------------------------------------------------------------------
# Canonical PBMC marker genes, used for manual sanity-checks in Stage 3.
# ---------------------------------------------------------------------------
MARKER_GENES = {
    "T cell (CD4)":    ["IL7R", "CD3D", "CD4"],
    "T cell (CD8)":    ["CD8A", "CD3D", "GZMK"],
    "NK cell":         ["GNLY", "NKG7", "KLRD1"],
    "B cell":          ["MS4A1", "CD79A", "CD79B"],
    "Monocyte CD14":   ["CD14", "LYZ", "S100A8"],
    "Monocyte FCGR3A": ["FCGR3A", "MS4A7"],
    "Dendritic":       ["FCER1A", "CST3"],
    "Platelet":        ["PPBP", "PF4"],
}

# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------
SMOKE_TEST_CELLS_PER_SAMPLE = 500
