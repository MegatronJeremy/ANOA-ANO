"""
Central configuration for the project.

Everything tunable lives here so the notebooks stay clean and so that
"why did you pick this threshold" has a single, documented answer.
Import in any notebook with:  from src import config as cfg
"""
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths. PROJECT_ROOT is resolved relative to this file so it works no matter
# where the notebook is launched from.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DIR        = PROJECT_ROOT / "data" / "raw"        # untouched downloaded .h5ad files
PROCESSED_DIR  = PROJECT_ROOT / "data" / "processed"  # intermediate .h5ad checkpoints
FIG_DIR        = PROJECT_ROOT / "results" / "figures"
TABLE_DIR      = PROJECT_ROOT / "results" / "tables"

for _d in (RAW_DIR, PROCESSED_DIR, FIG_DIR, TABLE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Samples. Map raw filenames -> friendly labels used throughout the analysis.
# !! Adjust the filenames after you download and unzip the Zenodo archive. !!
# ---------------------------------------------------------------------------
SAMPLES = {
    "PSNP_40nm":     "sample1_40nm.h5ad",      # 40 nm particles
    "PSNP_200nm":    "sample2_200nm.h5ad",     # 200 nm particles
    "PSNP_mixture":  "sample3_mixture.h5ad",   # 40 nm + 200 nm
    "control":       "sample4_control.h5ad",   # no exposure (reference)
}
CONTROL_LABEL = "control"   # every exposed sample is compared against this

# ---------------------------------------------------------------------------
# QC thresholds (step 1). Treat these as STARTING points. The QC notebook
# plots the distributions first; you then come back and set numbers you can
# defend from the plots. Justifying these is worth marks.
# ---------------------------------------------------------------------------
QC_MIN_GENES_PER_CELL = 200      # below this -> likely empty droplet / debris
QC_MAX_GENES_PER_CELL = 6000     # above this -> likely doublet (two cells in one droplet)
QC_MAX_PCT_MITO       = 15.0     # high mito % -> dying / stressed cell
QC_MIN_CELLS_PER_GENE = 3        # drop genes seen in fewer than this many cells

# ---------------------------------------------------------------------------
# Analysis parameters (steps 2+).
# ---------------------------------------------------------------------------
N_TOP_GENES   = 2000     # highly variable genes to keep
N_PCS         = 30       # principal components for the neighbour graph
N_NEIGHBORS   = 15       # kNN graph size for UMAP / Leiden
LEIDEN_RES    = 1.0      # clustering resolution (higher = more clusters)
BATCH_KEY     = "sample" # .obs column distinguishing the four samples
RANDOM_SEED   = 0        # set everywhere for reproducibility

# ---------------------------------------------------------------------------
# Canonical PBMC marker genes for manual annotation sanity-checks (step 3).
# Used to eyeball whether automated labels make sense.
# ---------------------------------------------------------------------------
MARKER_GENES = {
    "T cell (CD4)":  ["IL7R", "CD3D", "CD4"],
    "T cell (CD8)":  ["CD8A", "CD3D", "GZMK"],
    "NK cell":       ["GNLY", "NKG7", "KLRD1"],
    "B cell":        ["MS4A1", "CD79A", "CD79B"],
    "Monocyte CD14": ["CD14", "LYZ", "S100A8"],
    "Monocyte FCGR3A": ["FCGR3A", "MS4A7"],
    "Dendritic":     ["FCER1A", "CST3"],
    "Platelet":      ["PPBP", "PF4"],
}
