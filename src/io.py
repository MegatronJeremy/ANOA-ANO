"""
Loading the 4 raw samples and reading/writing per-stage checkpoints.
"""
import anndata as ad
import numpy as np
import scanpy as sc

from . import config as cfg
from .logging_utils import get_logger


def load_all_samples() -> ad.AnnData:
    """
    Read the four raw .h5ad files and concatenate them into one AnnData.

    Each source file is a Seurat object exported to AnnData: `.X` holds
    Seurat's own LogNormalize output and `.layers['counts']` holds the raw
    integer counts (see the note in config.py). We deliberately discard the
    Seurat-side `.X` and recompute everything from `.layers['counts']`, so
    every cell entering this pipeline goes through the *same* scanpy
    normalization regardless of what was done upstream.

    Adds `.obs['sample']` so the four conditions are distinguishable
    downstream and usable as the batch key for integration. Cell barcodes
    collide across samples (e.g. every sample has a cell "AAACCCAAGTGGACGT-1"),
    so `index_unique="-"` is required for ad.concat to keep them distinct.
    """
    log = get_logger()
    adatas = {}
    for label, fname in cfg.SAMPLES.items():
        path = cfg.RAW_DIR / fname
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {path}. Run src/download_data.py or check config.SAMPLES."
            )
        a = sc.read_h5ad(path)

        if "counts" not in a.layers:
            raise AssertionError(
                f"{fname}: expected a 'counts' layer with raw integer counts, found none. "
                f"Available layers: {list(a.layers.keys())}"
            )
        a.X = a.layers["counts"].copy()

        a.obs[cfg.BATCH_KEY] = label
        n_obs, n_vars = a.shape
        log.debug(f"loaded {fname}: {n_obs} cells x {n_vars} genes (label={label})")
        adatas[label] = a

    combined = ad.concat(adatas, label=cfg.BATCH_KEY, index_unique="-", join="inner")
    combined.obs_names_make_unique()
    combined.var_names_make_unique()

    log.debug(
        f"concatenated: {combined.n_obs} cells x {combined.n_vars} genes "
        f"(genes = intersection across all 4 samples)"
    )
    return combined


def subsample_for_smoke_test(adata: ad.AnnData, n_per_sample: int, seed: int) -> ad.AnnData:
    """
    Take a random subset of up to n_per_sample cells from each sample, so the
    whole pipeline can be exercised end-to-end in seconds.
    """
    log = get_logger()
    rng = np.random.default_rng(seed)
    keep_idx = []
    batch_values = adata.obs[cfg.BATCH_KEY].astype(str).values
    for label in sorted(set(batch_values)):
        sample_idx = np.where(batch_values == label)[0]
        n_take = min(n_per_sample, len(sample_idx))
        chosen = rng.choice(sample_idx, size=n_take, replace=False)
        keep_idx.append(chosen)
        log.debug(f"smoke-test subsample: {label} -> {n_take}/{len(sample_idx)} cells")
    keep_idx = np.concatenate(keep_idx)
    keep_idx.sort()
    out = adata[keep_idx].copy()
    log.info(f"smoke-test subsample: {adata.n_obs} -> {out.n_obs} cells total")
    return out


def save_checkpoint(adata: ad.AnnData, name: str):
    """Write an intermediate AnnData to data/processed/ for the next stage."""
    log = get_logger()
    out = cfg.PROCESSED_DIR / f"{name}.h5ad"
    adata.write_h5ad(out)
    log.info(f"saved checkpoint -> {out}  ({adata.n_obs} cells x {adata.n_vars} genes)")
    return out


def load_checkpoint(name: str) -> ad.AnnData:
    """Read an intermediate AnnData produced by a previous stage."""
    path = cfg.PROCESSED_DIR / f"{name}.h5ad"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing checkpoint {path}. Run the earlier stage(s) first."
        )
    return sc.read_h5ad(path)
