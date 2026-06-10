"""
Shared helper functions. Keep notebooks thin by putting reusable logic here.
"""
import scanpy as sc
import anndata as ad
from . import config as cfg


def setup_scanpy():
    """Apply consistent scanpy/figure settings. Call once at the top of a notebook."""
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=80, dpi_save=150, frameon=False)
    sc.settings.figdir = cfg.FIG_DIR
    sc.settings.seed = cfg.RANDOM_SEED


def load_all_samples():
    """
    Read the four raw .h5ad files and concatenate them into one AnnData object.

    Adds an .obs['sample'] column so we can tell the four conditions apart and
    use it later as the batch key. Cell barcodes are made unique across samples.
    """
    adatas = {}
    for label, fname in cfg.SAMPLES.items():
        path = cfg.RAW_DIR / fname
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {path}. Download from Zenodo and update SAMPLES in config.py."
            )
        a = sc.read_h5ad(path)
        a.obs[cfg.BATCH_KEY] = label
        adatas[label] = a

    combined = ad.concat(adatas, label=cfg.BATCH_KEY, index_unique="-")
    combined.obs_names_make_unique()
    return combined


def save_checkpoint(adata, name):
    """Write an intermediate AnnData to data/processed for the next notebook."""
    out = cfg.PROCESSED_DIR / f"{name}.h5ad"
    adata.write_h5ad(out)
    print(f"saved -> {out}  ({adata.n_obs} cells x {adata.n_vars} genes)")
    return out


def load_checkpoint(name):
    """Read an intermediate AnnData produced by a previous notebook."""
    return sc.read_h5ad(cfg.PROCESSED_DIR / f"{name}.h5ad")
