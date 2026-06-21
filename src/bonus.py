"""
Bonus analyses (deliverable: 3-5 additional insights).

Five extra analyses, each building on a finding from the main pipeline:

  1. stress / inflammation gene-module scoring (translates DE into programs)
  2. mixture additivity (is the mixture = 40nm + 200nm, or not?)
  3. robustness of clustering to the Leiden resolution
  4. lightweight ligand-receptor communication shift (exposed vs control)
  5. dose-response: transcriptional disruption magnitude vs particle size

Each writes a 07_* table + figure to results/ and is independently guarded, so
one failing analysis never sinks the rest. Registered as a single `--stage bonus`.
"""
import numpy as np

from . import config as cfg
from .logging_utils import get_logger, log_table
from .differential_expression import significant_genes
from .size_effects import load_de_results, _safe


# Curated, canonical PBMC-expressed signatures (score_genes uses the intersection).
GENE_MODULES = {
    "oxidative_stress": ["HMOX1", "NQO1", "GCLM", "GCLC", "TXN", "SOD1", "SOD2",
                         "PRDX1", "GPX1", "FTL", "FTH1"],
    "inflammation_nfkb": ["NFKB1", "NFKBIA", "TNF", "IL1B", "IL6", "CXCL8", "CCL3",
                          "CCL4", "PTGS2", "TNFAIP3"],
    "interferon": ["ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1", "STAT1", "IRF7",
                   "ISG20", "IFITM3"],
    "heat_shock": ["HSPA1A", "HSPA1B", "HSPB1", "DNAJB1", "HSPH1", "BAG3"],
}

# Ligand-receptor pairs relevant to monocyte / inflammatory signalling.
LR_PAIRS = [
    ("CCL3", "CCR1"), ("CCL4", "CCR5"), ("CXCL8", "CXCR1"), ("CXCL8", "CXCR2"),
    ("IL1B", "IL1R1"), ("TNF", "TNFRSF1A"), ("TNF", "TNFRSF1B"), ("CCL2", "CCR2"),
    ("IL15", "IL2RB"), ("CD40LG", "CD40"),
]


def _exposed():
    return [s for s in cfg.SAMPLES if s != cfg.CONTROL_LABEL]


# --- 1. gene-module scoring ------------------------------------------------
def module_scoring(adata, smoke=False):
    import scanpy as sc
    import pandas as pd
    log = get_logger()
    sfx = "_smoke" if smoke else ""
    for name, genes in GENE_MODULES.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        sc.tl.score_genes(adata, present, score_name=f"score_{name}", ctrl_size=50)
    score_cols = [f"score_{n}" for n in GENE_MODULES if f"score_{n}" in adata.obs.columns]
    tbl = adata.obs.groupby([cfg.BATCH_KEY, "lineage"])[score_cols].mean().reset_index()
    tbl.to_csv(cfg.TABLE_DIR / f"07_module_scores{sfx}.csv", index=False)

    # heatmap: module score per sample (averaged over cells), control first
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    per_sample = adata.obs.groupby(cfg.BATCH_KEY)[score_cols].mean()
    order = [cfg.CONTROL_LABEL] + [s for s in per_sample.index if s != cfg.CONTROL_LABEL]
    per_sample = per_sample.loc[order]
    fig, ax = plt.subplots(figsize=(7, 4))
    im = ax.imshow(per_sample.values, cmap="RdBu_r", aspect="auto",
                   vmin=-np.abs(per_sample.values).max(), vmax=np.abs(per_sample.values).max())
    ax.set_xticks(range(len(score_cols))); ax.set_xticklabels([c.replace("score_", "") for c in score_cols], rotation=30, ha="right")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
    fig.colorbar(im, ax=ax, label="mean module score")
    ax.set_title("Stress / inflammation module scores per sample")
    fig.tight_layout(); fig.savefig(cfg.FIG_DIR / f"07_module_scores{sfx}.png", dpi=120); plt.close(fig)
    log.info(f"[bonus 1] module scoring -> 07_module_scores{sfx}.csv/.png")


# --- 2. mixture additivity -------------------------------------------------
def expected_vs_observed(lfc40, lfc200, lfcmix):
    """Pure: additive null on the log scale is lfc40 + lfc200. Returns
    (expected, observed, residual=observed-expected) arrays."""
    expected = np.asarray(lfc40) + np.asarray(lfc200)
    observed = np.asarray(lfcmix)
    return expected, observed, observed - expected


def mixture_additivity(smoke=False):
    import pandas as pd
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    log = get_logger()
    sfx = "_smoke" if smoke else ""
    de = load_de_results(smoke=smoke)
    c40 = next(c for c in cfg.SAMPLES if c.endswith("40nm"))
    c200 = next(c for c in cfg.SAMPLES if c.endswith("200nm"))
    cmix = next(c for c in cfg.SAMPLES if "mix" in c.lower())
    lineages = sorted({lin for (lin, _) in de})

    rows = []
    fig, axes = plt.subplots(1, max(1, len(lineages)), figsize=(4 * max(1, len(lineages)), 4), squeeze=False)
    for ax, lin in zip(axes[0], lineages):
        try:
            d40 = de[(lin, c40)].set_index("names")["logfoldchanges"]
            d200 = de[(lin, c200)].set_index("names")["logfoldchanges"]
            dmix = de[(lin, cmix)].set_index("names")["logfoldchanges"]
        except KeyError:
            continue
        common = d40.index.intersection(d200.index).intersection(dmix.index)
        exp, obs, resid = expected_vs_observed(d40[common], d200[common], dmix[common])
        # focus on genes responding in at least one single size
        responsive = (np.abs(d40[common]) > cfg.DE_LFC_THRESHOLD) | (np.abs(d200[common]) > cfg.DE_LFC_THRESHOLD)
        n_resp = int(responsive.sum())
        sub = int((np.abs(obs[responsive]) < np.abs(exp[responsive])).sum())
        rows.append([lin, n_resp, sub, n_resp - sub,
                     f"{sub / n_resp:.1%}" if n_resp else "n/a"])
        ax.scatter(exp[responsive], obs[responsive], s=5, alpha=0.4)
        lim = np.nanpercentile(np.abs(np.concatenate([exp[responsive], obs[responsive]])) if n_resp else [1], 99)
        ax.plot([-lim, lim], [-lim, lim], "r--", lw=0.7)
        ax.set_title(lin, fontsize=9); ax.set_xlabel("expected (40+200)"); ax.set_ylabel("observed (mix)")
    fig.suptitle("Mixture additivity: observed vs additive expectation")
    fig.tight_layout(); fig.savefig(cfg.FIG_DIR / f"07_mixture_additivity{sfx}.png", dpi=120); plt.close(fig)

    pd.DataFrame(rows, columns=["lineage", "n_responsive", "n_sub_additive", "n_supra_additive",
                                "pct_sub_additive"]).to_csv(
        cfg.TABLE_DIR / f"07_mixture_additivity{sfx}.csv", index=False)
    log.info(f"[bonus 2] mixture additivity -> 07_mixture_additivity{sfx}.csv/.png")


# --- 3. clustering robustness ----------------------------------------------
def clustering_robustness(adata, smoke=False):
    import scanpy as sc
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    log = get_logger()
    sfx = "_smoke" if smoke else ""
    if "X_pca_harmony" not in adata.obsm:
        log.warning("[bonus 3] no X_pca_harmony -- skipping robustness"); return
    sc.pp.neighbors(adata, n_neighbors=min(cfg.N_NEIGHBORS, adata.n_obs - 1),
                    use_rep="X_pca_harmony", random_state=cfg.RANDOM_SEED)
    base = adata.obs["leiden"].astype(str).values if "leiden" in adata.obs else None
    rows = []
    for res in [0.5, 1.0, 1.5, 2.0]:
        key = f"leiden_r{res}"
        with __import__("warnings").catch_warnings():
            __import__("warnings").simplefilter("ignore")
            sc.tl.leiden(adata, resolution=res, key_added=key, random_state=cfg.RANDOM_SEED, flavor="leidenalg")
        labels = adata.obs[key].astype(str).values
        ari = adjusted_rand_score(base, labels) if base is not None else np.nan
        rows.append([res, adata.obs[key].nunique(), round(ari, 3)])
    df = pd.DataFrame(rows, columns=["resolution", "n_clusters", "ARI_vs_base(res1.0)"])
    df.to_csv(cfg.TABLE_DIR / f"07_clustering_robustness{sfx}.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df["resolution"], df["n_clusters"], "o-", label="n_clusters")
    ax.set_xlabel("Leiden resolution"); ax.set_ylabel("n_clusters")
    ax2 = ax.twinx(); ax2.plot(df["resolution"], df["ARI_vs_base(res1.0)"], "s--", color="green", label="ARI vs base")
    ax2.set_ylabel("ARI vs resolution=1.0")
    ax.set_title("Clustering robustness to resolution")
    fig.tight_layout(); fig.savefig(cfg.FIG_DIR / f"07_clustering_robustness{sfx}.png", dpi=120); plt.close(fig)
    log_table("Clustering robustness", list(df.columns), df.values.tolist(), debug_only=False)
    log.info(f"[bonus 3] clustering robustness -> 07_clustering_robustness{sfx}.csv/.png")


# --- 4. ligand-receptor communication --------------------------------------
def ligand_receptor(adata, smoke=False):
    import pandas as pd
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    log = get_logger()
    sfx = "_smoke" if smoke else ""

    def mean_expr(genes, mask):
        present = [g for g in genes if g in adata.var_names]
        if not present or mask.sum() == 0:
            return 0.0
        sub = adata[mask, present].X
        return float(np.asarray(sub.mean()))

    rows = []
    for sample in adata.obs[cfg.BATCH_KEY].unique():
        m = (adata.obs[cfg.BATCH_KEY] == sample).values
        for lig, rec in LR_PAIRS:
            score = mean_expr([lig], m) * mean_expr([rec], m)
            rows.append([sample, f"{lig}-{rec}", score])
    df = pd.DataFrame(rows, columns=["sample", "pair", "comm_score"])
    wide = df.pivot(index="pair", columns="sample", values="comm_score")
    if cfg.CONTROL_LABEL in wide.columns:
        rel = np.log2(wide.div(wide[cfg.CONTROL_LABEL], axis=0).replace(0, np.nan))
        rel.to_csv(cfg.TABLE_DIR / f"07_ligand_receptor_log2fc{sfx}.csv")
    wide.to_csv(cfg.TABLE_DIR / f"07_ligand_receptor_scores{sfx}.csv")

    exposed = [s for s in _exposed() if s in wide.columns]
    if cfg.CONTROL_LABEL in wide.columns and exposed:
        fig, ax = plt.subplots(figsize=(8, 5))
        rel[exposed].plot(kind="barh", ax=ax)
        ax.axvline(0, color="grey", lw=0.6)
        ax.set_xlabel("log2 fold-change vs control"); ax.set_title("Ligand-receptor communication shift")
        ax.legend(fontsize=7, frameon=False)
        fig.tight_layout(); fig.savefig(cfg.FIG_DIR / f"07_ligand_receptor{sfx}.png", dpi=120); plt.close(fig)
    log.info(f"[bonus 4] ligand-receptor -> 07_ligand_receptor_*{sfx}.csv/.png")


# --- 5. dose-response ------------------------------------------------------
def disruption_magnitude(df):
    """Pure: from one DE table, (n_significant, total |log2FC| of significant)."""
    sig = df[(df["pvals_adj"] < cfg.DE_PADJ_THRESHOLD) & (df["logfoldchanges"].abs() > cfg.DE_LFC_THRESHOLD)]
    return len(sig), float(sig["logfoldchanges"].abs().sum())


def dose_response(smoke=False):
    import pandas as pd
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    log = get_logger()
    sfx = "_smoke" if smoke else ""
    de = load_de_results(smoke=smoke)
    rows = []
    for (lin, cond), df in de.items():
        n_sig, total = disruption_magnitude(df)
        rows.append([lin, cond, n_sig, round(total, 1)])
    tbl = pd.DataFrame(rows, columns=["lineage", "condition", "n_significant", "total_abs_log2fc"])
    tbl.to_csv(cfg.TABLE_DIR / f"07_dose_response{sfx}.csv", index=False)

    order = [c for c in (next(c for c in cfg.SAMPLES if c.endswith("40nm")),
                         next(c for c in cfg.SAMPLES if c.endswith("200nm")),
                         next(c for c in cfg.SAMPLES if "mix" in c.lower()))]
    pivot = tbl.pivot(index="condition", columns="lineage", values="n_significant").reindex(order)
    fig, ax = plt.subplots(figsize=(8, 5))
    pivot.plot(kind="bar", ax=ax)
    ax.set_ylabel("# significant DE genes"); ax.set_title("Transcriptional disruption by particle size")
    ax.tick_params(axis="x", rotation=15); ax.legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(cfg.FIG_DIR / f"07_dose_response{sfx}.png", dpi=120); plt.close(fig)
    log_table("Dose-response (# significant DE genes)", list(tbl.columns), tbl.values.tolist(), debug_only=False)
    log.info(f"[bonus 5] dose-response -> 07_dose_response{sfx}.csv/.png")


def run(adata, debug: bool = False, smoke: bool = False):
    """Run all five bonus analyses, each guarded independently."""
    log = get_logger()
    for i, fn in enumerate([
        lambda: module_scoring(adata, smoke=smoke),
        lambda: mixture_additivity(smoke=smoke),
        lambda: clustering_robustness(adata, smoke=smoke),
        lambda: ligand_receptor(adata, smoke=smoke),
        lambda: dose_response(smoke=smoke),
    ], 1):
        try:
            fn()
        except Exception as exc:
            log.warning(f"[bonus {i}] failed ({type(exc).__name__}: {exc}) -- continuing")
    return adata
