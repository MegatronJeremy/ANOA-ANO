"""
Stage 5: Differential expression + pathway enrichment.

Ported from legacy/notebooks/05_differential_expression.ipynb. For each major
cell type (lineage) and each exposed sample vs control, find genes that change.

METHOD (one donor, NO biological replicates): cell-level Wilcoxon via scanpy
`rank_genes_groups`, NOT pseudobulk/DESeq2 (which estimate variance from
replicates that don't exist here). Results carry an explicit pseudoreplication
caveat -- cells are not independent units. There is no external gene-level
ground truth (CoDi_KLD is an annotation, not a DE reference), so DE is validated
INTERNALLY: genes recurring across cell types + pathway coherence.

Pathway enrichment uses Enrichr (needs network); it degrades gracefully to
"skipped" if offline rather than failing the stage.
"""
import numpy as np

from . import config as cfg
from .logging_utils import get_logger, describe_adata, log_table, require


def significant_genes(df, padj=None, lfc=None, direction="up"):
    """Pure helper: significant gene names from a scanpy DE dataframe
    (columns: names, logfoldchanges, pvals_adj). direction in {'up','down','any'}.
    Unit-tested -- it is the gate feeding enrichment and the recurrence check."""
    padj = cfg.DE_PADJ_THRESHOLD if padj is None else padj
    lfc = cfg.DE_LFC_THRESHOLD if lfc is None else lfc
    sig = df[df["pvals_adj"] < padj]
    if direction == "up":
        sig = sig[sig["logfoldchanges"] > lfc]
    elif direction == "down":
        sig = sig[sig["logfoldchanges"] < -lfc]
    else:
        sig = sig[sig["logfoldchanges"].abs() > lfc]
    return sig["names"].tolist()


def _safe(name):
    return str(name).replace(" ", "_").replace("/", "-")


def de_per_lineage(adata, smoke: bool = False):
    """Run Wilcoxon DE for every (lineage, exposed-vs-control) combo with enough
    cells in both groups. Returns {(lineage, cond): df}; saves one CSV each."""
    import scanpy as sc
    log = get_logger()
    require("lineage" in adata.obs.columns, "No 'lineage' column -- run Stage 3 first.")
    sfx = "_smoke" if smoke else ""
    exposed = [s for s in cfg.SAMPLES if s != cfg.CONTROL_LABEL]

    results, summary_rows = {}, []
    for lin in sorted(set(adata.obs["lineage"])):
        sub = adata[adata.obs["lineage"] == lin]
        for cond in exposed:
            s2 = sub[sub.obs[cfg.BATCH_KEY].isin([cond, cfg.CONTROL_LABEL])].copy()
            s2.obs[cfg.BATCH_KEY] = s2.obs[cfg.BATCH_KEY].astype(str).astype("category")
            n_exp = int((s2.obs[cfg.BATCH_KEY] == cond).sum())
            n_ctrl = int((s2.obs[cfg.BATCH_KEY] == cfg.CONTROL_LABEL).sum())
            if n_exp < cfg.DE_MIN_CELLS_PER_GROUP or n_ctrl < cfg.DE_MIN_CELLS_PER_GROUP:
                summary_rows.append([lin, cond, n_exp, n_ctrl, "skipped (too few cells)"])
                continue
            sc.tl.rank_genes_groups(s2, cfg.BATCH_KEY, groups=[cond],
                                    reference=cfg.CONTROL_LABEL, method="wilcoxon")
            df = sc.get.rank_genes_groups_df(s2, group=cond)
            results[(lin, cond)] = df
            df.to_csv(cfg.TABLE_DIR / f"05_DE_{_safe(lin)}_{cond}_vs_control{sfx}.csv", index=False)
            n_sig = len(significant_genes(df, direction="any"))
            summary_rows.append([lin, cond, n_exp, n_ctrl, f"{n_sig} sig genes"])

    log_table("DE per (lineage, condition)",
              ["lineage", "condition", "n_exposed", "n_control", "result"],
              summary_rows, debug_only=False)
    return results


def save_volcano(df, lineage, cond, suffix=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = df["logfoldchanges"].values
    y = -np.log10(df["pvals_adj"].clip(lower=1e-300).values)
    sig = (df["pvals_adj"] < cfg.DE_PADJ_THRESHOLD) & (df["logfoldchanges"].abs() > cfg.DE_LFC_THRESHOLD)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x[~sig.values], y[~sig.values], s=4, color="lightgrey", alpha=0.5)
    ax.scatter(x[sig.values], y[sig.values], s=6, color="crimson", alpha=0.7)
    ax.axvline(cfg.DE_LFC_THRESHOLD, ls="--", lw=0.6, color="grey")
    ax.axvline(-cfg.DE_LFC_THRESHOLD, ls="--", lw=0.6, color="grey")
    ax.set_xlabel("log2 fold-change"); ax.set_ylabel("-log10 adj. p")
    ax.set_title(f"{lineage} -- {cond} vs control")
    fig.tight_layout()
    fig.savefig(cfg.FIG_DIR / f"05_volcano_{_safe(lineage)}_{cond}{suffix}.png", dpi=110)
    plt.close(fig)


def recurrence_check(results):
    """Internal validation: genes significant across MULTIPLE (lineage, cond)
    combos. A coherent biological response (stress/inflammation) should recur,
    not appear in one comparison only. Returns a DataFrame, most-recurrent first."""
    import pandas as pd
    from collections import Counter
    counter = Counter()
    for (lin, cond), df in results.items():
        for g in significant_genes(df, direction="any"):
            counter[g] += 1
    rows = [[g, n] for g, n in counter.most_common() if n >= 2]
    return pd.DataFrame(rows, columns=["gene", "n_comparisons_significant"])


def run_enrichment(results, smoke: bool = False):
    """Best-effort Enrichr enrichment on up-regulated genes per combo. Degrades
    to 'skipped' (one warning) if offline -- never fails the stage."""
    import pandas as pd
    log = get_logger()
    sfx = "_smoke" if smoke else ""
    all_rows, network_ok = [], True
    for (lin, cond), df in results.items():
        if not network_ok:
            break
        genes = significant_genes(df, direction="up")
        if len(genes) < 10:
            continue
        try:
            import gseapy as gp
            enr = gp.enrichr(gene_list=genes, gene_sets=cfg.ENRICHR_GENE_SETS, outdir=None)
            res = enr.results.sort_values("Adjusted P-value").head(10)
            res.insert(0, "condition", cond); res.insert(0, "lineage", lin)
            all_rows.append(res)
        except Exception as exc:
            log.warning(f"pathway enrichment skipped (Enrichr unavailable: {type(exc).__name__}). "
                        f"DE tables are unaffected.")
            network_ok = False
    if all_rows:
        out = pd.concat(all_rows, ignore_index=True)
        out.to_csv(cfg.TABLE_DIR / f"05_pathway_enrichment{sfx}.csv", index=False)
        log.info(f"saved pathway enrichment -> 05_pathway_enrichment{sfx}.csv ({len(out)} rows)")
    else:
        log.info("no pathway enrichment table written (no combo had >=10 up genes, or offline)")


def run(adata, debug: bool = False, smoke: bool = False):
    """Run Stage 5 on a Stage-3 (annotated) checkpoint AnnData. Tables + figures
    only, no checkpoint."""
    describe_adata(adata, "de:input")
    sfx = "_smoke" if smoke else ""

    results = de_per_lineage(adata, smoke=smoke)
    for (lin, cond), df in results.items():
        save_volcano(df, lin, cond, suffix=sfx)

    recurrent = recurrence_check(results)
    recurrent.to_csv(cfg.TABLE_DIR / f"05_DE_recurrent_genes{sfx}.csv", index=False)
    get_logger().info(f"internal validation: {len(recurrent)} genes significant in >=2 comparisons")

    run_enrichment(results, smoke=smoke)
    return adata
