"""
Intent tests for Stage 5 (differential expression). Offline + synthetic.
The scanpy Wilcoxon call is not re-tested; what is pinned is OUR gating logic
(which genes count as significant) and the cross-cell-type recurrence check.
"""
import pandas as pd
import pytest

from src import differential_expression as de


def _df():
    return pd.DataFrame({
        "names":          ["g_up", "g_down", "g_small", "g_nsig"],
        "logfoldchanges": [2.0,    -2.0,     0.5,       3.0],
        "pvals_adj":      [0.01,   0.01,     0.001,     0.20],
    })


def test_significant_genes_up_down_any():
    df = _df()
    assert de.significant_genes(df, padj=0.05, lfc=1.0, direction="up") == ["g_up"]
    assert de.significant_genes(df, padj=0.05, lfc=1.0, direction="down") == ["g_down"]
    assert set(de.significant_genes(df, padj=0.05, lfc=1.0, direction="any")) == {"g_up", "g_down"}
    # g_small fails the LFC gate; g_nsig fails the padj gate
    assert "g_small" not in de.significant_genes(df, direction="any")
    assert "g_nsig" not in de.significant_genes(df, direction="any")


def test_recurrence_keeps_only_genes_in_multiple_comparisons():
    df1 = pd.DataFrame({"names": ["g1", "g2"], "logfoldchanges": [2.0, 2.0], "pvals_adj": [0.01, 0.01]})
    df2 = pd.DataFrame({"names": ["g1", "g3"], "logfoldchanges": [2.0, 2.0], "pvals_adj": [0.01, 0.01]})
    out = de.recurrence_check({("T cell", "PSNP_40nm"): df1, ("B cell", "PSNP_40nm"): df2})
    # g1 is significant in both -> kept (n=2); g2/g3 appear once -> dropped
    assert out["gene"].tolist() == ["g1"]
    assert out.loc[out["gene"] == "g1", "n_comparisons_significant"].iloc[0] == 2
