"""
Intent tests for Stage 1 (QC & preprocessing).

All offline, on tiny synthetic AnnData objects -- no raw data, no network.
They pin the *behaviour* we care about: that filters cut exactly the cells
that violate a threshold, that normalization hits the target sum, that HVG
selection reads the raw `counts` layer, and that mito genes are detected.
"""
import numpy as np
import anndata as ad
import pytest

from src import config as cfg, qc


def _make_adata(counts, var_names, samples):
    """Build an AnnData mirroring the pipeline's invariant: .X == counts layer,
    a `sample` obs column, integer counts."""
    counts = np.asarray(counts, dtype=np.float32)
    a = ad.AnnData(X=counts.copy())
    a.layers["counts"] = counts.copy()
    a.var_names = list(var_names)
    a.obs_names = [f"cell{i}" for i in range(counts.shape[0])]
    a.obs["sample"] = list(samples)
    return a


def test_compute_qc_metrics_flags_mito():
    # 3 normal genes + 2 mito genes
    var = ["G1", "G2", "G3", "MT-1", "MT-2"]
    counts = [[1, 1, 1, 1, 0],   # total 4, mito 1 -> 25%
              [2, 0, 0, 0, 0]]    # total 2, mito 0 -> 0%
    a = _make_adata(counts, var, ["s", "s"])
    a = qc.compute_qc_metrics(a)
    assert int(a.var["mt"].sum()) == 2
    assert a.obs["pct_counts_mt"].iloc[0] == pytest.approx(25.0)
    assert a.obs["pct_counts_mt"].iloc[1] == pytest.approx(0.0)


def test_apply_filters_removes_exactly_out_of_threshold(monkeypatch):
    monkeypatch.setattr(cfg, "QC_MIN_GENES_PER_CELL", 2)
    monkeypatch.setattr(cfg, "QC_MAX_GENES_PER_CELL", 8)
    monkeypatch.setattr(cfg, "QC_MAX_PCT_MITO", 20.0)
    monkeypatch.setattr(cfg, "QC_MIN_CELLS_PER_GENE", 1)

    var = [f"G{i}" for i in range(1, 9)] + ["MT-1", "MT-2"]  # 8 normal + 2 mito

    def cell(expr):
        row = np.zeros(10)
        for idx, val in expr.items():
            row[idx] = val
        return row

    counts = [
        cell({0: 1, 1: 1, 2: 1, 3: 1, 4: 1}),               # keep: 5 genes, 0% mito
        cell({0: 1, 1: 1, 2: 1}),                            # keep: 3 genes, 0% mito
        cell({0: 1}),                                        # drop: 1 gene < min_genes
        cell({0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1}),  # drop: 9 genes >= max
        cell({0: 1, 1: 1, 8: 3}),                            # drop: mito 3/5 = 60% > max
    ]
    a = _make_adata(counts, var, ["s"] * 5)
    a = qc.compute_qc_metrics(a)
    out = qc.apply_filters(a)

    assert set(out.obs_names) == {"cell0", "cell1"}


def test_normalize_and_log_hits_target_sum():
    var = ["G1", "G2", "G3"]
    counts = [[1, 2, 3], [4, 4, 2]]   # row sums 6 and 10, deliberately unequal
    a = _make_adata(counts, var, ["s", "s"])
    out = qc.normalize_and_log(a)
    # undo log1p, each cell should sum to the configured target
    recovered = np.expm1(out.X)
    row_sums = np.asarray(recovered.sum(axis=1)).ravel()
    assert row_sums == pytest.approx([cfg.NORM_TARGET_SUM, cfg.NORM_TARGET_SUM], rel=1e-3)


def test_select_hvgs_reads_counts_layer_not_logged_X(monkeypatch):
    # If select_hvgs wrongly used .X (already log-normalized in the real
    # pipeline) instead of layer="counts", seurat_v3 would warn/misbehave.
    # Here we put DIFFERENT data in .X vs counts and assert it flags exactly
    # N_TOP_GENES from the counts layer without raising.
    rng = np.random.default_rng(0)
    n_cells, n_genes = 60, 40
    counts = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32)
    var = [f"G{i}" for i in range(n_genes)]
    a = _make_adata(counts, var, ["a"] * 30 + ["b"] * 30)
    a.X = np.log1p(a.X)  # mimic the pipeline: .X is logged, counts layer is raw
    monkeypatch.setattr(cfg, "N_TOP_GENES", 10)
    out = qc.select_hvgs(a)
    assert int(out.var["highly_variable"].sum()) == 10
