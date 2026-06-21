"""
Intent tests for Stage 2 (integration & clustering).

Offline + synthetic. The heavy parts (PCA/Harmony/UMAP/Leiden) are scanpy/
harmonypy code we don't re-test here; what we pin is OUR logic: the per-cluster
sample-composition computation and the stage-boundary guards.
"""
import numpy as np
import anndata as ad
import pandas as pd
import pytest

from src import config as cfg, integration


def _adata_with_obs(leiden, sample):
    n = len(leiden)
    a = ad.AnnData(X=np.zeros((n, 2), dtype=np.float32))
    a.obs["leiden"] = pd.Categorical([str(x) for x in leiden])
    a.obs[cfg.BATCH_KEY] = list(sample)
    return a


def test_composition_fractions_sum_to_one_per_cluster():
    # cluster 0: mixed across 2 samples; cluster 1: single sample
    a = _adata_with_obs(
        leiden=[0, 0, 0, 0, 1, 1],
        sample=["s1", "s2", "s1", "s2", "s1", "s1"],
    )
    frac = integration.cluster_sample_composition(a)
    row_sums = frac.sum(axis=1).values
    assert row_sums == pytest.approx([1.0, 1.0])


def test_composition_single_sample_cluster_is_pure():
    a = _adata_with_obs(
        leiden=[0, 0, 1, 1],
        sample=["s1", "s2", "s2", "s2"],
    )
    frac = integration.cluster_sample_composition(a)
    # cluster 1 is 100% s2 -> max fraction 1.0 (the case the purity warning targets)
    assert frac.loc["1"].max() == pytest.approx(1.0)
    # cluster 0 is evenly split -> max fraction 0.5
    assert frac.loc["0"].max() == pytest.approx(0.5)


def test_run_pca_requires_hvg_from_stage1():
    a = ad.AnnData(X=np.random.default_rng(0).poisson(1.0, size=(20, 10)).astype(np.float32))
    # no 'highly_variable' column -> stage-boundary guard must fire
    with pytest.raises(AssertionError):
        integration.run_pca(a)
