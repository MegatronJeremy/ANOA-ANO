"""
Intent tests for Stage 4 (composition analysis). Offline + synthetic.
Pins the descriptive maths: proportions sum to 1 per sample, control's
fold-change vs itself is 0, and a known shift comes out with the right sign.
"""
import numpy as np
import anndata as ad
import pandas as pd
import pytest

from src import config as cfg, composition


def _adata(sample, lineage):
    n = len(sample)
    a = ad.AnnData(X=np.zeros((n, 2), dtype=np.float32))
    a.obs[cfg.BATCH_KEY] = list(sample)
    a.obs["lineage"] = list(lineage)
    return a


def test_proportions_sum_to_one_per_sample():
    a = _adata(
        sample=["control"] * 10 + ["PSNP_40nm"] * 10,
        lineage=(["T cell"] * 8 + ["B cell"] * 2) + (["T cell"] * 6 + ["B cell"] * 4),
    )
    _, proportions, _ = composition.composition_tables(a)
    assert proportions.sum(axis=1).values == pytest.approx([1.0, 1.0])


def test_control_log2fc_is_zero_and_shift_sign_correct():
    a = _adata(
        sample=["control"] * 10 + ["PSNP_40nm"] * 10,
        lineage=(["T cell"] * 8 + ["B cell"] * 2) + (["T cell"] * 6 + ["B cell"] * 4),
    )
    _, _, log2fc = composition.composition_tables(a)
    # control vs itself -> 0 everywhere
    assert np.allclose(log2fc.loc["control"].values, 0.0)
    # 40nm T cells dropped 0.8 -> 0.6 => negative log2fc; B cells rose => positive
    assert log2fc.loc["PSNP_40nm", "T cell"] == pytest.approx(np.log2(0.6 / 0.8))
    assert log2fc.loc["PSNP_40nm", "B cell"] > 0


def test_missing_control_raises():
    a = _adata(sample=["PSNP_40nm"] * 4, lineage=["T cell"] * 4)
    with pytest.raises(AssertionError):
        composition.composition_tables(a)
