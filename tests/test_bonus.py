"""Intent tests for the pure helpers behind the bonus analyses. Offline."""
import numpy as np
import pandas as pd
import pytest

from src import bonus


def test_expected_vs_observed_additive_null():
    exp, obs, resid = bonus.expected_vs_observed([1.0, 2.0], [1.0, 1.0], [2.0, 1.0])
    # additive expectation = 40nm + 200nm
    assert list(exp) == [2.0, 3.0]
    assert list(obs) == [2.0, 1.0]
    # gene 1 is exactly additive (resid 0); gene 2 is sub-additive (resid -2)
    assert list(resid) == [0.0, -2.0]


def test_disruption_magnitude_counts_and_sums_only_significant():
    df = pd.DataFrame({
        "names":          ["g1", "g2", "g3", "g4"],
        "logfoldchanges": [2.0,  -3.0, 0.5,  4.0],
        "pvals_adj":      [0.01, 0.001, 0.01, 0.50],
    })
    # g1, g2 significant; g3 fails LFC gate; g4 fails padj gate
    n_sig, total = bonus.disruption_magnitude(df)
    assert n_sig == 2
    assert total == pytest.approx(5.0)   # |2| + |-3|
