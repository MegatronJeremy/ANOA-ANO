"""
Intent tests for Stage 3 (cell-type annotation).

Offline. The load-bearing logic here is the label harmonization (three
different vocabularies -> one coarse lineage) and the agreement computation;
celltypist itself is not re-tested.
"""
import numpy as np
import pytest

from src import annotation as an


@pytest.mark.parametrize("label,expected", [
    # Azimuth / CoDi vocabulary (style "CD4+ T cell")
    ("CD4+ T cell", "T cell"),
    ("Cytotoxic T cell", "T cell"),
    ("B cell", "B cell"),
    ("CD14+ monocyte", "Monocyte"),
    ("Natural killer cell", "NK cell"),
    ("Plasmacytoid dendritic cell", "Dendritic cell"),
    # celltypist vocabulary (fine immune labels, plural)
    ("Tcm/Naive helper T cells", "T cell"),
    ("Regulatory T cells", "T cell"),
    ("Classical monocytes", "Monocyte"),
    ("NK cells", "NK cell"),
    ("Naive B cells", "B cell"),
    ("pDC", "Dendritic cell"),
    # tricky: plasma cell is B lineage, but pDC must NOT fall into B
    ("Plasma cells", "B cell"),
    ("Megakaryocytes/platelets", "Platelet"),
    # unknown / missing
    ("Mast cells", "Other"),
    (None, "Other"),
    (float("nan"), "Other"),
])
def test_to_lineage(label, expected):
    assert an.to_lineage(label) == expected


def test_three_way_agreement_perfect_and_partial():
    ours = np.array(["T cell", "B cell", "NK cell", "Monocyte"])

    # all identical -> everything 100%
    s = an.three_way_agreement(ours, ours, ours)
    assert s["all_three"] == pytest.approx(1.0)
    assert s["ours_vs_azimuth"] == pytest.approx(1.0)

    # azimuth differs on the last cell, codi agrees with ours everywhere
    azimuth = np.array(["T cell", "B cell", "NK cell", "T cell"])
    s2 = an.three_way_agreement(ours, azimuth, ours)
    assert s2["ours_vs_azimuth"] == pytest.approx(0.75)
    assert s2["ours_vs_codi"] == pytest.approx(1.0)
    assert s2["all_three"] == pytest.approx(0.75)


def test_to_lineage_empty_is_other():
    assert an.to_lineage("") == "Other"
    assert an.to_lineage("unassigned") == "Other"
