"""
Intent tests for Stage 6 (size-specific effects).

The set-membership logic is exactly the kind of thing that goes subtly wrong
(an off-by-one set operation, mixture-emergent accidentally including genes that
were already significant in a single size). So every one of the four categories
is pinned explicitly.
"""
from src import size_effects as se


def test_categorize_all_four_categories():
    g40 = ["A", "B", "C"]          # A,B,C significant at 40nm
    g200 = ["B", "C", "D"]         # B,C,D significant at 200nm
    gmix = ["C", "E", "F"]         # C,E,F significant in the mixture

    cats = se.categorize(g40, g200, gmix)

    assert cats["unique_40nm"] == ["A"]            # in 40nm only
    assert cats["unique_200nm"] == ["D"]           # in 200nm only
    assert cats["shared"] == ["B", "C"]            # in both single sizes
    # E and F are in the mixture and in NEITHER single size; C is in mixture but
    # also in both single sizes -> must NOT count as emergent.
    assert cats["mixture_emergent"] == ["E", "F"]


def test_mixture_emergent_excludes_single_size_genes():
    # every mixture gene is already significant in a single size -> nothing emergent
    cats = se.categorize(["A"], ["B"], ["A", "B"])
    assert cats["mixture_emergent"] == []


def test_empty_inputs_give_empty_categories():
    cats = se.categorize([], [], [])
    assert all(v == [] for v in cats.values())


def test_categories_are_disjoint_where_expected():
    cats = se.categorize(["A", "B"], ["B", "C"], ["D"])
    # unique_40nm and unique_200nm never overlap with shared
    assert set(cats["unique_40nm"]).isdisjoint(cats["shared"])
    assert set(cats["unique_200nm"]).isdisjoint(cats["shared"])
    assert cats["mixture_emergent"] == ["D"]
