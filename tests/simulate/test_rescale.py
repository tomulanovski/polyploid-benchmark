import pytest
from ete3 import Tree

from polybench.simulate.rescale import rescale


def test_rescale_height():
    t = Tree(rescale("((A,B),(C,(D,E)),F);", height=10_000_000), format=1)
    for leaf in t.iter_leaves():
        assert t.get_distance(leaf) == pytest.approx(1e7, abs=1e-3)


def test_rescale_resolves_polytomies_and_keeps_labels():
    t = Tree(rescale("(A,B,C,D);", seed=1), format=1)
    assert sorted(t.get_leaf_names()) == ["A", "B", "C", "D"]
    assert all(len(n.children) == 2 for n in t.traverse() if not n.is_leaf())


def test_rescale_repeated_labels_kept():
    t = Tree(rescale("((A,B),(A,C));"), format=1)
    assert sorted(t.get_leaf_names()) == ["A", "A", "B", "C"]
