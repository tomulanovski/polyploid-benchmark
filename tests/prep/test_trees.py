import pytest
from ete3 import Tree

from polybench.prep.trees import (
    clean_gene_trees,
    make_unique_labels,
    relabel_leaves,
    restore_labels,
    split_newick,
)


def test_substring_labels_roundtrip():
    taxa = ["Qlo", "Qlo03", "Abc"]
    mapping = make_unique_labels(taxa)
    unique = [mapping.get(t, t) for t in taxa]
    assert len(set(unique)) == 3
    assert not any(a in b for a in unique for b in unique if a != b)
    newick = "((Qlo:1,Qlo03:1):1,(Abc:1,Qlo:1):2);"
    assert restore_labels(relabel_leaves(newick, mapping), mapping) == newick


def test_unique_labels_untouched_when_no_substrings():
    assert make_unique_labels(["Aaa", "Bbb"]) == {}


def test_unique_label_avoids_existing_label():
    mapping = make_unique_labels(["Qlo", "QloX", "Qlo03"])
    unique = [mapping.get(t, t) for t in ["Qlo", "QloX", "Qlo03"]]
    assert not any(a in b for a in unique for b in unique if a != b)


def test_restore_only_whole_labels():
    mapping = {"Qlo": "QloX"}
    assert restore_labels("(QloX:1,QloXY:1);", mapping) == "(Qlo:1,QloXY:1);"


def test_strip_support_values():
    out = clean_gene_trees(["((A:1,B:1)100:1,(C:1,D:1)0.95:1)90:0;"])
    assert out == ["((A:1,B:1):1,(C:1,D:1):1);"]


def test_strip_internal_labels():
    out = clean_gene_trees(["((A:1,B:1)n110:1,(C:1,D:1)r98:1);"])
    assert out == ["((A:1,B:1):1,(C:1,D:1):1);"]


def test_resolve_polytomies_seeded():
    line = "(A:1,B:1,C:1,D:1,E:1);"
    first = clean_gene_trees([line], seed=3)
    assert first == clean_gene_trees([line], seed=3)
    tree = Tree(first[0], format=5)
    assert all(len(n.children) == 2 for n in tree.traverse() if not n.is_leaf())
    assert {leaf.name for leaf in tree} == set("ABCDE")


def test_bad_tree_fails_loud():
    with pytest.raises(ValueError):
        clean_gene_trees(["((A,B);"])


def test_split_newick_multiline():
    assert split_newick("(a,\n b);(c,d);\n\n") == ["(a,\n b);", "(c,d);"]
    assert split_newick("  \n") == []
