import itertools
import subprocess
import sys

from ete3 import Tree

from polybench.prep.copy_selection import reduce_tree, select_copies

T3 = "((a:1,a:1):5,(a:1,(a:1,a:1):1):5);"


def test_k1_keeps_farthest_from_root():
    t = Tree("((a:1,a:4):1,b:1);", format=1)
    sel = select_copies(t, "a", 1)
    assert len(sel) == 1 and sel[0].dist == 4


def test_k3_farthest_point():
    t = Tree(T3, format=1)
    leaves = t.get_leaves()
    sel = select_copies(t, "a", 3)
    assert len({id(n) for n in sel}) == 3

    def min_d(sub):
        return min(t.get_distance(x, y) for x, y in itertools.combinations(sub, 2))

    assert min_d(sel) == max(min_d(s) for s in itertools.combinations(leaves, 3))
    first = set(t.children[0].get_leaves())
    assert len([n for n in sel if n in first]) == 1


def test_selection_is_deterministic():
    picks = set()
    for _ in range(20):
        t = Tree(T3, format=1)
        leaves = t.get_leaves()
        picks.add(tuple(leaves.index(n) for n in select_copies(t, "a", 3)))
    assert len(picks) == 1


def test_fewer_copies_than_k_returns_all():
    t = Tree("(a:1,b:1);", format=1)
    assert len(select_copies(t, "a", 3)) == 1


def test_reduce_tree_prunes_and_keeps_lengths():
    nwk = "(((a:1,a:1):2,a:3):1,(b:1,c:2):1);"
    out = Tree(reduce_tree(nwk, {"a": 2, "b": 1, "c": 1}), format=1)
    assert sum(n.name == "a" for n in out) == 2
    assert len(out) == 4
    assert abs(out.get_distance("c") - 3.0) < 1e-9


def test_reduce_tree_keeps_all_when_ploidy_exceeds_copies():
    out = Tree(reduce_tree("((a:1,a:1):1,b:1);", {"a": 4, "b": 1}), format=1)
    assert len(out) == 3


def test_reduce_tree_genus_species_label():
    nwk = "((Ephedra_sinica:1,Ephedra_sinica:1):1,Ephedra_foeminea:1);"
    out = Tree(reduce_tree(nwk, {"Ephedra_sinica": 1, "Ephedra_foeminea": 1}), format=1)
    assert sorted(out.get_leaf_names()) == ["Ephedra_foeminea", "Ephedra_sinica"]


def test_cli_reduces_every_tree(tmp_path):
    (tmp_path / "trees.nwk").write_text("(((a:1,a:1):2,a:3):1,(b:1,c:2):1);\n((a:1,a:2):1,b:1);\n")
    (tmp_path / "ploidy.tsv").write_text("a\t1\nb\t1\nc\t1\n")
    args = [str(tmp_path / f) for f in ("trees.nwk", "ploidy.tsv", "out.nwk")]
    subprocess.run([sys.executable, "-m", "polybench.prep.copy_selection", *args], check=True)
    out = (tmp_path / "out.nwk").read_text().splitlines()
    assert len(out) == 2
    for line in out:
        assert sorted(Tree(line, format=1).get_leaf_names()).count("a") == 1
    (tmp_path / "ploidy.tsv").write_text("a\t1\nb\t1\n")
    r = subprocess.run([sys.executable, "-m", "polybench.prep.copy_selection", *args], capture_output=True, text=True)
    assert r.returncode != 0 and "taxa without a ploidy: ['c']" in r.stderr
