import networkx as nx

from polybench.evaluate.compare import compare_files, load_network
from polybench.measures.reticulate_tree import ReticulateTree


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def test_load_network_reads_utf8_regardless_of_locale(tmp_path):
    # a ground truth with an accented taxon must load, not raise
    rt = load_network(write(tmp_path, "gt.tre", "((A,jalapaënsis),(B,jalapaënsis));"))
    assert "jalapaënsis" in rt._mu_taxa()


def test_load_network_refolds_a_network_given_as_enewick(tmp_path):
    # MPAllopp writes extended Newick, it must be refolded or it is compared at another resolution
    resolved = ReticulateTree(_resolved_hexaploid(), refold=False)
    loaded = load_network(write(tmp_path, "inf.tre", resolved.to_enewick()))
    assert sorted((d for _, d in loaded.dag.in_degree() if d >= 2), reverse=True) == [3]


def test_load_failure_carries_its_reason(tmp_path):
    gt = write(tmp_path, "gt.tre", "this is not newick")
    inf = write(tmp_path, "inf.tre", "((A,B),C);")
    try:
        compare_files(gt, inf, partial_match=False)
    except RuntimeError as e:
        msg = str(e)
    else:
        raise AssertionError("expected a load failure")
    assert msg.startswith(f"failed to load {gt}: ")
    assert len(msg) > len(f"failed to load {gt}: ")


def _resolved_hexaploid():
    G = nx.DiGraph()
    G.add_edges_from([('root', 'nAB'), ('root', 'pC'), ('nAB', 'pA'), ('nAB', 'pB'),
                      ('pA', 'A'), ('pB', 'B'), ('pC', 'C'),
                      ('pA', 'r1'), ('pB', 'r1'), ('r1', 'r2'), ('pC', 'r2'),
                      ('r2', 'X')])
    for leaf in ['A', 'B', 'C', 'X']:
        G.nodes[leaf]['label'] = leaf
    return G
