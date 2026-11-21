import pytest

from polybench.measures.compare_reticulations import pairwise_compare
from polybench.measures.reticulate_tree import ReticulateTree

A = "(((q,p),(r,p)),(s,(t,u)));"
B = "(((q,p),(s,p)),(r,(t,u)));"
C = "(((q,u),(r,u)),(s,(t,p)));"
D = "(((q,p),(r,p)),((s,u),(t,u)));"


@pytest.mark.parametrize("newick, desc, sis, bias, mu", [
    (C, 1.0, 0.0, 0, 1.0),
    (B, 0.0, 0.5, 0, 0.43),
    (D, 0.5, 0.5, 1, 0.38),
])
def test_table3(newick, desc, sis, bias, mu):
    m = pairwise_compare(ReticulateTree(A, refold=True), ReticulateTree(newick, refold=True))
    assert m["ret_leaf_jaccard"]["dist"] == pytest.approx(desc)
    assert m["ret_sisters_jaccard"]["dist"] == pytest.approx(sis)
    assert m["num_rets_bias"] == bias
    # manuscript reports values rounded half-up to 2 decimals
    assert round(m["mu_distance"] + 1e-9, 2) == mu
