from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from ete3 import Tree

from polybench.simulate import fractionation as fr
from polybench.simulate.simphy import relabel_for_simphy

NETWORKS = Path(__file__).resolve().parent.parent / "data" / "model_networks"


def test_species_tree_leaf():
    assert fr.species_tree_leaf("Lamiumconfertum1_1_0_0") == "Lamiumconfertum1_1"
    assert fr.species_tree_leaf("Lamiumalbum2_0_0") == "Lamiumalbum2"
    assert fr.species_tree_leaf("Galeopsisspeciosa_0_0") == "Galeopsisspeciosa"


def test_split_copy():
    assert fr.split_copy("Lamiumhybridum3_1") == ("Lamiumhybridum3", 1)
    assert fr.split_copy("Lamiumhybridum3_2") == ("Lamiumhybridum3", 2)
    assert fr.split_copy("Lamiumalbum2") == ("Lamiumalbum2", None)


def _sides(event):
    return {frozenset(s) for s in event["sides"]}


def test_diploid_only_no_events():
    assert fr.build_events("((A,X),(B,Y));") == []


def test_rejects_duplicate_leaf_labels():
    with pytest.raises(ValueError):
        fr.build_events("((Alchemillakiwuensis,Alchemillakiwuensis),Y);")


def test_autopolyploid_duplicated_clade():
    events = fr.build_events("(((A_1,B_1),(A_2,B_2)),C);")
    assert len(events) == 1
    assert _sides(events[0]) == {frozenset({"A_1", "B_1"}), frozenset({"A_2", "B_2"})}
    assert events[0]["autopolyploidy"] is True


def test_allopolyploid_copies_under_different_parents():
    events = fr.build_events("((A_1,X),(A_2,Y));")
    assert len(events) == 1
    assert _sides(events[0]) == {frozenset({"A_1"}), frozenset({"A_2"})}
    assert events[0]["autopolyploidy"] is False


def test_flat_triple_copy():
    events = fr.build_events("((X_1,X_2,X_3),Y);")
    assert len(events) == 1
    assert _sides(events[0]) == {frozenset({"X_1"}), frozenset({"X_2"}), frozenset({"X_3"})}


def test_nested_autopolyploid_outer_and_mirrored_inner():
    sides = [_sides(e) for e in fr.build_events("(((X_1,X_2),(X_3,X_4)),Y);")]
    assert {frozenset({"X_1", "X_2"}), frozenset({"X_3", "X_4"})} in sides
    assert {frozenset({"X_1"}), frozenset({"X_2"})} in sides
    assert {frozenset({"X_3"}), frozenset({"X_4"})} in sides


def test_bendiksby_six_binary_allo_events():
    nwk, _ = relabel_for_simphy((NETWORKS / "Bendiksby_2011.tre").read_text().strip())
    events = fr.build_events(nwk)
    assert len(events) == 6
    species = set()
    for e in events:
        assert len(e["sides"]) == 2 and e["autopolyploidy"] is False
        bases = {fr.split_copy(n)[0] for side in e["sides"] for n in side}
        assert len(bases) == 1
        species.add(bases.pop())
    assert species == {
        "Lamiumamplexicaulevarorientale2", "Lamiumconfertum1", "Lamiumgaleobdolonsubspargentatum",
        "Lamiumgaleobdolonsubspmontanum", "Lamiumhybridum1", "Lamiumhybridum3"}


AUTO_TREE = "(((A_1,B_1),(A_2,B_2)),C);"
ALLO_TREE = "((A_1,X),(A_2,Y));"
NESTED_TREE = "(((X_1,X_2),(X_3,X_4)),Y);"


def _present(newick):
    return {leaf.name for leaf in Tree(newick, format=1).iter_leaves()}


def test_q1_removes_nothing():
    events = fr.build_events(AUTO_TREE)
    for seed in range(20):
        assert fr.plan_removals(events, [1.0] * len(events), _present(AUTO_TREE), np.random.default_rng(seed)) == set()


def test_q0_collapses_binary_event_to_single_copy():
    events = fr.build_events(ALLO_TREE)
    removed = fr.plan_removals(events, [0.0] * len(events), _present(ALLO_TREE), np.random.default_rng(0))
    assert removed in ({"A_1"}, {"A_2"})


def test_q0_keeps_at_least_one_copy_per_species_nested():
    events = fr.build_events(NESTED_TREE)
    present = _present(NESTED_TREE)
    for seed in range(20):
        removed = fr.plan_removals(events, [0.0] * len(events), present, np.random.default_rng(seed))
        assert any(fr.split_copy(leaf)[0] == "X" for leaf in present - removed)


def test_nested_inner_applies_to_whichever_outer_block_survives():
    events = fr.build_events(NESTED_TREE)
    present = _present(NESTED_TREE)
    for seed in range(30):
        removed = fr.plan_removals(events, [1.0, 0.0], present, np.random.default_rng(seed))
        assert len({leaf for leaf in present - removed if fr.split_copy(leaf)[0] == "X"}) == 2


def test_unbiased_drop_is_roughly_fifty_fifty():
    events = fr.build_events(ALLO_TREE)
    present = _present(ALLO_TREE)
    rng = np.random.default_rng(123)
    n = 4000
    dropped = sum(fr.plan_removals(events, [0.0], present, rng) == {"A_1"} for _ in range(n))
    assert abs(dropped / n - 0.5) < 0.05


def test_realized_retention_tracks_q():
    events = fr.build_events(ALLO_TREE)
    present = _present(ALLO_TREE)
    rng = np.random.default_rng(99)
    n, q = 8000, 0.6
    both = sum(not fr.plan_removals(events, [q], present, rng) for _ in range(n))
    assert abs(both / n - q) < 0.03


GENE = "((A_1_0_0:1,A_2_0_0:1):1,(B_0_0:1,C_0_0:1):1);"


def test_prune_removes_matching_tips():
    tips = Tree(fr.prune_gene_tree(GENE, {"A_1"}), format=1).get_leaf_names()
    assert sorted(tips) == ["A_2_0_0", "B_0_0", "C_0_0"]


def test_prune_nothing_keeps_all():
    assert len(Tree(fr.prune_gene_tree(GENE, set()), format=1).get_leaf_names()) == 4


def test_prune_keeps_branch_lengths():
    out = Tree(fr.prune_gene_tree(GENE, {"A_1"}), format=1)
    assert all(leaf.dist > 0 for leaf in out.iter_leaves())


def test_prune_handles_duplicate_tip_names():
    g = "((A_1_0_0:1,A_2_0_0:1):1,(B_0_0:1,B_0_0:1):1);"
    tips = sorted(Tree(fr.prune_gene_tree(g, {"A_1"}), format=1).get_leaf_names())
    assert tips == ["A_2_0_0", "B_0_0", "B_0_0"]


SP = "((A_1:1,A_2:1):1,B:2);"
GT = "((A_1_0_0:1,A_2_0_0:1):1,B_0_0:2);"


def test_apply_fractionation_zero_collapses_polyploid():
    out = fr.apply_fractionation([GT] * 10, SP, retention=0.0, seed=0)
    for line in out:
        tips = Tree(line, format=1).get_leaf_names()
        assert len([t for t in tips if t.startswith("A_")]) == 1 and "B_0_0" in tips


def test_apply_fractionation_one_keeps_everything():
    out = fr.apply_fractionation([GT] * 10, SP, retention=1.0, seed=0)
    assert all(Counter(Tree(x, format=1).get_leaf_names()) == Counter(["A_1_0_0", "A_2_0_0", "B_0_0"]) for x in out)


def test_apply_fractionation_deterministic_and_seeded():
    a = fr.apply_fractionation([GT] * 30, SP, 0.5, seed=4)
    assert a == fr.apply_fractionation([GT] * 30, SP, 0.5, seed=4)
    assert a != fr.apply_fractionation([GT] * 30, SP, 0.5, seed=5)


def test_apply_fractionation_rejects_bad_retention():
    with pytest.raises(ValueError):
        fr.apply_fractionation([GT], SP, retention=1.5, seed=0)
