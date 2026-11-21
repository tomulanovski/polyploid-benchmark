import pickle
from pathlib import Path
from collections import Counter

import pytest
from ete3 import Tree

from polybench.simulate.simphy import (CONFIGS, prepare_species_tree, relabel_for_simphy, restore_simphy_labels,
                                       sample_subst_rate, simphy_command, write_replicates)


@pytest.fixture
def rates_pkl(tmp_path):
    p = tmp_path / "rates.pkl"
    p.write_bytes(pickle.dumps({"d1": [{"rate": 1e-9}, {"rate": 2e-9}], "d2": [{"rate": 3e-9}]}))
    return p


def test_simphy_command_rates(tmp_path, rates_pkl):
    c = CONFIGS["conf_dup_loss_high_10M_ne2M"]
    cmd = simphy_command(tmp_path / "s.nex", tmp_path / "out", c["ne"], c["dup"], c["loss"],
                         seed=7, subst_rates_pkl=rates_pkl)
    s = " ".join(cmd)
    assert "-lb f:1e-07 -ld f:1e-07 -sp f:2000000" in s
    assert "-rl f:1000" in s and "-rs 5" in s and "-cs 7" in s
    assert "-lt f:0 -lg f:0" in s and "-sg f:1" in s


def test_simphy_command_no_dup_loss(tmp_path, rates_pkl):
    cmd = simphy_command(tmp_path / "s.nex", tmp_path / "out", 200_000, 0, 0, seed=1, subst_rates_pkl=rates_pkl)
    assert "-lb f:0 -ld f:0 -sp f:200000" in " ".join(cmd)


def test_subst_rate_sampled_from_pool_deterministically(rates_pkl):
    assert sample_subst_rate(rates_pkl, 3) == sample_subst_rate(rates_pkl, 3)
    assert sample_subst_rate(rates_pkl, 3) in (1e-9, 2e-9, 3e-9)


def test_configs_grid():
    base = {k: v for k, v in CONFIGS.items() if "retention" not in v}
    assert len(base) == 12
    assert {(v["ne"], v["dup"]) for v in base.values()} == {
        (ne, d) for ne in (200_000, 1_000_000, 2_000_000) for d in (0, 1e-9, 1e-8, 1e-7)}
    assert all(v["dup"] == v["loss"] for v in CONFIGS.values())
    assert CONFIGS["conf_ils_high_10M"] == {"ne": 2_000_000, "dup": 0, "loss": 0}


def test_retention_configs():
    fix = {k: v for k, v in CONFIGS.items() if "retention" in v}
    assert {k: v["retention"] for k, v in fix.items()} == {
        "conf_dup_loss_medium_10M_ne1M_fix075": 0.75,
        "conf_dup_loss_medium_10M_ne1M_fix050": 0.50,
        "conf_dup_loss_medium_10M_ne1M_fix025": 0.25}
    for v in fix.values():
        assert (v["ne"], v["dup"], v["loss"]) == (1_000_000, 1e-8, 1e-8)
        assert v["base"] == "conf_dup_loss_medium_10M_ne1M"


def test_relabel_for_simphy():
    tree, mapping = relabel_for_simphy("((A,B),(A,C),D);")
    assert sorted(Tree(tree, format=1).get_leaf_names()) == ["A_1", "A_2", "B", "C", "D"]
    assert mapping == {"A_1": "A", "A_2": "A", "B": "B", "C": "C", "D": "D"}


def test_relabel_rejects_underscore():
    with pytest.raises(ValueError):
        relabel_for_simphy("((A_x,B),C);")


def test_restore_simphy_labels():
    _, mapping = relabel_for_simphy("((A,B),(A,C),D);")
    gene = "(((A_1_0_0:1,B_0_0:1):1,(A_2_0_0:1,C_0_0:1):1):1,D_0_0:3);"
    out = Tree(restore_simphy_labels(gene, mapping), format=1)
    assert Counter(out.get_leaf_names()) == Counter({"A": 2, "B": 1, "C": 1, "D": 1})


def test_restore_unknown_label_fails():
    with pytest.raises(ValueError):
        restore_simphy_labels("(Z_0_0:1,Y_0_0:1);", {"A": "A"})


def test_prepare_species_tree_is_singly_labelled_ultrametric():
    nex, mapping = prepare_species_tree("((A,B),(A,C),D);", height=1000, seed=0)
    assert nex.startswith("#NEXUS") and "tree 1=" in nex
    t = Tree(nex.split("tree 1=")[1].split("\n")[0], format=1)
    assert sorted(t.get_leaf_names()) == ["A_1", "A_2", "B", "C", "D"]
    assert {round(t.get_distance(leaf), 6) for leaf in t.iter_leaves()} == {1000.0}
    assert mapping["A_2"] == "A"


def _fake_simphy(raw, replicates, n_genes, stree, gtree):
    for r in range(1, replicates + 1):
        d = raw / str(r)
        d.mkdir(parents=True)
        (d / "s_tree.trees").write_text(stree + "\n")
        for g in range(1, n_genes + 1):
            (d / f"g_trees{g:04d}.trees").write_text(gtree + "\n")


S = "((A_1:1,B:1):1,(A_2:1,C:1):1);"
G = "((A_1_0_0:1,B_0_0:1):1,(A_2_0_0:1,C_0_0:1):1);"
MAP = {"A_1": "A", "A_2": "A", "B": "B", "C": "C"}


def test_write_replicates_layout(tmp_path):
    raw = tmp_path / "raw"
    _fake_simphy(raw, 2, 3, S, G)
    write_replicates(raw, MAP, tmp_path / "data", "conf_x", "Net", {"ne": 200000}, replicates=2, n_trees=3)
    rep = tmp_path / "data" / "simulated" / "conf_x" / "Net" / "rep2"
    lines = (rep / "gene_trees.nwk").read_text().splitlines()
    assert len(lines) == 3
    assert Counter(Tree(lines[0], format=1).get_leaf_names()) == Counter({"A": 2, "B": 1, "C": 1})
    assert "ne=200000" in (rep / "simphy_params.txt").read_text()
    assert (rep / "simphy_species_tree.nwk").is_file() and (rep / "simphy_gene_trees.nwk").is_file()


def test_write_replicates_wrong_gene_count_fails(tmp_path):
    raw = tmp_path / "raw"
    _fake_simphy(raw, 1, 2, S, G)
    with pytest.raises(ValueError):
        write_replicates(raw, MAP, tmp_path / "d", "c", "N", {}, replicates=1, n_trees=3)


def test_write_replicates_with_fractionation(tmp_path):
    raw = tmp_path / "raw"
    _fake_simphy(raw, 1, 5, "((A_1:1,X:1):1,(A_2:1,Y:1):1);", "((A_1_0_0:1,X_0_0:1):1,(A_2_0_0:1,Y_0_0:1):1);")
    mapping = {"A_1": "A", "A_2": "A", "X": "X", "Y": "Y"}
    write_replicates(raw, mapping, tmp_path / "d", "c", "N", {}, replicates=1, n_trees=5, retention=0.0, seed=0)
    lines = (tmp_path / "d/simulated/c/N/rep1/gene_trees.nwk").read_text().splitlines()
    assert len(lines) == 5
    for line in lines:
        assert Counter(Tree(line, format=1).get_leaf_names())["A"] == 1


def test_simulate_network_runs_simphy_once_per_replicate(tmp_path, rates_pkl, monkeypatch):
    from polybench.simulate import simphy

    (tmp_path / "model_networks").mkdir()
    (tmp_path / "model_networks" / "N.tre").write_text("((A,B),(A,C));\n")
    calls = []

    def fake_run(cmd):
        calls.append(cmd)
        out = Path(cmd[cmd.index("-o") + 1]) / "1"
        out.mkdir(parents=True)
        (out / "s_tree.trees").write_text(S + chr(10))
        for g in (1, 2):
            (out / f"g_trees{g:04d}.trees").write_text(G + chr(10))

    monkeypatch.setattr(simphy, "run_simphy", fake_run)
    simphy.simulate_network(tmp_path, tmp_path / "work", "conf_ils_low_10M", "N", simphy_bin="simphy",
                            subst_rates_pkl=rates_pkl, seed=5, n_trees=2, replicates=3)
    assert len(calls) == 3
    assert all(c[c.index("-rs") + 1] == "1" for c in calls)
    seeds = [c[c.index("-cs") + 1] for c in calls]
    assert len(set(seeds)) == 3
    assert len({c[c.index("-su") + 1] for c in calls}) == 1
    for k in (1, 2, 3):
        params = (tmp_path / "simulated/conf_ils_low_10M/N" / f"rep{k}" / "simphy_params.txt").read_text()
        assert f"seed={seeds[k - 1]}" in params


def test_base_config_records_base_subst_rate(tmp_path, rates_pkl):
    from polybench.simulate import simphy

    (tmp_path / "model_networks").mkdir()
    (tmp_path / "model_networks" / "N.tre").write_text("((A,B),(A,C));" + chr(10))
    base = "conf_dup_loss_medium_10M_ne1M"
    _fake_simphy(tmp_path / "work" / base / "N", 1, 2, S, G)
    write_replicates(tmp_path / "work" / base / "N", MAP, tmp_path, base, "N", {"subst_rate": "1.5e-09"},
                     replicates=1, n_trees=2)
    simphy.simulate_network(tmp_path, tmp_path / "work", base + "_fix050", "N", simphy_bin="simphy",
                            subst_rates_pkl=rates_pkl, seed=5, n_trees=2, replicates=1)
    params = (tmp_path / "simulated" / (base + "_fix050") / "N" / "rep1" / "simphy_params.txt").read_text()
    assert "subst_rate=1.5e-09" in params and "base_config=" + base in params


def test_simulate_network_batches(tmp_path, rates_pkl, monkeypatch):
    from polybench.simulate import simphy

    (tmp_path / "model_networks").mkdir()
    (tmp_path / "model_networks" / "N.tre").write_text("((A,B),(A,C));\n")
    calls = []

    def fake_run(cmd):
        # one tree per batch, its branch length marks the batch so the order can be checked
        calls.append(cmd)
        out = Path(cmd[cmd.index("-o") + 1]) / "1"
        out.mkdir(parents=True)
        assert cmd[cmd.index("-rl") + 1] == "f:1"
        b = int(out.parent.name.split("_")[1])
        (out / "s_tree.trees").write_text(S + "\n")
        (out / "g_trees1.trees").write_text(G.replace("C_0_0:1", f"C_0_0:{b}") + "\n")

    monkeypatch.setattr(simphy, "run_simphy", fake_run)
    simphy.simulate_network(tmp_path, tmp_path / "work", "conf_ils_low_10M", "N", simphy_bin="simphy",
                            subst_rates_pkl=rates_pkl, seed=5, n_trees=12, replicates=2, batches=12)
    assert len(calls) == 24
    seeds = [int(c[c.index("-cs") + 1]) for c in calls]
    assert len(set(seeds)) == 24
    rep = tmp_path / "simulated/conf_ils_low_10M/N/rep2"
    lengths = [Tree(t, format=1).search_nodes(name="C")[0].dist
               for t in (rep / "gene_trees.nwk").read_text().splitlines()]
    assert lengths == list(range(1, 13))  # batch_2 before batch_10
    params = (rep / "simphy_params.txt").read_text()
    assert "batches=12" in params and f"batch_seeds={','.join(map(str, seeds[12:]))}" in params
    with pytest.raises(ValueError, match="divide"):
        simphy.simulate_network(tmp_path, tmp_path / "work2", "conf_ils_low_10M", "N", simphy_bin="simphy",
                                subst_rates_pkl=rates_pkl, seed=5, n_trees=10, replicates=1, batches=3)


def test_batch_seeds_unique_across_replicates(tmp_path, rates_pkl, monkeypatch):
    from polybench.simulate import simphy

    (tmp_path / "model_networks").mkdir()
    (tmp_path / "model_networks" / "N.tre").write_text("((A,B),(A,C));\n")
    seeds = []

    def fake_run(cmd):
        seeds.append(int(cmd[cmd.index("-cs") + 1]))
        out = Path(cmd[cmd.index("-o") + 1]) / "1"
        out.mkdir(parents=True)
        (out / "s_tree.trees").write_text(S + "\n")
        (out / "g_trees1.trees").write_text(G + "\n")

    monkeypatch.setattr(simphy, "run_simphy", fake_run)
    simphy.simulate_network(tmp_path, tmp_path / "work", "conf_ils_low_10M", "N", simphy_bin="simphy",
                            subst_rates_pkl=rates_pkl, seed=5, n_trees=1000, replicates=3, batches=1000)
    assert len(seeds) == 3000 and len(set(seeds)) == 3000
