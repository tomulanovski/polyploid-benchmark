import pickle

import pytest

from polybench.simulate.alisim import alisim_command, iqtree_model, sample_gtr_params

P = {"AC": 1.0, "AG": 2.0, "AT": 0.5, "CG": 0.8, "CT": 2.5, "GT": 1.0,
     "pi_A": 0.3, "pi_C": 0.2, "pi_G": 0.2, "pi_T": 0.3, "alpha": 0.9, "alignment_length": 700}


def test_iqtree_model_gamma():
    assert iqtree_model(P) == ("GTR{1.000000/2.000000/0.500000/0.800000/2.500000/1.000000}"
                               "+FU{0.300000/0.200000/0.200000/0.300000}+G{0.900000}")


def test_iqtree_model_large_alpha_has_no_gamma():
    assert "+G" not in iqtree_model({**P, "alpha": 3.5})


def test_alisim_command():
    cmd = alisim_command("g.tree", "out/gene_0001", iqtree_model(P), 700, 5, iqtree_bin="iqtree3")
    assert cmd[0] == "iqtree3" and cmd[1:3] == ["--alisim", "out/gene_0001"]
    assert cmd[cmd.index("-t") + 1] == "g.tree"
    assert cmd[cmd.index("--length") + 1] == "700" and cmd[cmd.index("-seed") + 1] == "5"
    assert cmd[cmd.index("-af") + 1] == "fasta"


def test_sample_gtr_params_pools_datasets(tmp_path):
    p = tmp_path / "gtr.pkl"
    p.write_bytes(pickle.dumps({"a": [P], "b": [{**P, "alignment_length": 900}]}))
    assert sample_gtr_params(p, 1)["alignment_length"] in (700, 900)
    assert sample_gtr_params(p, 1) == sample_gtr_params(p, 1)


def test_sample_gtr_params_empty_fails(tmp_path):
    p = tmp_path / "gtr.pkl"
    p.write_bytes(pickle.dumps({}))
    with pytest.raises(ValueError):
        sample_gtr_params(p, 1)


@pytest.mark.filterwarnings("ignore:.*fewer than 10:UserWarning")  # small fixtures on purpose
def test_write_taxon_map_and_validate_simulated_replicate(tmp_path):
    from polybench.layout import find_datasets
    from polybench.simulate.alisim import write_taxon_map
    from polybench.simulate.simphy import write_replicates
    from polybench.validate import validate_dataset

    raw = tmp_path / "raw" / "1"
    raw.mkdir(parents=True)
    (raw / "s_tree.trees").write_text("((A_1:1,B:1):1,(A_2:1,C:1):1);\n")
    gene = "((A_1_0_0:1,B_0_0:1):1,(A_2_0_0:1,C_0_0:1):1);"
    for g in (1, 2):
        (raw / f"g_trees{g:04d}.trees").write_text(gene + "\n")
    mapping = {"A_1": "A", "A_2": "A", "B": "B", "C": "C"}
    write_replicates(tmp_path / "raw", mapping, tmp_path, "conf_x", "N", {}, replicates=1, n_trees=2)
    rep = tmp_path / "simulated" / "conf_x" / "N" / "rep1"
    aln = rep / "alignments"
    aln.mkdir()
    for g in (1, 2):
        (aln / f"gene_{g:04d}.fasta").write_text(">A_1_0_0\nAC\n>A_2_0_0\nAG\n>B_0_0\nAT\n>C_0_0\nAA\n".replace("\n", "\n"))
    simphy_genes = (rep / "simphy_gene_trees.nwk").read_text().splitlines()
    write_taxon_map(simphy_genes, mapping, aln / "taxon_map.tsv")
    assert (aln / "taxon_map.tsv").read_text().splitlines() == ["A_1_0_0\tA", "A_2_0_0\tA", "B_0_0\tB", "C_0_0\tC"]
    (tmp_path / "model_networks").mkdir()
    (tmp_path / "model_networks" / "N.tre").write_text("((A,B),(A,C));" + chr(10))
    (ds,) = find_datasets(tmp_path, "simulated/conf_x/N/rep1")
    assert validate_dataset(ds, ["alignments"]) == []


def test_simulate_alignments_one_param_set_per_dataset(tmp_path, monkeypatch):
    from polybench.simulate import alisim

    pool = [{**P, "alignment_length": 100 * n, "alpha": 0.5 + n / 10} for n in range(1, 20)]
    pkl = tmp_path / "gtr.pkl"
    pkl.write_bytes(pickle.dumps({"a": pool}))
    cmds = []

    def fake_run(cmd, check):
        cmds.append(cmd)
        (tmp_path / "aln" / (cmd[2].rsplit("/", 1)[-1].rsplit(chr(92), 1)[-1] + ".fa")).write_text(">x" + chr(10) + "A" + chr(10))

    monkeypatch.setattr(alisim.subprocess, "run", fake_run)
    genes = ["((A_1_0_0:1,B_0_0:1):1,C_0_0:1);"] * 4
    mapping = {"A_1": "A", "B": "B", "C": "C"}
    alisim.simulate_alignments(genes, mapping, tmp_path / "aln", pkl, seed=7)
    assert len({(c[c.index("-m") + 1], c[c.index("--length") + 1]) for c in cmds}) == 1
    assert [int(c[c.index("-seed") + 1]) for c in cmds] == [700_001, 700_002, 700_003, 700_004]
    text = (tmp_path / "aln" / "alisim_params.txt").read_text()
    assert "model=GTR" in text and "alignment_length=" in text
    assert (tmp_path / "aln" / "gene_0004.fasta").is_file()


@pytest.mark.filterwarnings("ignore:.*fewer than 10:UserWarning")  # small fixtures on purpose
def test_simulate_cli_writes_alignments_per_replicate(tmp_path, monkeypatch):
    from pathlib import Path

    from polybench import cli
    from polybench.layout import find_datasets
    from polybench.simulate import alisim, simphy
    from polybench.validate import validate_dataset

    (tmp_path / "model_networks").mkdir()
    (tmp_path / "model_networks" / "N.tre").write_text("((A,B),(A,C));\n")
    (tmp_path / "gtr.pkl").write_bytes(pickle.dumps({"a": [P]}))
    (tmp_path / "rates.pkl").write_bytes(pickle.dumps({"a": [{"rate": 1e-9}]}))
    (tmp_path / "c.yaml").write_text("data_root: .\nresults_root: res\ntools: {simphy: simphy, iqtree: iqtree3}\n"
                                     "executor: {type: local}\n")

    def fake_simphy(cmd):
        out = Path(cmd[cmd.index("-o") + 1]) / "1"
        out.mkdir(parents=True)
        (out / "s_tree.trees").write_text("((A_1:1,B:1):1,(A_2:1,C:1):1);\n")
        for g in (1, 2):
            (out / f"g_trees{g:04d}.trees").write_text("((A_1_0_0:1,B_0_0:1):1,(A_2_0_0:1,C_0_0:1):1);\n")

    seeds = []

    def fake_alisim(cmd, check):
        assert cmd[0] == "iqtree3"
        seeds.append(int(cmd[cmd.index("-seed") + 1]))
        Path(cmd[2] + ".fa").write_text(">A_1_0_0\nAC\n>A_2_0_0\nAG\n>B_0_0\nAT\n>C_0_0\nAA\n")

    monkeypatch.setattr(simphy, "run_simphy", fake_simphy)
    monkeypatch.setattr(alisim.subprocess, "run", fake_alisim)
    cli.main(["simulate", "--config", str(tmp_path / "c.yaml"), "--network", "N", "--sim-config", "conf_ils_low_10M",
              "--work-dir", str(tmp_path / "work"), "--subst-rates", str(tmp_path / "rates.pkl"), "--seed", "3",
              "--n-trees", "2", "--replicates", "2", "--alignments", "--gtr-params", str(tmp_path / "gtr.pkl")])
    assert len(seeds) == 4 and len(set(seeds)) == 4
    for ds in find_datasets(tmp_path, "simulated/conf_ils_low_10M/N/rep*"):
        assert (ds.path / "alignments" / "taxon_map.tsv").is_file()
        assert sorted(f.name for f in (ds.path / "alignments").glob("*.fasta")) == ["gene_0001.fasta", "gene_0002.fasta"]
        assert validate_dataset(ds, ["alignments"]) == []


def test_simulate_cli_alignments_need_gtr_params(tmp_path):
    from polybench import cli

    (tmp_path / "c.yaml").write_text("data_root: .\nresults_root: res\ntools: {simphy: s, iqtree: i}\n"
                                     "executor: {type: local}\n")
    with pytest.raises(SystemExit) as e:
        cli.main(["simulate", "--config", str(tmp_path / "c.yaml"), "--network", "N", "--sim-config",
                  "conf_ils_low_10M", "--work-dir", "w", "--subst-rates", "r.pkl", "--seed", "1", "--alignments"])
    assert e.value.code == 1
