import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from ete3 import Tree

from polybench.evaluate.compare import build_inventory, write_summary
from polybench.executors import bash
from polybench.layout import find_datasets, sequence_taxa
from polybench.prep.formats import taxa_table
from polybench.prep.copy_selection import select_copies
from polybench.simulate.simphy import restore_simphy_labels
from polybench.validate import validate_or_raise

ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tests" / "fixtures" / "legacy"
CONFIG = "conf_ils_low_10M"
LABELS = ["alloppnet", "grampa", "grampa_iter", "mpsugar", "padre", "polyphest_p50", "polyphest_p70"]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "dryad" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pack = load("pack_dryad")
copy_fix = load("affected_by_copy_fix")


@pytest.fixture(scope="module")
def packed(tmp_path_factory):
    out = tmp_path_factory.mktemp("dryad") / "pkg"
    env = {**os.environ, "POLYBENCH_PYTHON": Path(sys.executable).as_posix()}
    r = subprocess.run([bash(), (ROOT / "dryad" / "pack_dryad.sh").as_posix(), "--gene2net-root", LEGACY.as_posix(),
                        "--out", out.as_posix(), "--configs", CONFIG, "--n-trees", "10",
                        "--empirical-sources", (LEGACY / "empirical_sources.yaml").as_posix()],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stdout + r.stderr
    return out, r.stdout


def test_grampa_iter_copied_rule(tmp_path):
    (tmp_path / "grandma_result.tre").write_text("(a,b);\n")
    assert pack.grampa_iter_copied(tmp_path)
    (tmp_path / "final_multree.tre").write_text("(a,b);\n")
    assert not pack.grampa_iter_copied(tmp_path)


def test_packed_layout_validates(packed):
    out, stdout = packed
    data = out / "data"
    datasets = find_datasets(data, "simulated/*/*/rep*") + find_datasets(data, "empirical/*")
    assert [d.id for d in datasets] == [f"simulated/{CONFIG}/NetA_2001/rep1", f"simulated/{CONFIG}/NetB_2002/rep1",
                                        "empirical/Emp_2020"]
    validate_or_raise(datasets)
    assert sorted(p.name for p in (out / "distributions").iterdir()) == ["gtr_parameters_all.pkl",
                                                                          "substitution_rates.pkl"]
    assert sorted(p.name for p in (out / "archives").iterdir()) == [f"results_{CONFIG}.tar.gz",
                                                                     f"simulated_{CONFIG}.tar.gz"]
    assert "1 runs had grandma_result.tre without GRANDMA's final_multree.tre" in stdout


def test_gene_trees_restored_as_simulate_writes(packed):
    out, _ = packed
    rep = out / "data" / "simulated" / CONFIG / "NetB_2002" / "rep1"
    raw = rep.joinpath("simphy_gene_trees.nwk").read_text().splitlines()
    mapping = {"Eve_1": "Eve", "Eve_2": "Eve", "Fay": "Fay", "Gus": "Gus", "Hal": "Hal", "Ivy": "Ivy"}
    assert rep.joinpath("gene_trees.nwk").read_text().splitlines() == [restore_simphy_labels(g, mapping) for g in raw]
    # legacy prep order: batch_1, batch_10, batch_2, ... (shell glob); the stray top-level 1/ is ignored
    assert [g.split(":")[1].split(",")[0] for g in raw] == [str(1 + i / 10) for i in [1, 10, 2, 3, 4, 5, 6, 7, 8, 9]]
    params = dict(line.split("=", 1) for line in rep.joinpath("simphy_params.txt").read_text().splitlines())
    assert params["ne"] == "200000" and params["subst_rate"] == "f:1.234500000000000e-09"


def test_alignments_and_taxon_map(packed):
    out, _ = packed
    aln = out / "data" / "simulated" / CONFIG / "NetA_2001" / "rep1" / "alignments"
    assert len(list(aln.glob("gene_*.fasta"))) == 10
    assert aln.joinpath("taxon_map.tsv").read_text().splitlines()[:2] == ["Ana_1_0_0\tAna\tAna_L0_I0:A",
                                                                         "Ana_2_0_0\tAna\tAna_L0_I0:B"]
    assert not (out / "data" / "simulated" / CONFIG / "NetB_2002" / "rep1" / "alignments").exists()


def test_inventory_statuses(packed):
    out, _ = packed
    inv = build_inventory(out / "results", out / "data", LABELS)
    st = {(r.network, r.label): r.status for r in inv.itertuples()}
    assert st[("NetA_2001", "grampa_iter")] == "completed"
    assert st[("NetB_2002", "grampa_iter")] == "failed"
    assert not (out / "results" / "grampa_iter" / "simulated" / CONFIG / "NetB_2002" / "rep1" / "network.nwk").exists()
    assert st[("NetA_2001", "mpsugar")] == "failed" and st[("NetB_2002", "mpsugar")] == "missing"
    assert st[("NetB_2002", "padre")] == "failed" and st[("NetB_2002", "polyphest_p50")] == "failed"
    assert st[("Emp_2020", "polyphest_p50")] == "completed" and st[("Emp_2020", "padre")] == "failed"
    assert st[("NetB_2002", "alloppnet")] == "missing"
    status = json.loads((out / "results" / "polyphest_p70" / "simulated" / CONFIG / "NetA_2001" / "rep1"
                         / "status.json").read_text())
    assert status["method"] == "polyphest" and status["params"]["percentile"] == 70 and status["runtime_s"] is None


def test_evaluate_consumes_package(packed, tmp_path):
    out, _ = packed
    pkg = tmp_path / "pkg"
    shutil.copytree(out, pkg, ignore=shutil.ignore_patterns("archives"))
    summary, errors = write_summary(pkg / "results", pkg / "data")
    assert errors == 0
    cmp = pd.read_csv(summary / "comparisons.csv")
    assert (cmp["status"] == "SUCCESS").all()
    assert set(cmp["label"]) == {"alloppnet", "grampa", "grampa_iter", "padre", "polyphest_p50", "polyphest_p70"}
    pw = pd.read_csv(summary / "empirical_pairwise.csv")
    assert {"published", "grampa_iter", "polyphest_p50"} <= set(pw["label_a"]) | set(pw["label_b"])


def test_homeolog_column_reproduces_legacy_pairing(tmp_path):
    # a duplicated locus: legacy pairs within (locus, individual), the generic rule pairs in sorted id order
    aln = tmp_path / "rep" / "1" / "alignments"
    aln.mkdir(parents=True)
    ids = ["Ana_1_0_0", "Ana_1_1_0", "Ana_2_0_0", "Bel_0_0"]
    aln.joinpath("alignment_1.phy").write_text(" 4 4\n" + "".join(f"{i}  ACGT\n" for i in ids))
    legacy_in = tmp_path / "alloppnet_input"
    legacy_in.mkdir()
    legacy_in.joinpath("ploidy_level.json").write_text(json.dumps({"Ana": 4, "Bel": 2}))
    legacy = ["ID species individual genome", "Ana_1_0_0 Ana Ana_L0_I0 A", "Ana_2_0_0 Ana Ana_L0_I0 B",
              "Ana_1_1_0 Ana Ana_L1_I0 A", "Ana_1_1_0_miss Ana Ana_L1_I0 B", "Bel_0_0 Bel Bel_0_0_ind A"]
    legacy_in.joinpath("taxa_table.txt").write_text("\n".join(legacy))
    rep = pack.SimRep(CONFIG, "Net", 1, tmp_path / "rep", [], None, None, legacy_in)
    out = tmp_path / "packed"
    pack.write_alignments(rep, {"Ana_1": "Ana", "Ana_2": "Ana", "Bel": "Bel"}, out)
    assert "Ana_1_1_0\tAna\tAna_L1_I0:A" in out.joinpath("taxon_map.tsv").read_text().splitlines()
    genes = {Path(f).stem: m for f, m in sequence_taxa(out).items()}
    assert sorted(taxa_table(genes, {"Ana": 2, "Bel": 1}).splitlines()[1:]) != sorted(legacy[1:])  # generic differs
    legacy_in.joinpath("taxa_table.txt").write_text("\n".join(legacy[:3] + ["Ana_1_1_0 Ana Ana_L0_I0 A"] + legacy[5:]))
    with pytest.raises(ValueError, match="two sequences for genome A"):
        pack.write_alignments(rep, {"Ana_1": "Ana", "Ana_2": "Ana", "Bel": "Bel"}, tmp_path / "packed2")


def test_fill_placeholders_raise(tmp_path):
    with pytest.raises(ValueError, match="FILL"):
        pack.plan_empirical(LEGACY, ROOT / "dryad" / "empirical_sources.yaml", pack.Plan())


K3 = "((a:1,a:1):5,(a:1,(a:1,a:1):1):5);"


def test_legacy_selection_differs_for_k3_when_id_order_disagrees():
    tree = Tree(K3, format=1)
    leaves = list(tree.iter_leaves())
    position = {leaf: i for i, leaf in enumerate(leaves)}
    fixed = set(select_copies(tree, "a", 3))
    assert set(copy_fix.legacy_select(tree, leaves, 3, order=position.get)) == fixed
    assert set(copy_fix.legacy_select(tree, leaves, 3, order=lambda n: -position[n])) != fixed
    assert copy_fix.tree_changes(tree, {"a": 3}, order=lambda n: -position[n])
    assert not copy_fix.tree_changes(tree, {"a": 2}, order=lambda n: -position[n])


def test_affected_lists_datasets_with_three_copies(tmp_path, capsys):
    three = "((a:1,a:1):5,(a:1,b:1):5);"
    (tmp_path / "High").mkdir()
    five = "((a:1,a:1):5,(a:1,(a:1,(a:1,b:1):1):1):5);"
    (tmp_path / "High" / "gene_trees.nwk").write_text("\n".join([three] * 3 + [five]))
    (tmp_path / "Low").mkdir()
    (tmp_path / "Low" / "gene_trees.nwk").write_text("((a:1,a:1):1,b:1);\n")
    (tmp_path / "Empty").mkdir()
    copy_fix.main([str(tmp_path)])
    text = capsys.readouterr().out
    assert "High: taxa with >= 3 copies {'a': 3};" in text and "of 4 gene trees change" in text
    assert "Low:" not in text and "no gene trees (gene_trees.nwk) in: Empty" in text


def test_affected_with_legacy_taxa_copies_json(tmp_path, capsys):
    # labels carry accessions; k and leaf -> taxon come from the JSON (longest prefix, as reduce_diverse_k.py)
    d = tmp_path / "Acc_2020"
    d.mkdir()
    tree = "((sin_K_1:1,sin_K_2:1):5,(sin_K_3:1,(sin_K_4:1,(sin_1:1,b_x:1):1):1):5);"
    (d / "Acc_2020_rooted.tre").write_text(tree + "\n")
    (d / "taxa_copies.json").write_text(json.dumps({"sin_K": 3, "sin": 1}))
    copy_fix.main([str(tmp_path), "--trees", "{dataset}_rooted.tre", "--ploidy-json", "taxa_copies.json"])
    assert "Acc_2020: taxa with >= 3 copies {'sin_K': 3}; " in capsys.readouterr().out
    assert copy_fix.legacy_taxon("sin_K_1", {"sin_K": 3, "sin": 1}) == "sin_K"
    assert copy_fix.legacy_taxon("b_x", {"sin": 1}) == "b"
    trees_map = tmp_path / "trees.yaml"
    trees_map.write_text("Acc_2020: Acc_2020_rooted.tre\n")
    (tmp_path / "Other").mkdir()
    copy_fix.main([str(tmp_path), "--trees", str(trees_map)])
    text = capsys.readouterr().out
    assert "Acc_2020:" not in text and "in: Other" in text  # labels as taxa: every leaf is its own taxon
    with pytest.raises(FileNotFoundError, match="taxa_copies"):
        copy_fix.main([str(tmp_path), "--trees", str(trees_map), "--ploidy-json", "missing.json"])


def test_dry_run_compares_copied_runs_with_grampa(tmp_path, capsys):
    legacy = tmp_path / "legacy"
    shutil.copytree(LEGACY, legacy)
    args = ["--gene2net-root", str(legacy), "--out", str(tmp_path / "out"), "--configs", CONFIG, "--n-trees", "10",
            "--empirical-sources", str(LEGACY / "empirical_sources.yaml"), "--dry-run"]
    pack.main(args)
    assert "copied GRAMPAIter runs vs grampa_result.tre: 1 match, 0 mismatch" in capsys.readouterr().out
    res = legacy / "simulations" / "simulations" / "NetB_2002" / "results" / CONFIG
    (res / "grandma_split" / "replicate_1" / "grandma_result.tre").write_text("((Eve,Ivy),(Fay,(Gus,(Eve,Hal))));\n")
    pack.main(args)
    out = capsys.readouterr().out
    assert "0 match, 1 mismatch" in out and f"mismatch: simulated/{CONFIG}/NetB_2002/rep1" in out
    assert not (tmp_path / "out").exists()


def test_empty_final_multree_counts_as_copied(tmp_path):
    (tmp_path / "grandma_result.tre").write_text("(a,b);\n")
    (tmp_path / "final_multree.tre").write_text("\n")
    assert pack.grampa_iter_copied(tmp_path)


def test_config_map_and_zero_replicates(tmp_path, capsys):
    legacy = tmp_path / "legacy"
    shutil.copytree(LEGACY, legacy)
    for net in ("NetA_2001", "NetB_2002"):
        for sub in ("data", "results", "processed"):
            d = legacy / "simulations" / "simulations" / net / sub
            (d / CONFIG).rename(d / "conf_ils_low_10M_old")
    args = ["--gene2net-root", str(legacy), "--out", str(tmp_path / "out"), "--configs", CONFIG, "--n-trees", "10",
            "--empirical-sources", str(LEGACY / "empirical_sources.yaml"), "--dry-run"]
    with pytest.raises(FileNotFoundError, match="no replicates"):
        pack.main(args)
    pack.main(args + ["--config-map", f"conf_ils_low_10M_old={CONFIG}"])
    out = capsys.readouterr().out
    assert f"simulated {CONFIG}: 2 replicates" in out and "1 match" in out
    with pytest.raises(ValueError, match="LEGACY=NEW"):
        pack.main(args + ["--config-map", "conf_ils_low_10M_old=unknown"])


def test_grampa_from_history_as_legacy_script(tmp_path):
    history = tmp_path / "history.json"
    history.write_text(json.dumps({"(0, 0)": {"best_mt": "((QloX*<1>,c+)<2>,(QloX,Qlo03)<3>)<4>"}}))
    taxa_map = tmp_path / "taxa_map.txt"
    taxa_map.write_text("Qlo\tQloX\n")
    assert pack.grampa_from_history(history, taxa_map) == "((Qlo,c),(Qlo,Qlo03));"
    assert pack.grampa_from_history(history, tmp_path / "none.txt") == "((QloX,c),(QloX,Qlo03));"
    history.write_text(json.dumps({"(1, 0)": {"best_mt": "(a,b);"}}))
    with pytest.raises(ValueError, match=r"iteration \(0, 0\)"):
        pack.grampa_from_history(history, taxa_map)
