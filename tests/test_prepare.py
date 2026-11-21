import json
import os
import shutil
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path

import pytest
from ete3 import Tree

from conftest import canonical, cluster_dir
from polybench import runner
from polybench.layout import Dataset
from polybench.runner import load_method, run_one

FIX = Path(__file__).parent / "fixtures"
GOOD = FIX / "bad_datasets" / "good"
FAKE = FIX / "fake_tools"

# every file each run.sh reads from <out_dir>/work/inputs
RUN_INPUTS = {
    "padre": ["padre_trees.tre"],
    "polyphest": ["polyphest_trees.tre", "multi_set.txt"],
    "grampa": ["grampa_trees.tre", "species_tree.tre", "taxa_map.tsv"],
    "grampa_iter": ["grampa_trees.tre", "species_tree.tre", "taxa_map.tsv"],
    "mpsugar": ["mpsugar_trees.nex", "taxon_map.json"],
    "mpallopp_phylonet": ["mpallopp.nex"],
}


def make_dataset(tmp_path, src=GOOD, name="good"):
    root = tmp_path / "data"
    d = root / "empirical" / name
    shutil.copytree(src, d)
    return Dataset("empirical", None, name, None, d, root)


def cfg_for(tmp_path, tools):
    return {"data_root": str(tmp_path / "data"), "results_root": str(tmp_path / "res"), "tools": tools,
            "executor": {"type": "local", "slurm_header": [], "conda_activate": ""}}


def prepared(name, dataset, tmp_path):
    spec = load_method(name)
    work = tmp_path / "work"
    work_ds = runner.prepare_inputs(spec, dataset, cfg_for(tmp_path, {}), work)
    runner.method_prepare(spec, work_ds, work, spec.params)
    return work / "inputs"


@pytest.mark.parametrize("name", RUN_INPUTS)
def test_prepare_writes_every_run_input(name, tmp_path):
    inputs = prepared(name, make_dataset(tmp_path), tmp_path)
    run_sh = (runner.METHODS_DIR / name / "run.sh").read_text()
    for f in RUN_INPUTS[name]:
        assert f in run_sh
        assert (inputs / f).stat().st_size > 0, f


def test_alloppnet_prepare_writes_every_run_input(tmp_path):
    # run.sh hands the inputs folder to generate_beast_xml.r, which reads every .nex and taxa_table.txt
    inputs = prepared("alloppnet", alloppnet_dataset(tmp_path), tmp_path)
    assert "taxa_table.txt" in (runner.METHODS_DIR / "alloppnet" / "generate_beast_xml.r").read_text()
    for f in ["taxa_table.txt", "g1.nex", "g2.nex"]:
        assert (inputs / f).stat().st_size > 0, f


def test_polyphest_inputs(tmp_path):
    inputs = prepared("polyphest", make_dataset(tmp_path), tmp_path)
    assert Counter((inputs / "multi_set.txt").read_text().split()) == {"a": 2, "b": 2, "c": 4, "d": 2}
    trees = (inputs / "polyphest_trees.tre").read_text().splitlines()
    assert trees == [l for l in (GOOD / "gene_trees.nwk").read_text().splitlines() if l.strip()]


def test_grampa_inputs_number_copies(tmp_path):
    inputs = prepared("grampa", make_dataset(tmp_path), tmp_path)
    first = (inputs / "grampa_trees.tre").read_text().splitlines()[1]
    assert sorted(l.name for l in Tree(first, format=1)) == ["1_a", "1_b", "1_c", "1_d", "2_c"]
    species = Tree((inputs / "species_tree.tre").read_text(), format=1)
    assert sorted(species.get_leaf_names()) == ["a", "b", "c", "d"]


def test_polyphest_needs_filter_strategy(tmp_path):
    from polybench.validate import validate_or_raise
    spec = load_method("polyphest")
    with pytest.raises(ValueError, match="filter_strategy"):
        validate_or_raise([make_dataset(tmp_path)], replace(spec, requires=["gene_trees"]))


def test_grampa_rejects_unknown_mode(tmp_path):
    from polybench.validate import validate_or_raise
    spec = load_method("grampa")
    with pytest.raises(ValueError, match="mode"):
        validate_or_raise([make_dataset(tmp_path)], replace(spec, requires=["gene_trees"], params={"mode": "x"}))


def _python_on_path(monkeypatch):
    # run.sh calls `python` for tools that run in their own environment
    monkeypatch.setenv("PATH", str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"])


def test_polyphest_run_passes_params(tmp_path, monkeypatch):
    _python_on_path(monkeypatch)
    spec = replace(load_method("polyphest"), env=None)
    out = tmp_path / "out"
    params = {"filter_strategy": "percentile", "percentile": 50, "use_near_isomorphic": True,
              "isomorphic_threshold": 0.2, "use_clusters": False}
    st = run_one(spec, make_dataset(tmp_path), out, params, "polyphest_p50", 120,
                 cfg_for(tmp_path, {"polyphest": str(FAKE / "Polyphest.py")}))
    assert st["status"] == "completed", (out / "log.txt").read_text()
    argv = json.loads((out / "work" / "polyphest" / "argv.json").read_text())
    tail = argv[argv.index("--filter_strategy"):]
    assert tail == ["--filter_strategy", "percentile", "--percentile", "50", "--use_near_isomorphic", "True",
                    "--isomorphic_threshold", "0.2"]
    nwk = (out / "network.nwk").read_text()
    assert Counter(l.name for l in Tree(nwk.strip(), format=1)) == {"a": 2, "b": 2, "c": 4, "d": 2}


SUBSTRING_TREES = ["((Qlo:1,Qlo03:1):1,(c:1,d:1):1);", "(((Qlo:1,c:1):1,Qlo03:1):1,(Qlo:1,d:1):1);"] * 5


@pytest.mark.parametrize("name,params,flags", [
    ("grampa", {}, ["--mode", "single", "--cap", "15"]),
    ("grampa", {"mode": "original"}, []),
    ("grampa_iter", {}, ["--mode", "split", "--cap", "15", "--iter", "50"]),
])
def test_grampa_substring_labels(name, params, flags, tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    (src / "gene_trees.nwk").write_text("\n".join(SUBSTRING_TREES) + "\n")
    (src / "species_tree.nwk").write_text("((Qlo,Qlo03),(c,d));\n")
    out = tmp_path / "out"
    st = run_one(load_method(name), make_dataset(tmp_path, src, "subs"), out, params, name, 120,
                 cfg_for(tmp_path, {"grandma": str(FAKE / "gran.py"), "grampa": str(FAKE / "grampa")}))
    assert st["status"] == "completed", (out / "log.txt").read_text()
    if params.get("mode") == "original":
        assert (out / "work" / "grampa" / "grampa-scores.txt").is_file()
        assert not list((out / "work").rglob("argv.json"))
    grampa_trees = (out / "work" / "inputs" / "grampa_trees.tre").read_text()
    assert "1_QloX" in grampa_trees and "2_QloX" in grampa_trees and "1_Qlo03" in grampa_trees
    if flags:
        argv = json.loads(next((out / "work").rglob("argv.json")).read_text())
        assert all(f in argv for f in flags)
    leaves = [l.name for l in Tree((out / "network.nwk").read_text().strip(), format=1)]
    assert sorted(leaves) == ["Qlo", "Qlo03", "c", "d"]


# legacy prepared input file -> (our file, how to compare after parsing)
def _trees(path):
    return [canonical(l) for l in Path(path).read_text().splitlines() if l.strip()]


def _multiset(path):
    return Counter(Path(path).read_text().split())


def _nexus_trees(path):
    return [canonical(l.split("=", 1)[1]) for l in Path(path).read_text().splitlines() if l.startswith("Tree ")]


def _json(path):
    return json.loads(Path(path).read_text())


def _table_species(path):
    # ID -> species. The pairing into individuals is not compared, legacy paired by SimPhy locus ids.
    rows = [l.split() for l in Path(path).read_text().splitlines()[1:] if l.strip()]
    return {r[0]: r[1] for r in rows if not r[0].endswith("_miss")}


def _nexus_matrices(folder):
    # {file name: {sequence id: sequence}} for every .nex in the folder
    out = {}
    for f in sorted(Path(folder).glob("*.nex")):
        rows = f.read_text().split("MATRIX", 1)[1].split(";", 1)[0].split()
        out[f.name] = dict(zip(rows[::2], rows[1::2]))
    return out


LEGACY_INPUTS = {
    "padre": [("padre_trees.tre", "padre_trees.tre", _trees)],
    "polyphest": [("polyphest_trees.tre", "polyphest_trees.tre", _trees),
                  ("multi_set.txt", "multi_set.txt", _multiset)],
    "grampa": [("grampa_trees.tre", "grampa_trees.tre", _trees), ("species.tre", "species_tree.tre", _trees)],
    "grampa_iter": [("grampa_trees.tre", "grampa_trees.tre", _trees), ("species.tre", "species_tree.tre", _trees)],
    "mpsugar": [("mpsugar_trees.nex", "mpsugar_trees.nex", _nexus_trees), ("taxon_map.json", "taxon_map.json", _json)],
    # legacy alignment_<n>.phy became dataset alignments/alignment_<n>.fasta, so the .nex names match
    "alloppnet": [("taxa_table.txt", "taxa_table.txt", _table_species), (".", ".", _nexus_matrices)],
}


@pytest.mark.needs_cluster_fixture
@pytest.mark.parametrize("name", LEGACY_INPUTS)
def test_prepared_inputs_match_legacy(name, tmp_path):
    d = cluster_dir(name)
    inputs = prepared(name, make_dataset(tmp_path, d / "dataset", "legacy"), tmp_path)
    for legacy, ours, parse in LEGACY_INPUTS[name]:
        assert parse(inputs / ours) == parse(d / "legacy_input" / legacy), legacy


def _fake_bin_on_path(monkeypatch):
    # java and Rscript come from the method's environment, the fakes stand in for them
    monkeypatch.setenv("PATH", str(FAKE / "bin") + os.pathsep + os.environ["PATH"])


def test_mpsugar_inputs(tmp_path):
    inputs = prepared("mpsugar", make_dataset(tmp_path), tmp_path)
    assert json.loads((inputs / "taxon_map.json").read_text()) == {
        "a": ["a_1"], "b": ["b_1"], "c": ["c_1", "c_2"], "d": ["d_1"]}
    nexus = (inputs / "mpsugar_trees.nex").read_text()
    assert nexus.count("Tree geneTree") == 10 and "Tree geneTree2 = (((a_1,b_1),c_1),(c_2,d_1));" in nexus


def test_mpsugar_run(tmp_path, monkeypatch):
    _python_on_path(monkeypatch)
    monkeypatch.setenv("PYTHONPATH", str(FAKE))
    spec = replace(load_method("mpsugar"), env=None)
    out = tmp_path / "out"
    st = run_one(spec, make_dataset(tmp_path), out, {"chains": 2, "seed": 5}, "mpsugar", 120, cfg_for(tmp_path, {}))
    assert st["status"] == "completed", (out / "log.txt").read_text()
    calls = json.loads((out / "work" / "inputs" / "mp_sugar_calls.json").read_text())
    assert [(c["iter_ct"], c["seed"]) for c in calls] == [(500, 5), (500, 6)]
    assert (out / "network.nwk").read_text().strip() == "(a,b,c,d);"


def test_mpallopp_phylonet_inputs(tmp_path):
    inputs = prepared("mpallopp_phylonet", make_dataset(tmp_path), tmp_path)
    nexus = (inputs / "mpallopp.nex").read_text()
    assert "InferNetwork_MP_Allopp (all) 1 -x 1 -pl {threads} -di -a <a:a; b:b; c:c,c_1,c_2; d:d>;" in nexus
    assert "Tree geneTree2 = (((a,b),c_1),(c_2,d));" in nexus


def test_mpallopp_phylonet_run(tmp_path, monkeypatch):
    _fake_bin_on_path(monkeypatch)
    out = tmp_path / "out"
    spec = load_method("mpallopp_phylonet")
    spec = replace(spec, resources={**spec.resources, "threads": 4})
    st = run_one(spec, make_dataset(tmp_path), out, {"command": "InferNetwork_MP_Allopp (all) 2 -x 3 -pl {threads} -di"},
                 "mpallopp_phylonet", 120, cfg_for(tmp_path, {"phylonet_jar": "PhyloNet.jar"}))
    assert st["status"] == "completed", (out / "log.txt").read_text()
    assert "InferNetwork_MP_Allopp (all) 2 -x 3 -pl 4 -di -a <" in (out / "work" / "mpallopp_run.nex").read_text()
    assert (out / "network.nwk").read_text().strip() == "(a,b,c,d);"


ALLOPP_TREES = ["((a,b),(c,(c,d)));"] * 10
ALLOPP_ALIGNMENTS = {
    "g1.fasta": ">s1\nACGT\n>s2\nACGA\n>s3\nACGG\n>s4 second homeolog\nACGC\n>s5\nACTT\n",
    "g2.fasta": ">s1\nAC\n>s2\nAA\n>s3\nAG\n>s5\nAT\n",
}


def alloppnet_dataset(tmp_path, ploidy="a\t1\nb\t1\nc\t2\nd\t1\n", taxon_map="s1\ta\ns2\tb\ns3\tc\ns4\tc\ns5\td\n"):
    src = tmp_path / "src"
    (src / "alignments").mkdir(parents=True)
    (src / "gene_trees.nwk").write_text("\n".join(ALLOPP_TREES) + "\n")
    (src / "ploidy.tsv").write_text(ploidy)
    for name, text in ALLOPP_ALIGNMENTS.items():
        (src / "alignments" / name).write_text(text)
    (src / "alignments" / "taxon_map.tsv").write_text(taxon_map)
    return make_dataset(tmp_path, src, "allopp")


def test_alloppnet_homeolog_column_overrides_pairing(tmp_path):
    # as the legacy SimPhy-id rule for a duplicated locus: s3 and s4 are the B genomes of two individuals
    taxon_map = "s1\ta\ns2\tb\ns3\tc\tc_L0:B\ns4\tc\tc_L1:B\ns5\td\n"
    inputs = prepared("alloppnet", alloppnet_dataset(tmp_path, taxon_map=taxon_map), tmp_path)
    assert (inputs / "taxa_table.txt").read_text().splitlines() == [
        "ID species individual genome",
        "s1 a s1_ind A",
        "s2 b s2_ind A",
        "s3 c c_L0 B",
        "s3_miss c c_L0 A",
        "s4 c c_L1 B",
        "s4_miss c c_L1 A",
        "s5 d s5_ind A",
    ]


@pytest.mark.parametrize("taxon_map, message", [
    ("s1\ta\ta1:A\ns2\tb\ns3\tc\tc1:A\ns4\tc\tc1:B\ns5\td\n", "not a tetraploid"),
    ("s1\ta\ns2\tb\ns3\tc\tc1:A\ns4\tc\ns5\td\n", "some sequences only"),
    ("s1\ta\ns2\tb\ns3\tc\tc1:A\ns4\tc\tc1:A\ns5\td\n", "two sequences for genome A"),
])
def test_alloppnet_homeolog_column_checked(tmp_path, taxon_map, message):
    from polybench.validate import validate_or_raise
    with pytest.raises(ValueError, match=message):
        validate_or_raise([alloppnet_dataset(tmp_path, taxon_map=taxon_map)], load_method("alloppnet"))


def test_alloppnet_maps_sequences_with_taxon_map(tmp_path):
    inputs = prepared("alloppnet", alloppnet_dataset(tmp_path), tmp_path)
    assert (inputs / "taxa_table.txt").read_text().splitlines() == [
        "ID species individual genome",
        "s1 a s1_ind A",
        "s2 b s2_ind A",
        "s3 c c_I1 A",
        "s4 c c_I1 B",
        "s5 d s5_ind A",
    ]
    assert sorted(p.name for p in inputs.glob("*.nex")) == ["g1.nex", "g2.nex"]
    assert "NTAX=4 NCHAR=2" in (inputs / "g2.nex").read_text()
    assert "s4  ACGC" in (inputs / "g1.nex").read_text()


def test_alloppnet_unpaired_homeolog_gets_missing_partner():
    from polybench.prep.formats import taxa_table
    table = taxa_table({"g1": {"x": "T", "y": "T", "z": "T"}, "g2": {"x": "T", "w": "D"}}, {"T": 2, "D": 1})
    assert table.splitlines()[1:] == ["w D w_ind A", "x T T_I1 A", "y T T_I1 B", "z T T_I2 A", "z_miss T T_I2 B"]


def test_alloppnet_needs_a_tetraploid(tmp_path):
    with pytest.raises(ValueError, match="tetraploid"):
        prepared("alloppnet", alloppnet_dataset(tmp_path, "a\t1\nb\t1\nc\t1\nd\t1\n"), tmp_path)


def test_alloppnet_rejects_hexaploid(tmp_path):
    from polybench.validate import validate_or_raise
    with pytest.raises(ValueError, match="ploidy <= 2"):
        validate_or_raise([alloppnet_dataset(tmp_path, "a\t1\nb\t1\nc\t3\nd\t1\n")], load_method("alloppnet"))


def alloppnet_cfg(tmp_path):
    return cfg_for(tmp_path, {"beast": str(FAKE / "beast"), "treeannotator": str(FAKE / "treeannotator"),
                              "alloppdt": str(tmp_path / "alloppdt")})


def test_alloppnet_run(tmp_path, monkeypatch):
    _fake_bin_on_path(monkeypatch)
    spec = replace(load_method("alloppnet"), env=None)
    out = tmp_path / "out"
    # 1e8 arrives as the float 100000000.0, BEAST needs the integer
    st = run_one(spec, alloppnet_dataset(tmp_path), out, {"chain_length": 1e8, "log_every": 1000.0}, "alloppnet",
                 120, alloppnet_cfg(tmp_path))
    assert st["status"] == "completed", (out / "log.txt").read_text()
    run_dir = out / "work" / "alloppnet"
    assert (run_dir / "alloppnet.XML").read_text().splitlines()[0].endswith(" 100000000 1000 10000")
    assert (run_dir / "burnin.txt").read_text().strip() == "1"  # 11 trees, STATE_0 left out, 10%
    leaves = [l.name for l in Tree((out / "network.nwk").read_text().strip(), format=1)]
    assert sorted(leaves) == ["a", "b", "c", "c", "d"]


def test_alloppnet_summarizes_when_beast_runs_out_of_time(tmp_path, monkeypatch):
    _fake_bin_on_path(monkeypatch)
    monkeypatch.setenv("FAKE_BEAST_SLEEP", "60")
    spec = replace(load_method("alloppnet"), env=None)
    out = tmp_path / "out"
    st = run_one(spec, alloppnet_dataset(tmp_path), out, {"beast_hours": 0.001}, "alloppnet", 120,
                 alloppnet_cfg(tmp_path))
    assert st["status"] == "completed", (out / "log.txt").read_text()
    assert st["runtime_s"] < 40
    # the half-written last tree is dropped before TreeAnnotator, which here keeps the last tree
    assert "BEAST exit code 124 (stopped after beast_hours), 11 sampled trees" in (out / "log.txt").read_text()
    assert (out / "work" / "alloppnet" / "burnin.txt").read_text().strip() == "1"
    leaves = [l.name for l in Tree((out / "network.nwk").read_text().strip(), format=1)]
    assert sorted(leaves) == ["a", "b", "c", "c", "d"]
