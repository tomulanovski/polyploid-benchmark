import json
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import canonical, cluster_dir
from polybench.measures.reticulate_tree import ReticulateTree
from polybench.prep import outputs
from polybench.prep.formats import write_grampa_inputs

OUT = Path(__file__).parent / "fixtures" / "outputs"
GRAMPA_TAXA = {"Qlo", "Qlo03", "Abc", "Def"}


def taxa(newick: str) -> set[str]:
    return set(ReticulateTree(newick, is_multree=True, refold=True).get_leaf_counts())


def test_parse_padre():
    nwk = outputs.parse_padre(OUT / "padre" / "padre_trees-result.tre")
    assert taxa(nwk) == {"A", "B", "C", "D"}


def test_parse_polyphest_takes_multree_line():
    nwk = outputs.parse_polyphest(OUT / "polyphest" / "polyphest_trees-polyphest.txt")
    assert nwk == "((A,(B,C)),(C,D));"
    assert taxa(nwk) == {"A", "B", "C", "D"}


def test_parse_grampa_scores_takes_first_row():
    nwk = outputs.parse_grampa_scores(OUT / "grampa" / "grampa-scores.txt", OUT / "grampa" / "taxa_map.tsv")
    assert nwk == "(((Qlo,Abc),Qlo03),(Qlo,Def));"
    assert taxa(nwk) == GRAMPA_TAXA


def test_parse_grandma_first_iteration():
    nwk = outputs.parse_grandma_first_iteration(OUT / "grampa" / "history.json", OUT / "grampa" / "taxa_map.tsv")
    assert nwk == "(((Qlo,Abc),Qlo03),(Qlo,Def));"
    assert taxa(nwk) == GRAMPA_TAXA


def test_parse_grandma_split():
    nwk = outputs.parse_grandma_split(OUT / "grampa_iter", OUT / "grampa_iter" / "taxa_map.tsv")
    assert nwk == "(((Qlo,Abc),Qlo03),(Qlo,Def));"
    assert taxa(nwk) == GRAMPA_TAXA


def network_taxa(newick: str) -> set[str]:
    return set(ReticulateTree(newick, refold=True).get_leaf_counts())


def test_parse_mpsugar_takes_first_network():
    nwk = outputs.parse_mpsugar(OUT / "mpsugar" / "mpsugar_results.txt")
    assert nwk == "((A,(B)#UID_3)UID_1,(#UID_3,(C,D)UID_2)UID_4)UID_0;"
    assert network_taxa(nwk) == {"A", "B", "C", "D"}


def test_parse_phylonet_takes_dendroscope_line():
    nwk = outputs.parse_phylonet(OUT / "mpallopp_phylonet" / "phylonet_out.txt")
    assert nwk == "((A,(B)#H1),(#H1,(C,D)));"
    assert network_taxa(nwk) == {"A", "B", "C", "D"}


def test_parse_alloppnet_mcc_translates_and_strips_homeolog_suffix():
    d = OUT / "alloppnet"
    nwk = outputs.parse_alloppnet_mcc(d / "mcc.tre", d / "ploidy.tsv")
    assert nwk == "((A,(B,C)),(C,D));"
    assert taxa(nwk) == {"A", "B", "C", "D"}


def test_parse_alloppnet_mcc_without_translate(tmp_path):
    f = tmp_path / "mcc.tre"
    f.write_text("#NEXUS\nBegin trees;\ntree TREE1 = [&R] ((A[&h=0]:1,C0:1):1,(C1:1,D:1)[&posterior=1]:1);\nEnd;\n")
    assert outputs.parse_alloppnet_mcc(f, OUT / "alloppnet" / "ploidy.tsv") == "((A,C),(C,D));"


def test_alloppnet_label_outside_ploidy_fails(tmp_path):
    f = tmp_path / "mcc.tre"
    f.write_text("tree TREE1 = ((A,C0),(C1,Zzz));\n")
    with pytest.raises(ValueError, match="Zzz"):
        outputs.parse_alloppnet_mcc(f, OUT / "alloppnet" / "ploidy.tsv")


def test_mpsugar_without_newick_fails(tmp_path):
    f = tmp_path / "mpsugar_results.txt"
    f.write_text("MP-SUGAR Results\nChain 1 (seed: 1)\nNetworks found: 0\n")
    with pytest.raises(ValueError, match="Newick"):
        outputs.parse_mpsugar(f)


def test_phylonet_without_dendroscope_line_fails(tmp_path):
    f = tmp_path / "phylonet_out.txt"
    f.write_text("Inferred Network #1:\n((A,B),C);\n")
    with pytest.raises(ValueError, match="Dendroscope"):
        outputs.parse_phylonet(f)


def test_polyphest_without_multree_fails(tmp_path):
    f = tmp_path / "x-polyphest.txt"
    f.write_text("network: (A,B);\n")
    with pytest.raises(ValueError, match="multree"):
        outputs.parse_polyphest(f)


def test_history_without_first_iteration_fails(tmp_path):
    f = tmp_path / "history.json"
    f.write_text(json.dumps({"(1, 0)": {"best_mt": "(A,B);"}}))
    with pytest.raises(ValueError, match=r"\(0, 0\)"):
        outputs.parse_grandma_first_iteration(f, OUT / "grampa" / "taxa_map.tsv")


def test_grampa_label_outside_map_fails(tmp_path):
    (tmp_path / "final_multree.tre").write_text("((QloX,Abc),Zzz);")
    with pytest.raises(ValueError, match="Zzz"):
        outputs.parse_grandma_split(tmp_path, OUT / "grampa" / "taxa_map.tsv")


def test_two_trees_fail(tmp_path):
    f = tmp_path / "padre_trees-result.tre"
    f.write_text("(A,B);\n(A,C);\n")
    with pytest.raises(ValueError, match="one Newick"):
        outputs.parse_padre(f)


def test_empty_output_fails(tmp_path):
    f = tmp_path / "padre_trees-result.tre"
    f.write_text("\n")
    with pytest.raises(ValueError, match="no tree"):
        outputs.parse_padre(f)


def test_cli_prints_newick():
    res = subprocess.run([sys.executable, "-m", "polybench.prep.outputs", "grandma_split", str(OUT / "grampa_iter"),
                          str(OUT / "grampa_iter" / "taxa_map.tsv")], capture_output=True, text=True, check=True)
    assert res.stdout.strip() == "(((Qlo,Abc),Qlo03),(Qlo,Def));"


def test_cli_unknown_parser_fails():
    res = subprocess.run([sys.executable, "-m", "polybench.prep.outputs", "nope", "x"], capture_output=True, text=True)
    assert res.returncode != 0 and "usage" in res.stderr


def _grampa_map(d: Path, tmp_path: Path) -> Path:
    write_grampa_inputs(d / "dataset", tmp_path)
    return tmp_path / "taxa_map.tsv"


CLUSTER_PARSE = {
    "padre": lambda d, tmp: outputs.parse_padre(d / "legacy_output" / "padre_trees-result.tre"),
    "polyphest": lambda d, tmp: outputs.parse_polyphest(d / "legacy_output" / "polyphest_trees-polyphest.txt"),
    "grampa": lambda d, tmp: outputs.parse_grandma_first_iteration(d / "legacy_output" / "history.json",
                                                                   _grampa_map(d, tmp)),
    "grampa_iter": lambda d, tmp: outputs.parse_grandma_split(d / "legacy_output", _grampa_map(d, tmp)),
    "mpsugar": lambda d, tmp: outputs.parse_mpsugar(d / "legacy_output" / "mpsugar_results.txt"),
    "mpallopp_phylonet": lambda d, tmp: outputs.parse_phylonet(d / "legacy_output" / "phylonet_out.txt"),
    "alloppnet": lambda d, tmp: outputs.parse_alloppnet_mcc(d / "legacy_output" / "alloppnet_consensus.tre",
                                                            d / "dataset" / "ploidy.tsv"),
}
NETWORK_OUTPUT = ("mpsugar", "mpallopp_phylonet")


@pytest.mark.needs_cluster_fixture
@pytest.mark.parametrize("method", CLUSTER_PARSE)
def test_parser_matches_legacy_postprocess(method, tmp_path):
    d = cluster_dir(method)
    nwk = CLUSTER_PARSE[method](d, tmp_path)
    assert canonical(nwk) == canonical((d / "legacy_result.tre").read_text())
    network_taxa(nwk) if method in NETWORK_OUTPUT else taxa(nwk)
