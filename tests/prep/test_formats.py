import json

import pytest

from polybench.prep.formats import (
    alignment_to_nexus,
    clean_grampa_tree,
    read_fasta,
    strip_alloppnet_suffixes,
    strip_copy_suffixes,
    to_grampa,
    to_mpsugar_nexus,
    to_phylonet_nexus,
)
from polybench.prep.species_tree import astral_command


def test_to_grampa_numbers_copies():
    out = to_grampa(["((Rch:1,Rch:1):1,Abc:1);"], {})
    assert out == ["((1_Rch:1,2_Rch:1):1,1_Abc:1);"]


def test_to_grampa_applies_mapping():
    out = to_grampa(["(Qlo:1,Qlo03:1);"], {"Qlo": "QloX"})
    assert out == ["(1_QloX:1,1_Qlo03:1);"]


def test_clean_grampa_tree():
    assert clean_grampa_tree("((Lamium*,Lamium+)<5>,Abc)<9>;") == "((Lamium,Lamium),Abc);"


def test_mpsugar_taxon_map():
    nexus, taxon_map = to_mpsugar_nexus(["((A:1,A:1):1,B:1);", "(A:1,B:1);"])
    assert json.loads(taxon_map) == {"A": ["A_1", "A_2"], "B": ["B_1"]}
    assert nexus.startswith("#NEXUS\nBEGIN TREES;\n")
    assert "Tree geneTree1 = ((A_1,A_2),B_1);" in nexus
    assert "Tree geneTree2 = (A_1,B_1);" in nexus


def test_mpsugar_taxon_map_orders_copies_numerically():
    tree = "(" + ",".join(["A"] * 11) + ",B);"
    _, taxon_map = to_mpsugar_nexus([tree])
    assert json.loads(taxon_map)["A"] == [f"A_{i}" for i in range(1, 12)]


def test_phylonet_nexus_defaults_and_override():
    nexus = to_phylonet_nexus(["((A:1,A:1):1,B:1);"])
    assert "Tree geneTree1 = ((A_1,A_2),B);" in nexus
    assert "InferNetwork_MP_Allopp (all) 1 -x 1 -pl 12 -di -a <A:A_1,A_2; B:B>;" in nexus
    nexus = to_phylonet_nexus(["(A:1,B:1);"], command="InferNetwork_MP_Allopp (all) 2 -x 5 -pl 4 -di;")
    assert "InferNetwork_MP_Allopp (all) 2 -x 5 -pl 4 -di -a <A:A; B:B>;" in nexus


def test_numbered_trees_accept_support_values():
    _, taxon_map = to_mpsugar_nexus(["((A:1,A:1)90:1,B:1);"])
    assert json.loads(taxon_map) == {"A": ["A_1", "A_2"], "B": ["B_1"]}
    assert "Tree geneTree1 = ((A_1,A_2),B);" in to_phylonet_nexus(["((A,A)0.9,B);"])


def test_strip_copy_suffixes():
    assert strip_copy_suffixes("((A_1,A_2),B);") == "((A,A),B);"


def test_strip_copy_suffixes_drops_annotations_and_lengths():
    out = strip_copy_suffixes("((A_1:0.1[&x=1],A_2:0.2):1,B:2);")
    assert out == "((A,A),B);"


def test_astral_command():
    cmd = astral_command("g.nwk", "s.nwk", "astral4", 8)
    assert cmd == ["astral4", "-i", "g.nwk", "-o", "s.nwk", "-t", "8"]


ALIGNMENT = {"A_0_1_0": "ACGT", "A_1_1_0": "ACGA", "B_1_0": "ACGG"}


def test_alignment_to_nexus():
    nex = alignment_to_nexus(ALIGNMENT)
    assert "NTAX=3 NCHAR=4" in nex
    assert "A_0_1_0  ACGT" in nex


def test_ragged_alignment_fails():
    with pytest.raises(ValueError):
        alignment_to_nexus({"A": "ACGT", "B": "AC"})


def test_strip_alloppnet_suffixes_only_polyploids():
    out = strip_alloppnet_suffixes("((A0,A1),(B1,C0));", {"A": 4, "B": 2, "C": 2})
    assert out == "((A,A),(B1,C0));"


def test_read_fasta(tmp_path):
    path = tmp_path / "g.fasta"
    path.write_text(">A_1_0 description\nAC\nGT\n>B_1_0\nACGA\n")
    assert read_fasta(path) == {"A_1_0": "ACGT", "B_1_0": "ACGA"}


def test_grampa_output_restores_original_names():
    from polybench.prep.trees import make_unique_labels, restore_labels

    mapping = make_unique_labels(["Qlo", "Qlo03"])
    raw = f"(({mapping['Qlo']}*,{mapping['Qlo']}+)<5>,Qlo03)<9>;"
    assert restore_labels(clean_grampa_tree(raw), mapping) == "((Qlo,Qlo),Qlo03);"


def test_clean_grampa_tree_grandma_labels():
    assert clean_grampa_tree("((A**,B)<P*>,(A+,C)<5+>)<1>;") == "((A,B),(A,C));"
