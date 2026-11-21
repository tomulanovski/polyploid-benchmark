import shutil
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest

from polybench.layout import Dataset, find_datasets, read_taxon_map, sequence_taxa
from polybench.validate import validate_dataset, validate_or_raise

FIXTURES = Path(__file__).parent / "fixtures" / "bad_datasets"
THREE = ("line 5: tree is not rooted", "ploidy.tsv has taxa not in the gene trees", "not single-labelled")


def ds(name):
    return Dataset("empirical", None, name, None, FIXTURES / name, FIXTURES)


def problems(name, requires=(), method=None):
    return validate_dataset(ds(name), list(requires), method)


def has(name, text, **kw):
    found = problems(name, **kw)
    assert any(text in p for p in found), found


def test_good_dataset_is_valid():
    assert problems("good", requires=["ploidy", "species_tree", "alignments"]) == []


def test_good_mapped_alignments_are_valid():
    assert problems("good_mapped", requires=["alignments"]) == []


def test_empty_gene_trees():
    has("empty_gene_trees", "gene_trees.nwk is empty")


def test_missing_gene_trees(tmp_path):
    dataset = Dataset("empirical", None, "x", None, tmp_path, tmp_path)
    assert any("gene_trees.nwk is missing" in p for p in validate_dataset(dataset, []))


@pytest.mark.filterwarnings("ignore:.*fewer than 10:UserWarning")  # small fixtures on purpose
def test_bad_newick_names_line():
    has("bad_newick", "line 5")


def test_unrooted_tree():
    has("unrooted", "line 4: tree is not rooted")


def test_single_leaf_tree():
    has("one_leaf", "line 3: fewer than 2 leaves")


def test_bad_label():
    has("bad_label", "invalid taxon label")


def test_one_taxon():
    has("one_taxon", "need at least 2")


def test_few_trees_warns_not_errors():
    with pytest.warns(UserWarning, match="fewer than 10"):
        found = problems("few_trees")
    assert found == []


def test_ploidy_columns():
    has("ploidy_columns", "2 tab-separated columns")


def test_ploidy_zero():
    has("ploidy_zero", "integer >= 1")


def test_ploidy_missing_taxon():
    has("ploidy_missing_taxon", "missing taxa: ['d']")


def test_ploidy_extra_taxon():
    has("ploidy_extra_taxon", "not in the gene trees: ['e']")


def test_species_tree_repeated_labels():
    has("species_repeated", "not single-labelled")


def test_species_tree_taxa_differ():
    has("species_taxa", "species_tree.nwk")


def test_alignments_missing_dir():
    has("alignments_missing", "no .fasta files", requires=["alignments"])


def test_alignments_not_checked_when_not_required():
    assert problems("alignments_missing") == []


def test_alignments_duplicate_id():
    has("alignments_duplicate_id", "duplicate sequence ids", requires=["alignments"])


def test_alignments_unknown_taxon():
    has("alignments_unknown_taxon", "['zzz']", requires=["alignments"])


def test_sequence_taxa_with_map():
    got = sequence_taxa(FIXTURES / "good_mapped" / "alignments")
    assert got == {"g1.fasta": {"a_1": "a", "a_2": "a", "b_1": "b", "c_1": "c", "d_1": "d"}}


def test_taxon_map_homeolog_column(tmp_path):
    (tmp_path / "g1.fasta").write_text(">x\nAC\n>y\nAG\n>z\nAT\n")
    (tmp_path / "taxon_map.tsv").write_text("x\tT\tT1:A\ny\tT\tT1:B\nz\tD\n")
    assert read_taxon_map(tmp_path) == ({"x": "T", "y": "T", "z": "D"}, {"x": ("T1", "A"), "y": ("T1", "B")})
    assert sequence_taxa(tmp_path) == {"g1.fasta": {"x": "T", "y": "T", "z": "D"}}
    (tmp_path / "taxon_map.tsv").write_text("x\tT\tT1:C\n")
    with pytest.raises(ValueError, match="individual:A"):
        read_taxon_map(tmp_path)


def test_alloppnet_rejects_given_ploidy_above_two():
    has("alloppnet_polyploid", "ploidy <= 2", method="alloppnet")


def test_alloppnet_rejects_inferred_ploidy_above_two():
    has("alloppnet_inferred", "ploidy <= 2", method="alloppnet")


def test_alloppnet_accepts_diploid():
    assert problems("good_mapped", method="alloppnet") == []


def test_polyphest_missing_ploidy_is_ok():
    assert problems("good_mapped", requires=["ploidy"], method="polyphest") == []


def test_grampa_prefix_overlap_warns():
    with pytest.warns(UserWarning, match="rename labels"):
        found = problems("grampa_prefix", method="grampa")
    assert found == []


def make_sim(tmp_path, network_text):
    root = tmp_path / "data"
    d = root / "simulated" / "conf" / "Net" / "rep1"
    shutil.copytree(FIXTURES / "good", d)
    (root / "model_networks").mkdir()
    if network_text is not None:
        (root / "model_networks" / "Net.tre").write_text(network_text)
    return find_datasets(root, "simulated/*/*/*")[0]


def test_simulated_true_network_ok(tmp_path):
    assert validate_dataset(make_sim(tmp_path, "((a,b),(c,(c,d)));"), []) == []


def test_simulated_true_network_missing(tmp_path):
    found = validate_dataset(make_sim(tmp_path, None), [])
    assert any("true network" in p and "missing" in p for p in found)


def test_simulated_true_network_taxa_differ(tmp_path):
    found = validate_dataset(make_sim(tmp_path, "((a,b),(c,e));"), [])
    assert any("true network" in p and "['e']" in p for p in found)


def test_three_problems_all_reported():
    found = problems("three_problems")
    assert len(found) == 3, found
    for text in THREE:
        assert any(text in p for p in found), found


def test_validate_or_raise_groups_by_dataset():
    spec = SimpleNamespace(name="padre", path=FIXTURES, requires=[], params={})
    with pytest.raises(ValueError) as err:
        validate_or_raise([ds("good"), ds("unrooted"), ds("ploidy_zero")], spec)
    msg = str(err.value)
    assert "empirical/unrooted:" in msg and "empirical/ploidy_zero:" in msg
    assert "empirical/good:" not in msg


def test_validate_or_raise_passes_on_valid():
    validate_or_raise([ds("good")], SimpleNamespace(name="padre", path=FIXTURES, requires=[], params={}))


def test_method_validate_hook(tmp_path):
    hook = tmp_path / "demo"
    hook.mkdir()
    (hook / "validate.py").write_text(
        "def extra_checks(dataset_dir, params):\n    return [f'demo says {params[\"k\"]} for {dataset_dir.name}']\n"
    )
    spec = SimpleNamespace(name="demo", path=hook, requires=[], params={"k": 7})
    with pytest.raises(ValueError, match="demo says 7 for good"):
        validate_or_raise([ds("good")], spec)
