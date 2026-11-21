from polybench.prep.ploidy import (
    copy_number_distribution,
    infer_ploidy,
    smooth_copy_number,
    write_ploidy_tsv,
    write_polyphest_multiset,
)


def test_paper_example():
    assert smooth_copy_number({1: 3, 2: 2, 3: 2}) == 2


def test_two_distinct_values_use_most_frequent():
    assert smooth_copy_number({1: 5, 2: 3}) == 1
    assert smooth_copy_number({1: 3, 2: 5}) == 2


def test_two_distinct_tie_takes_smaller():
    assert smooth_copy_number({1: 4, 3: 4}) == 1


def test_single_value():
    assert smooth_copy_number({4: 7}) == 4


def test_kernel_prefers_neighbour_mass():
    assert smooth_copy_number({1: 3, 2: 3, 3: 3, 6: 4}) == 2


def test_distribution_counts_only_trees_with_taxon():
    trees = [
        "((a:1,a:1):1,b:1);",
        "((a:1,a:1):1,b:1);",
        "(a:1,b:1);",
        "(b:1,c:1);",
    ]
    assert copy_number_distribution(trees, "a") == {2: 2, 1: 1}
    assert copy_number_distribution(trees, "b") == {1: 4}


def test_infer_ploidy():
    trees = ["((a:1,a:1):1,b:1);"] * 3
    assert infer_ploidy(trees) == {"a": 2, "b": 1}


def test_genus_species_label_is_one_taxon():
    trees = ["((Ephedra_sinica:1,Ephedra_sinica:1):1,Ephedra_foeminea:1);"] * 2
    assert infer_ploidy(trees) == {"Ephedra_foeminea": 1, "Ephedra_sinica": 2}


def test_writers(tmp_path):
    write_ploidy_tsv({"b": 1, "a": 2}, tmp_path / "p.tsv")
    assert (tmp_path / "p.tsv").read_text() == "a\t2\nb\t1\n"
    write_polyphest_multiset({"b": 1, "a": 2}, tmp_path / "m.txt")
    assert (tmp_path / "m.txt").read_text() == "a\na\nb\n"
