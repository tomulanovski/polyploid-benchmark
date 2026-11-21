import importlib.util
import re
import warnings
from pathlib import Path

from ete3 import Tree
from ete3.parser.newick import NewickError

from polybench.layout import Dataset, read_taxon_map, sequence_taxa, true_network_path
from polybench.prep.formats import homeolog_problems
from polybench.prep.ploidy import infer_ploidy
from polybench.prep.trees import make_unique_labels

METHODS_DIR = Path(__file__).resolve().parents[1] / "methods"
GRAMPA_METHODS = ("grampa", "grampa_iter")
BAD_LABEL = re.compile(r"[\s:,()\[\];]")


def _parse(text: str) -> Tree:
    return Tree(text, format=1)


def _read_gene_trees(path: Path, problems: list[str]) -> tuple[list[str], set[str]]:
    # returns the lines that parsed and the union of their leaf labels
    if not path.is_file():
        problems.append("gene_trees.nwk is missing")
        return [], set()
    lines = path.read_text().splitlines()
    if not any(l.strip() for l in lines):
        problems.append("gene_trees.nwk is empty")
        return [], set()
    parsed, taxa = [], set()
    for n, line in enumerate(lines, 1):
        if not line.strip():
            problems.append(f"gene_trees.nwk line {n}: empty line")
            continue
        try:
            tree = _parse(line)
        except NewickError as err:
            problems.append(f"gene_trees.nwk line {n}: not valid Newick ({err})")
            continue
        leaves = [l.name for l in tree.iter_leaves()]
        if len(tree.children) != 2:
            problems.append(f"gene_trees.nwk line {n}: tree is not rooted (root has {len(tree.children)} children), root the gene trees first")
        if len(leaves) < 2:
            problems.append(f"gene_trees.nwk line {n}: fewer than 2 leaves")
        bad = sorted({l for l in leaves if not l or BAD_LABEL.search(l)})
        if bad:
            problems.append(f"gene_trees.nwk line {n}: invalid taxon label {bad} (no whitespace or : , ( ) ; [ ], leaf label = taxon name)")
        parsed.append(line)
        taxa.update(l for l in leaves if l)
    return parsed, taxa


def _compare_taxa(what: str, taxa: set[str], expected: set[str], problems: list[str]) -> None:
    if taxa - expected:
        problems.append(f"{what} has taxa not in the gene trees: {sorted(taxa - expected)}")
    if expected - taxa:
        problems.append(f"{what} is missing taxa: {sorted(expected - taxa)}")


def _read_ploidy(path: Path, problems: list[str]) -> dict[str, int]:
    ploidy = {}
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) != 2:
            problems.append(f"ploidy.tsv line {n}: expected 2 tab-separated columns, got {len(cols)}")
        elif not cols[1].strip().isdigit() or int(cols[1]) < 1:
            problems.append(f"ploidy.tsv line {n}: copies must be an integer >= 1, got '{cols[1]}'")
        else:
            ploidy[cols[0]] = int(cols[1])
    return ploidy


def _tree_leaves(path: Path, what: str, problems: list[str]) -> list[str] | None:
    try:
        tree = _parse(path.read_text().strip())
    except NewickError as err:
        problems.append(f"{what} is not valid Newick ({err})")
        return None
    return [l.name for l in tree.iter_leaves()]


def validate_dataset(dataset: Dataset, requires: list[str], method: str | None = None) -> list[str]:
    d = dataset.path
    problems: list[str] = []
    lines, taxa = _read_gene_trees(d / "gene_trees.nwk", problems)

    if lines and len(taxa) < 2:
        problems.append(f"only {len(taxa)} taxon in the gene trees, need at least 2")
    if 0 < len(lines) < 10:
        warnings.warn(f"{dataset.id}: only {len(lines)} gene trees (fewer than 10)")

    ploidy = None
    if (d / "ploidy.tsv").is_file():
        ploidy = _read_ploidy(d / "ploidy.tsv", problems)
        if taxa:
            _compare_taxa("ploidy.tsv", set(ploidy), taxa, problems)

    if (d / "species_tree.nwk").is_file():
        names = _tree_leaves(d / "species_tree.nwk", "species_tree.nwk", problems)
        if names is not None:
            if len(names) != len(set(names)):
                repeated = sorted({n for n in names if names.count(n) > 1})
                problems.append(f"species_tree.nwk is not single-labelled, repeated labels {repeated}")
            if taxa:
                _compare_taxa("species_tree.nwk", set(names), taxa, problems)
    elif "species_tree" in requires:
        problems.append("species_tree.nwk is required but missing")

    if "alignments" in requires:
        try:
            seq_taxa = sequence_taxa(d / "alignments")
            for fname, mapping in seq_taxa.items():
                unknown = sorted(set(mapping.values()) - taxa)
                if taxa and unknown:
                    problems.append(f"alignments/{fname}: taxa not in the gene trees (map them in alignments/taxon_map.tsv): {unknown}")
            _, homeologs = read_taxon_map(d / "alignments")
            if homeologs and lines:
                copies = ploidy if ploidy is not None else infer_ploidy(lines)
                problems += [f"alignments/taxon_map.tsv: {p}" for p in homeolog_problems(seq_taxa, copies, homeologs)]
        except (ValueError, FileNotFoundError) as err:
            problems.append(f"alignments: {err}")

    # a missing ploidy.tsv is fine for methods that need ploidy, the runner infers it
    if method == "alloppnet" and lines:
        counts = ploidy if ploidy is not None else infer_ploidy(lines)
        high = sorted(t for t, k in counts.items() if k > 2)
        if high:
            problems.append(f"alloppnet only supports ploidy <= 2, taxa with higher ploidy: {high}")
    if method in GRAMPA_METHODS and taxa:
        try:
            renamed = make_unique_labels(sorted(taxa))
        except ValueError as err:
            problems.append(f"{method}: {err}")
        else:
            if renamed:
                warnings.warn(f"{dataset.id}: {method} will rename labels that are substrings of other labels: {sorted(renamed)}")

    if dataset.kind == "simulated":
        net = true_network_path(dataset.root, dataset)
        if not net.is_file():
            problems.append(f"true network {net} is missing")
        elif taxa:
            names = _tree_leaves(net, "true network", problems)
            if names is not None:
                _compare_taxa("true network", {n for n in names if not n.startswith("#")}, taxa, problems)
    return problems


def _extra_checks(dataset: Dataset, spec) -> list[str]:
    hook = spec.path / "validate.py"
    if not hook.is_file():
        return []
    mod_spec = importlib.util.spec_from_file_location(f"validate_{spec.name}", hook)
    module = importlib.util.module_from_spec(mod_spec)
    mod_spec.loader.exec_module(module)
    return module.extra_checks(dataset.path, spec.params)


def validate_or_raise(datasets: list[Dataset], spec=None) -> None:
    # spec needs .name, .path (method folder), .requires (list[str]) and .params (dict); None checks the data only
    blocks = []
    for ds in datasets:
        if spec is None:
            problems = validate_dataset(ds, [])
        else:
            problems = validate_dataset(ds, spec.requires, spec.name) + _extra_checks(ds, spec)
        if problems:
            blocks.append(f"{ds.id}:\n" + "\n".join(f"  - {p}" for p in problems))
    if blocks:
        raise ValueError("invalid input datasets\n" + "\n".join(blocks))
