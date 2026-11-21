import json
import re
import sys
from pathlib import Path

from ete3 import Tree
from ete3.parser.newick import NewickError

from polybench.prep.formats import clean_grampa_tree, strip_alloppnet_suffixes
from polybench.prep.ploidy import read_ploidy_tsv
from polybench.prep.trees import restore_labels


def _one_newick(text: str, source) -> str:
    newick = text.strip()
    if not newick:
        raise ValueError(f"{source}: no tree found")
    if not newick.endswith(";"):
        newick += ";"
    if newick.count(";") != 1:
        raise ValueError(f"{source}: expected one Newick tree, found {newick.count(';')}")
    try:
        Tree(newick, format=1)
    except NewickError as err:
        raise ValueError(f"{source}: not valid Newick ({err})") from err
    return newick


def read_taxa_map(path) -> dict[str, str]:
    # original<TAB>label the tool saw, one line per taxon (written by write_grampa_inputs)
    mapping = {}
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        cols = line.split("\t")
        if len(cols) != 2 or not all(cols):
            raise ValueError(f"{path} line {n}: expected original<TAB>label")
        mapping[cols[0]] = cols[1]
    if not mapping:
        raise ValueError(f"{path}: empty taxa map")
    return mapping


def _grampa_tree(text: str, taxa_map, source) -> str:
    mapping = read_taxa_map(taxa_map)
    newick = restore_labels(_one_newick(clean_grampa_tree(text), source), mapping)
    unknown = sorted({leaf.name for leaf in Tree(newick, format=1)} - set(mapping))
    if unknown:
        raise ValueError(f"{source}: labels not among the input taxa: {unknown}")
    return newick


def parse_padre(path) -> str:
    # <input>-result.tre holds the MUL-tree
    return _one_newick(Path(path).read_text(), path)


def parse_polyphest(path) -> str:
    # <input>-polyphest.txt has lines "multree: ...", "network: ...", "Time taken: ..."
    lines = [l for l in Path(path).read_text().splitlines() if l.startswith("multree:")]
    if len(lines) != 1:
        raise ValueError(f"{path}: expected one 'multree:' line, found {len(lines)}")
    return _one_newick(lines[0].split(":", 1)[1], path)


def parse_grampa_scores(path, taxa_map) -> str:
    # grampa-scores.txt: header, then rows best first, labeled tree in the last tab column
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2 or not lines[1].strip():
        raise ValueError(f"{path}: no score rows")
    cols = lines[1].split("\t")
    if len(cols) < 5:
        raise ValueError(f"{path}: expected 5 tab-separated columns in the first row, got {len(cols)}")
    return _grampa_tree(cols[-1], taxa_map, path)


def parse_grandma_first_iteration(history_json, taxa_map) -> str:
    # iteration (0, 0) of a GRANDMA run is the GRAMPA result
    first = json.loads(Path(history_json).read_text()).get("(0, 0)")
    if not first or not first.get("best_mt"):
        raise ValueError(f"{history_json}: no best_mt for iteration (0, 0)")
    return _grampa_tree(first["best_mt"], taxa_map, history_json)


def parse_grandma_split(out_dir, taxa_map) -> str:
    path = Path(out_dir) / "final_multree.tre"
    return _grampa_tree(path.read_text(), taxa_map, path)


def parse_mpsugar(path) -> str:
    # run_mpsugar.py output, the first "Newick: ..." line is network 1 of chain 1 (as the paper took it)
    lines = [l for l in Path(path).read_text().splitlines() if l.startswith("Newick:")]
    if not lines:
        raise ValueError(f"{path}: no 'Newick:' line")
    return _one_newick(lines[0].split(":", 1)[1], path)


def parse_phylonet(path) -> str:
    # PhyloNet stdout, with -di every inferred network has a "Visualize in Dendroscope : <eNewick>" line
    lines = [l for l in Path(path).read_text().splitlines() if l.startswith("Visualize in Dendroscope")]
    if not lines:
        raise ValueError(f"{path}: no 'Visualize in Dendroscope' line, keep -di in the PhyloNet command")
    return _one_newick(lines[0].split(":", 1)[1], path)


def parse_alloppnet_mcc(path, ploidy_tsv) -> str:
    # TreeAnnotator NEXUS: optional Translate block, one tree with [&...] annotations.
    # AlloppNET names the two homeologs of a tetraploid <taxon>0 and <taxon>1.
    text = Path(path).read_text()
    trees = re.findall(r"^\s*tree\s+\S+\s*=\s*(.+)$", text, flags=re.M | re.I)
    if len(trees) != 1:
        raise ValueError(f"{path}: expected one tree, found {len(trees)}")
    tree = Tree(_one_newick(re.sub(r"\[[^\]]*\]", "", trees[0]), path), format=1)
    block = re.search(r"\btranslate\s+(.*?);", text, flags=re.S | re.I)
    if block:
        names = dict(entry.split() for entry in block.group(1).split(",") if entry.strip())
        for leaf in tree:
            leaf.name = names[leaf.name]
    copies = read_ploidy_tsv(ploidy_tsv)
    newick = strip_alloppnet_suffixes(tree.write(format=9), {t: 2 * k for t, k in copies.items()})
    unknown = sorted({leaf.name for leaf in Tree(newick, format=9)} - set(copies))
    if unknown:
        raise ValueError(f"{path}: labels not in {ploidy_tsv}: {unknown}")
    return newick


PARSERS = {
    "padre": parse_padre,
    "polyphest": parse_polyphest,
    "grampa_scores": parse_grampa_scores,
    "grandma_first_iteration": parse_grandma_first_iteration,
    "grandma_split": parse_grandma_split,
    "mpsugar": parse_mpsugar,
    "phylonet": parse_phylonet,
    "alloppnet_mcc": parse_alloppnet_mcc,
}

if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in PARSERS:
        sys.exit(f"usage: python -m polybench.prep.outputs {{{','.join(PARSERS)}}} <input> [taxa_map | ploidy.tsv]")
    sys.stdout.reconfigure(newline="\n")
    print(PARSERS[sys.argv[1]](*sys.argv[2:]))
