#!/usr/bin/env python3
# Lists the empirical datasets with a taxon of inferred copy number >= 3 and, for each, how many gene trees keep
# other copies under the fixed copy selection (polybench.prep.copy_selection) than under the legacy one.
# Usage: affected_by_copy_fix.py <papers_dir> [--trees PATTERN|MAP.yaml] [--ploidy-json PATTERN|MAP.yaml]
#        [--taxon label|before-last-underscore]
import argparse
import json
import sys
from pathlib import Path

import yaml
from ete3 import Tree

from polybench.prep.copy_selection import select_copies
from polybench.prep.ploidy import infer_ploidy
from polybench.prep.trees import split_newick


# LEGACY, only for this comparison: select_most_diverse_leaves of gene2net scripts/reduce_diverse_k.py, bug kept.
# Distances are stored under (a, b) in leaf order but looked up under the id()-ordered pair, so a lookup can miss
# and the leaf then counts as infinitely far from the selected ones. `order` stands in for id() in the tests.
def legacy_select(tree, leaves, k, order=id):
    if len(leaves) <= k:
        return leaves
    if k == 1:
        return [max(leaves, key=lambda leaf: leaf.get_distance(tree))]
    dist = {(a, b): a.get_distance(b) for i, a in enumerate(leaves) for b in leaves[i + 1:]}
    selected = list(max(dist.items(), key=lambda x: x[1])[0])
    while len(selected) < k:
        best, best_d = None, -1
        for leaf in leaves:
            if leaf in selected:
                continue
            d = float("inf")
            for s in selected:
                key = (leaf, s) if order(leaf) < order(s) else (s, leaf)
                if key in dist:
                    d = min(d, dist[key])
            if d > best_d:
                best, best_d = leaf, d
        selected.append(best)
    return selected


def tree_changes(tree, ploidy: dict[str, int], order=id) -> bool:
    for taxon, k in ploidy.items():
        leaves = [leaf for leaf in tree.iter_leaves() if leaf.name == taxon]
        if k < 3 or len(leaves) <= k:
            continue
        if set(select_copies(tree, taxon, k)) != set(legacy_select(tree, leaves, k, order)):
            return True
    return False


def legacy_taxon(name: str, taxa) -> str:
    # LEGACY extract_taxon_name of reduce_diverse_k.py: longest taxon of the JSON that the label starts with,
    # else everything before the last underscore
    for taxon in sorted(taxa, key=len, reverse=True):
        if name.startswith(taxon + "_") or name == taxon:
            return taxon
    return name.rsplit("_", 1)[0] if "_" in name else name


def read_trees(path: Path, rename) -> list[Tree]:
    trees = []
    for i, text in enumerate(split_newick(path.read_text(encoding="utf-8")), 1):
        try:
            tree = Tree(text, format=1)
        except Exception as e:
            raise ValueError(f"{path} tree {i}: {e}") from e
        for leaf in tree.iter_leaves():
            leaf.name = rename(leaf.name)
        trees.append(tree)
    if not trees:
        raise ValueError(f"{path}: no trees")
    return trees


def resolver(spec: str):
    # a YAML map {dataset: path} or a pattern with {dataset}; paths are relative to the dataset folder
    if spec.endswith((".yaml", ".yml")):
        paths = yaml.safe_load(Path(spec).read_text())
        return lambda d: d / paths[d.name] if d.name in paths else None
    return lambda d: d / spec.format(dataset=d.name)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Empirical datasets affected by the copy selection fix")
    ap.add_argument("papers_dir", type=Path, help="one folder per dataset")
    ap.add_argument("--trees", default="gene_trees.nwk", help="gene trees before copy selection: a path inside "
                    "each dataset folder, may contain {dataset}, or a YAML map dataset -> path")
    ap.add_argument("--ploidy-json", help="legacy taxa_copies.json per dataset (pattern or YAML map as --trees): "
                    "copy numbers and leaf -> taxon as reduce_diverse_k.py, instead of inferring them")
    ap.add_argument("--taxon", choices=["label", "before-last-underscore"], default="label",
                    help="without --ploidy-json, leaf label -> taxon: the label, or everything before the last "
                         "underscore")
    a = ap.parse_args(argv)
    trees_of = resolver(a.trees)
    json_of = resolver(a.ploidy_json) if a.ploidy_json else None
    datasets = sorted(p for p in a.papers_dir.iterdir() if p.is_dir())
    if not datasets:
        raise FileNotFoundError(f"no dataset folders in {a.papers_dir}")
    print("The legacy selection depends on object addresses (id()), so its counts can vary between runs; "
          "a dataset with any taxon of 3+ copies should be rerun.")
    without = []
    for d in datasets:
        path = trees_of(d)
        if path is None or not path.is_file():
            without.append(d.name)
            continue
        if json_of:
            jpath = json_of(d)
            if jpath is None or not jpath.is_file():
                raise FileNotFoundError(f"{d.name}: no taxa_copies JSON ({jpath})")
            copies = {t: int(k) for t, k in json.loads(jpath.read_text()).items()}
            trees = read_trees(path, lambda name: legacy_taxon(name, copies))
            taxa = {leaf.name for t in trees for leaf in t.iter_leaves()}
            ploidy = {t: copies.get(t, 1) for t in taxa}  # legacy default: keep 1 copy of an unlisted taxon
        else:
            last = a.taxon == "before-last-underscore"
            trees = read_trees(path, lambda name: name.rsplit("_", 1)[0] if last else name)
            ploidy = infer_ploidy([t.write(format=1) for t in trees])
        high = {t: k for t, k in ploidy.items() if k >= 3}
        if not high:
            continue
        changed = sum(tree_changes(t, high) for t in trees)
        print(f"{d.name}: taxa with >= 3 copies {high}; {changed} of {len(trees)} gene trees change")
    if without:
        print(f"no gene trees ({a.trees}) in: {', '.join(without)}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError) as e:
        print(f"affected_by_copy_fix: error: {e}", file=sys.stderr)
        sys.exit(1)
