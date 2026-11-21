import sys
from itertools import combinations
from pathlib import Path

from ete3 import Tree, TreeNode

from polybench.prep.ploidy import read_ploidy_tsv


def select_copies(tree: Tree, taxon: str, k: int) -> list[TreeNode]:
    leaves = [leaf for leaf in tree.iter_leaves() if leaf.name == taxon]
    if len(leaves) <= k:
        return leaves
    if k == 1:
        return [max(leaves, key=lambda leaf: leaf.get_distance(tree))]
    dist = {frozenset((a, b)): a.get_distance(b) for a, b in combinations(leaves, 2)}
    first = max(dist, key=dist.get)  # ties: first pair in leaf order
    selected = [leaf for leaf in leaves if leaf in first]
    while len(selected) < k:
        rest = [leaf for leaf in leaves if leaf not in selected]
        selected.append(
            max(rest, key=lambda leaf: min(dist[frozenset((leaf, s))] for s in selected))
        )
    return selected


def reduce_tree(newick: str, ploidy: dict[str, int]) -> str:
    tree = Tree(newick, format=1)
    taxa = {leaf.name for leaf in tree.iter_leaves()}
    if taxa - set(ploidy):
        raise ValueError(f"taxa without a ploidy: {sorted(taxa - set(ploidy))}")
    keep = []
    for taxon in sorted(taxa):
        keep.extend(select_copies(tree, taxon, ploidy[taxon]))
    tree.prune(keep, preserve_branch_length=True)
    return tree.write(format=1)


def reduce_file(trees_path, ploidy_path, out_path) -> None:
    ploidy = read_ploidy_tsv(ploidy_path)
    lines = [t for t in Path(trees_path).read_text(encoding="utf-8").splitlines() if t.strip()]
    if not lines:
        raise ValueError(f"{trees_path}: no gene trees")
    Path(out_path).write_text("".join(reduce_tree(t, ploidy) + "\n" for t in lines), encoding="utf-8", newline="\n")


# python -m polybench.prep.copy_selection <gene_trees.nwk> <ploidy.tsv> <out.nwk>
if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: python -m polybench.prep.copy_selection <gene_trees.nwk> <ploidy.tsv> <out.nwk>")
    reduce_file(*sys.argv[1:])
