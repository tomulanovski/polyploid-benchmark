# Regenerates the toy data: gene trees are perturbed copies of a model MUL-tree (NNI moves, branch length
# jitter, occasional loss of one copy of a polyploid). Not SimPhy, just small inputs to exercise the pipeline.
# Run from the repository root: python examples/toy/make_toy.py
import random
import shutil
from collections import Counter
from pathlib import Path

from ete3 import Tree

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(__file__).resolve().parent / "data"
SIMULATED = ["Shahrestani_2015", "Wisecaver_2023"]
EMPIRICAL = "Ding_2023"
N_TREES = 50


def load(network: str) -> Tree:
    t = Tree((ROOT / "tests" / "data" / "model_networks" / f"{network}.tre").read_text().strip(), format=1)
    for n in t.traverse():
        n.dist /= 1e8  # generations -> about 0.1 substitutions per site root to tip
    return t


def drop_unary_root(t: Tree) -> Tree:
    while len(t.children) == 1:
        t = t.children[0].detach()
    return t


def nni(t: Tree, rng: random.Random) -> None:
    # swap a child of an internal node with that node's sibling
    inner = [n for n in t.traverse() if not n.is_leaf() and not n.is_root()]
    n = rng.choice(inner)
    sib = n.get_sisters()[0]
    child = rng.choice(n.children)
    parent = n.up
    child.detach()
    sib.detach()
    n.add_child(sib)
    parent.add_child(child)


def gene_tree(multree: Tree, rng: random.Random) -> str:
    t = multree.copy()
    for _ in range(rng.choice([0, 0, 1, 1, 2])):
        nni(t, rng)
    copies = Counter(t.get_leaf_names())
    polyploids = sorted(k for k, c in copies.items() if c > 1)
    if polyploids and rng.random() < 0.2:
        lost = rng.choice(polyploids)
        rng.choice([l for l in t.iter_leaves() if l.name == lost]).delete(preserve_branch_length=True)
        t = drop_unary_root(t)
    for n in t.traverse():
        n.dist = round(n.dist * rng.lognormvariate(0, 0.3), 6)
    return t.write(format=5)


def species_tree(multree: Tree) -> str:
    # single-labelled: keep the first copy of every taxon, as an ASTRAL tree would have one leaf per taxon
    t = multree.copy()
    seen = set()
    for leaf in list(t.iter_leaves()):
        if leaf.name in seen:
            leaf.delete(preserve_branch_length=True)
        seen.add(leaf.name)
    return drop_unary_root(t).write(format=9)


def write_dataset(out: Path, multree: Tree, rng: random.Random) -> None:
    out.mkdir(parents=True)
    (out / "gene_trees.nwk").write_text("".join(gene_tree(multree, rng) + "\n" for _ in range(N_TREES)))
    (out / "species_tree.nwk").write_text(species_tree(multree) + "\n")


def main() -> None:
    rng = random.Random(2026)
    for sub in ("model_networks", "simulated", "empirical"):
        if (DATA / sub).exists():
            shutil.rmtree(DATA / sub)
    (DATA / "model_networks").mkdir(parents=True)
    for net in SIMULATED:
        shutil.copy(ROOT / "tests" / "data" / "model_networks" / f"{net}.tre", DATA / "model_networks")
        write_dataset(DATA / "simulated" / "toy" / net / "rep1", load(net), rng)
    multree = load(EMPIRICAL)
    out = DATA / "empirical" / "toy_empirical"
    write_dataset(out, multree, rng)
    copies = Counter(multree.get_leaf_names())
    (out / "ploidy.tsv").write_text("".join(f"{k}\t{copies[k]}\n" for k in sorted(copies)))
    (out / "published_network.nwk").write_text(multree.write(format=9) + "\n")


if __name__ == "__main__":
    main()
