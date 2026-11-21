from collections import defaultdict

import numpy as np
from ete3 import Tree

from polybench.measures.reticulate_tree import ReticulateTree


def species_tree_leaf(gene_tip: str) -> str:
    # SimPhy tips are <species-tree leaf>_<locus>_<individual>
    return gene_tip.rsplit("_", 2)[0]


def split_copy(leaf: str):
    # "X_2" -> ("X", 2); labels never contain underscores, so any underscore is the copy suffix
    if "_" in leaf:
        base, idx = leaf.rsplit("_", 1)
        return base, int(idx)
    return leaf, None


def _iso_leaf_map(a, b, canonical):
    if a.is_leaf():
        return {a.uniq: b.uniq}
    mapping = {}
    for ca, cb in zip(sorted(a.children, key=lambda n: canonical[n]), sorted(b.children, key=lambda n: canonical[n])):
        mapping.update(_iso_leaf_map(ca, cb, canonical))
    return mapping


def build_events(species_tree) -> list[dict]:
    # Fold the uniquely suffixed species tree into WGD events (Holm grouping, as in ReticulateTree.multree_to_dag).
    # Inner events are mirrored onto the other blocks of their enclosing event so fractionation applies
    # whichever outer block survives. Mirrored events share the group (hence the rate) of their source.
    t = (species_tree if isinstance(species_tree, Tree) else Tree(species_tree, format=1)).copy()
    names = [leaf.name for leaf in t.iter_leaves()]
    if len(names) != len(set(names)):
        dupes = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(f"species tree has duplicate leaf labels, copies must be suffixed (X_1, X_2): {dupes[:5]}")
    for leaf in t.iter_leaves():
        leaf.add_feature("uniq", leaf.name)
        leaf.name = split_copy(leaf.name)[0]

    height_map = ReticulateTree._compute_height_map(t)
    canonical = ReticulateTree._compute_canonical_forms(t)

    removed = set()
    events = []
    group = 0
    for h in sorted(height_map, reverse=True):
        by_form = defaultdict(list)
        for n in height_map[h]:
            if id(n) not in removed:
                by_form[canonical[n]].append(n)
        for members in by_form.values():
            if len(members) <= 1:
                continue
            members.sort(key=lambda m: sorted(leaf.uniq for leaf in m.iter_leaves()))
            parents = [m.up for m in members]
            events.append({
                "sides": [frozenset(leaf.uniq for leaf in m.iter_leaves()) for m in members],
                "autopolyploidy": any(parents.count(p) > 1 for p in parents),
                "height": h, "group": group, "_nodes": members,
            })
            group += 1
            for iso in members[1:]:
                removed.update(id(d) for d in iso.traverse())

    for inner in list(events):
        inner_leaves = frozenset().union(*inner["sides"])
        enclosing = [e for e in events if e is not inner and e["height"] > inner["height"]
                     and inner_leaves <= e["sides"][0]]
        if not enclosing:
            continue
        outer = min(enclosing, key=lambda e: len(e["sides"][0]))
        for other in outer["_nodes"][1:]:
            imap = _iso_leaf_map(outer["_nodes"][0], other, canonical)
            if not all(leaf in imap for side in inner["sides"] for leaf in side):
                continue
            events.append({"sides": [frozenset(imap[leaf] for leaf in side) for side in inner["sides"]],
                           "autopolyploidy": inner["autopolyploidy"], "height": inner["height"],
                           "group": inner["group"], "_nodes": None})
    for e in events:
        del e["_nodes"]
    return events


def plan_removals(events, q_by_group, present_leaves, rng) -> set:
    # Per event and species, keep one random copy block and drop each other block with probability 1 - q.
    removed = set()
    for event in events:
        q = q_by_group[event["group"]]
        per_species = defaultdict(list)
        for side in event["sides"]:
            by_species = defaultdict(set)
            for leaf in side:
                by_species[split_copy(leaf)[0]].add(leaf)
            for species, block in by_species.items():
                per_species[species].append(block)
        for blocks in per_species.values():
            present = [{leaf for leaf in b if leaf in present_leaves and leaf not in removed} for b in blocks]
            present = [b for b in present if b]
            if len(present) < 2:
                continue
            order = list(rng.permutation(len(present)))
            for j in order[1:]:
                if rng.random() < 1.0 - q:
                    removed.update(present[j])
    return removed


def prune_gene_tree(newick: str, removed_leaves) -> str:
    # prune by node, not name: duplications can give repeated tip labels
    t = Tree(newick, format=1)
    keep = [leaf for leaf in t.iter_leaves() if species_tree_leaf(leaf.name) not in removed_leaves]
    t.prune(keep, preserve_branch_length=True)
    return t.write(format=5)


def apply_fractionation(gene_trees: list[str], species_tree: str, retention: float, seed: int) -> list[str]:
    # gene_trees carry SimPhy labels (<leaf>_<locus>_<individual>), species_tree the suffixed species tree
    if not 0.0 <= retention <= 1.0:
        raise ValueError(f"retention must be in [0, 1], got {retention}")
    events = build_events(species_tree)
    q_by_group = [retention] * (1 + max((e["group"] for e in events), default=-1))
    rng = np.random.default_rng(seed)
    out = []
    for newick in gene_trees:
        present = {species_tree_leaf(leaf.name) for leaf in Tree(newick, format=1).iter_leaves()}
        out.append(prune_gene_tree(newick, plan_removals(events, q_by_group, present, rng)))
    return out
