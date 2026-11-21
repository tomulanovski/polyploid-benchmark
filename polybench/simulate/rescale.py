import random

from ete3 import Tree


def rescale(newick: str, height: float = 10_000_000, seed: int = 0) -> str:
    # random polytomy resolution, then node heights spread evenly by depth so every tip sits at `height`
    random.seed(seed)
    tree = Tree(newick, format=1)
    tree.resolve_polytomy(recursive=True, default_dist=1)
    depth = {}
    for node in tree.traverse("postorder"):
        depth[node] = 0 if node.is_leaf() else 1 + max(depth[c] for c in node.children)
    if depth[tree] == 0:
        raise ValueError("tree has a single node")
    node_height = {n: height * (depth[tree] - d) / depth[tree] for n, d in depth.items()}
    for node in tree.traverse():
        if not node.is_root() and not node.is_leaf():
            node.dist = node_height[node] - node_height[node.up]
    for leaf in tree.iter_leaves():
        above = 0.0
        node = leaf.up
        while not node.is_root():
            above += node.dist
            node = node.up
        leaf.dist = max(0.0, round(height - above, 12))
    return _newick(tree) + ";"


def _newick(node) -> str:
    s = node.name if node.is_leaf() else "(" + ",".join(_newick(c) for c in node.children) + ")"
    return s if node.is_root() else f"{s}:{node.dist:.12f}"


def to_nexus(newick: str) -> str:
    return f"#NEXUS\nbegin trees;\ntree 1={newick}\nend;\n"
