import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

from polybench.measures.compare_reticulations import pairwise_compare
from polybench.measures.reticulate_tree import ReticulateTree
from polybench.paper.common import save, write_csv

# Fig 3a: polyploid p is a hybrid of q and r
TRUTH = "(((q,p),(r,p)),(s,(t,u)));"
ERRORS = {  # Table 3 column -> (Fig 3 panel title, inferred network)
    "Wrong polyploid": ("(b) Wrong polyploid species", "(((q,u),(r,u)),(s,(t,p)));"),
    "Wrong parental lineage": ("(c) Wrong parental lineage", "(((q,p),(s,p)),(r,(t,u)));"),
    "Extra reticulation": ("(d) Extra reticulation", "(((q,p),(r,p)),((s,u),(t,u)));"),
}
MEASURES = {"ret_leaf_jaccard": "Reticulation descendants", "ret_sisters_jaccard": "Reticulation sisters",
            "num_rets_bias": "Difference in reticulation count (signed)", "mu_distance": "mu-distance"}


def table3() -> pd.DataFrame:
    truth = ReticulateTree(TRUTH, refold=True)
    cols = {}
    for col, (_, newick) in ERRORS.items():
        m = pairwise_compare(truth, ReticulateTree(newick, refold=True))
        if not m["mu_distance"] > 0:
            raise ValueError(f"{col}: mu-distance is {m['mu_distance']}, the paper says it is nonzero in every case")
        cols[col] = [m[k]["dist"] if isinstance(m[k], dict) else m[k] for k in MEASURES]
    return pd.DataFrame({"measure": list(MEASURES.values()), **cols})


def draw(tree: ReticulateTree, ax, title: str) -> None:
    g = tree.dag
    for node, d in ReticulateTree.compute_depths(g).items():
        g.nodes[node]["layer"] = d
    pos = nx.multipartite_layout(g, subset_key="layer")
    retic = [n for n in g if ReticulateTree.is_reticulation_node(g, n)]
    leaves = [n for n in g if g.out_degree(n) == 0 and n not in retic]
    inner = [n for n in g if n not in retic and n not in leaves]
    for nodes, shape, color in ((retic, "D", "orange"), (leaves, "o", "lightgreen"), (inner, "o", "lightblue")):
        nx.draw_networkx_nodes(g, pos, nodelist=nodes, node_shape=shape, node_color=color, node_size=600, ax=ax)
    nx.draw_networkx_edges(g, pos, arrows=True, ax=ax)
    nx.draw_networkx_labels(g, pos, labels={n: g.nodes[n].get("label", "") or "" for n in g}, font_size=9, ax=ax)
    ax.set_title(title, fontsize=11)
    ax.set_axis_off()


def make(summary_dir, networks_dir, out_dir) -> list:
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.6), constrained_layout=True)
    draw(ReticulateTree(TRUTH), axes[0], "(a) True network")
    for ax, (title, newick) in zip(axes[1:], ERRORS.values()):
        draw(ReticulateTree(newick), ax, title)
    return save(fig, out_dir, "fig3") + [write_csv(table3(), out_dir, "table3.csv")]
