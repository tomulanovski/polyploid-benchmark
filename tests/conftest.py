from pathlib import Path

import pytest
from ete3 import Tree

CLUSTER = Path(__file__).parent / "fixtures" / "cluster"


def canonical(newick: str) -> str:
    # topology and leaf labels only, children sorted, so child order and branch lengths do not matter
    def walk(node):
        if node.is_leaf():
            return node.name
        return "(" + ",".join(sorted(walk(c) for c in node.children)) + ")"
    return walk(Tree(newick.strip(), format=1))


def cluster_dir(method: str) -> Path:
    d = CLUSTER / method
    if not d.is_dir():
        pytest.skip(f"no cluster fixture for {method}: copy dataset/, legacy_input/, legacy_output/ "
                    f"and legacy_result.tre from one finished cluster run into {d}")
    return d
