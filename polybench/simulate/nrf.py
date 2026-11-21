from pathlib import Path

import pandas as pd
from ete3 import Tree

from polybench.layout import find_datasets
from polybench.simulate.fractionation import species_tree_leaf
from polybench.simulate.simphy import CONFIGS, rename_leaves


def nrf(gene_tree: str, species_tree: str) -> float:
    # unrooted Robinson-Foulds distance normalised by its maximum
    g, s = Tree(gene_tree, format=1), Tree(species_tree, format=1)
    if sorted(g.get_leaf_names()) != sorted(s.get_leaf_names()):
        raise ValueError("gene tree and species tree have different leaf sets")
    rf, max_rf = g.robinson_foulds(s, unrooted_trees=True)[:2]
    return rf / max_rf if max_rf else 0.0


def nrf_table(data_root, configs, out_csv=None) -> pd.DataFrame:
    # mean gene tree vs species tree nRF per replicate, only meaningful without duplication and loss
    rows = []
    for config in configs:
        if config not in CONFIGS:
            raise ValueError(f"unknown config {config}")
        if CONFIGS[config]["dup"] or CONFIGS[config]["loss"]:
            raise ValueError(f"nRF is only defined for configs without duplication and loss, not {config}")
        for ds in find_datasets(data_root, f"simulated/{config}/*/rep*"):
            species = (ds.path / "simphy_species_tree.nwk").read_text().strip()
            genes = [rename_leaves(g, species_tree_leaf)
                     for g in (ds.path / "simphy_gene_trees.nwk").read_text().splitlines()]
            rows.append({"config": config, "network": ds.network, "replicate": ds.replicate,
                         "mean_nrf": sum(nrf(g, species) for g in genes) / len(genes)})
    df = pd.DataFrame(rows, columns=["config", "network", "replicate", "mean_nrf"])
    if out_csv is not None:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
    return df
