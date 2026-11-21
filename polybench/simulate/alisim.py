import random
import subprocess
import tempfile
from pathlib import Path

from ete3 import Tree

from polybench.simulate.fractionation import species_tree_leaf
from polybench.simulate.simphy import _pool, relabel_for_simphy, replicate_seed


def sample_gtr_params(gtr_pkl, seed: int) -> dict:
    # one set of GTR rates, base frequencies, gamma shape and alignment length from the empirical pool
    return random.Random(seed).choice(_pool(gtr_pkl))


def iqtree_model(p: dict) -> str:
    # alpha > 3 is treated as no rate heterogeneity
    rates = "/".join(f"{p[k]:.6f}" for k in ("AC", "AG", "AT", "CG", "CT", "GT"))
    freqs = "/".join(f"{p[k]:.6f}" for k in ("pi_A", "pi_C", "pi_G", "pi_T"))
    model = f"GTR{{{rates}}}+FU{{{freqs}}}"
    return model if p["alpha"] > 3.0 else f"{model}+G{{{p['alpha']:.6f}}}"


def alisim_command(gene_tree, out_prefix, model: str, length: int, seed: int, iqtree_bin: str = "iqtree3") -> list[str]:
    return [str(iqtree_bin), "--alisim", str(out_prefix), "-m", model, "-t", str(gene_tree),
            "--length", str(length), "-seed", str(seed), "-af", "fasta", "--quiet"]


def write_taxon_map(gene_trees: list[str], mapping: dict, path) -> None:
    # SimPhy-labelled sequence id -> restored taxon name (copies map to the same taxon)
    rows = {}
    for newick in gene_trees:
        for leaf in Tree(newick, format=1).iter_leaves():
            leaf_key = species_tree_leaf(leaf.name)
            if leaf_key not in mapping:
                raise ValueError(f"gene tree tip '{leaf.name}' is not a species tree leaf")
            rows[leaf.name] = mapping[leaf_key]
    Path(path).write_text("".join(f"{i}\t{t}\n" for i, t in sorted(rows.items())))


def simulate_alignments(gene_trees: list[str], mapping: dict, out_dir, gtr_pkl, seed: int,
                        iqtree_bin: str = "iqtree3") -> None:
    # one alignment per gene tree (SimPhy-labelled); one parameter set sampled per dataset and shared by all genes.
    # writes <out_dir>/gene_0001.fasta ..., taxon_map.tsv and alisim_params.txt
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    p = sample_gtr_params(gtr_pkl, seed)
    model, length = iqtree_model(p), int(p["alignment_length"])
    (out_dir / "alisim_params.txt").write_text(f"model={model}\nalignment_length={length}\nalpha={p['alpha']}\nseed={seed}\n")
    with tempfile.TemporaryDirectory() as tmp:
        for i, tree in enumerate(gene_trees, 1):
            tree_file = Path(tmp) / f"gene_{i:04d}.tree"
            tree_file.write_text(tree + "\n")
            prefix = out_dir / f"gene_{i:04d}"
            subprocess.run(alisim_command(tree_file, prefix, model, length, seed * 100_000 + i, iqtree_bin),
                           check=True)
            produced = prefix.with_suffix(".fa")
            if not produced.is_file():
                raise FileNotFoundError(f"AliSim did not write {produced}")
            produced.replace(prefix.with_suffix(".fasta"))
    write_taxon_map(gene_trees, mapping, out_dir / "taxon_map.tsv")


def simulate_replicate_alignments(data_root, config: str, network: str, gtr_pkl, seed: int, replicates: int,
                                  iqtree_bin: str = "iqtree3") -> None:
    # alignments/ for every replicate written by simulate_network, from its SimPhy-labelled gene trees
    multree = (Path(data_root) / "model_networks" / f"{network}.tre").read_text().strip()
    _, mapping = relabel_for_simphy(multree)
    for k in range(1, replicates + 1):
        rep = Path(data_root) / "simulated" / config / network / f"rep{k}"
        genes = [g for g in (rep / "simphy_gene_trees.nwk").read_text().splitlines() if g.strip()]
        simulate_alignments(genes, mapping, rep / "alignments", gtr_pkl, replicate_seed(seed, k), iqtree_bin)
