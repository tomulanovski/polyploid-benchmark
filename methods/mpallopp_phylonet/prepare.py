from pathlib import Path

from polybench.prep.formats import to_phylonet_nexus


def prepare(dataset_dir, work_dir, params) -> None:
    # repeated taxa numbered <taxon>_<copy>, PHYLONET block with the command and the taxa map
    inputs = Path(work_dir) / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    lines = [l for l in (Path(dataset_dir) / "gene_trees.nwk").read_text().splitlines() if l.strip()]
    (inputs / "mpallopp.nex").write_text(to_phylonet_nexus(lines, params["command"]))
