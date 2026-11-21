from pathlib import Path


def prepare(dataset_dir, work_dir, params) -> None:
    # PADRE reads the gene trees as given, one Newick per line
    inputs = Path(work_dir) / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    lines = [l for l in (Path(dataset_dir) / "gene_trees.nwk").read_text().splitlines() if l.strip()]
    (inputs / "padre_trees.tre").write_text("".join(l + "\n" for l in lines))
