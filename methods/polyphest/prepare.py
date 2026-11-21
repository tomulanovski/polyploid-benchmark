from pathlib import Path

from polybench.prep.ploidy import read_ploidy_tsv, write_polyphest_multiset


def prepare(dataset_dir, work_dir, params) -> None:
    inputs = Path(work_dir) / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    lines = [l for l in (Path(dataset_dir) / "gene_trees.nwk").read_text().splitlines() if l.strip()]
    (inputs / "polyphest_trees.tre").write_text("".join(l + "\n" for l in lines))
    write_polyphest_multiset(read_ploidy_tsv(Path(dataset_dir) / "ploidy.tsv"), inputs / "multi_set.txt")
