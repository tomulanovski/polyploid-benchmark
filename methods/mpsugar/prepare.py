from pathlib import Path

from polybench.prep.formats import to_mpsugar_nexus


def prepare(dataset_dir, work_dir, params) -> None:
    # every leaf numbered <taxon>_<copy>, taxon map JSON {taxon: [copy labels]}
    inputs = Path(work_dir) / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    lines = [l for l in (Path(dataset_dir) / "gene_trees.nwk").read_text().splitlines() if l.strip()]
    nexus, taxon_map = to_mpsugar_nexus(lines)
    (inputs / "mpsugar_trees.nex").write_text(nexus)
    (inputs / "taxon_map.json").write_text(taxon_map)
