from pathlib import Path

from polybench.prep.formats import write_grampa_inputs


def prepare(dataset_dir, work_dir, params) -> None:
    write_grampa_inputs(dataset_dir, Path(work_dir) / "inputs")
