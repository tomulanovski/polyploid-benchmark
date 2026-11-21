import argparse
import importlib.util
import json
import shutil
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass, replace
from pathlib import Path

import yaml

from polybench.config import load_config
from polybench.executors import (parse_slurm_time, render_slurm, run_command, run_with_timeout,
                                 time_budget, tool_env)
from polybench.layout import Dataset, find_datasets, result_dir
from polybench.prep.params import BOOL_FLAGS
from polybench.prep.ploidy import infer_ploidy, write_ploidy_tsv
from polybench.prep.species_tree import run_astral
from polybench.validate import METHODS_DIR, validate_or_raise

REQUIRES = ("gene_trees", "ploidy", "species_tree", "alignments")
KEYS = ("name", "version", "env", "requires", "output", "partial_match", "resources", "params")


@dataclass
class MethodSpec:
    name: str
    version: str
    env: str | None
    requires: list[str]
    output: str
    partial_match: bool
    resources: dict
    params: dict
    path: Path
    bool_flags: str = "value"


def load_method(name_or_path) -> MethodSpec:
    path = Path(name_or_path)
    if not path.is_dir():
        path = METHODS_DIR / str(name_or_path)
    if not (path / "method.yaml").is_file():
        raise FileNotFoundError(f"no method.yaml for '{name_or_path}' (looked in {path})")
    raw = yaml.safe_load((path / "method.yaml").read_text())
    missing = [k for k in KEYS if k not in raw]
    if missing:
        raise ValueError(f"{path / 'method.yaml'}: missing keys {missing}")
    if set(raw["requires"]) - set(REQUIRES):
        raise ValueError(f"{path / 'method.yaml'}: requires must be a subset of {list(REQUIRES)}")
    if raw["output"] not in ("multree", "network"):
        raise ValueError(f"{path / 'method.yaml'}: output must be 'multree' or 'network'")
    if not (path / "run.sh").is_file():
        raise FileNotFoundError(f"{path / 'run.sh'} not found")
    bool_flags = raw.get("bool_flags", "value")
    if bool_flags not in BOOL_FLAGS:
        raise ValueError(f"{path / 'method.yaml'}: bool_flags must be one of {list(BOOL_FLAGS)}")
    return MethodSpec(**{k: raw[k] for k in KEYS}, path=path, bool_flags=bool_flags)


def prepare_inputs(spec: MethodSpec, dataset: Dataset, cfg: dict, work_dir: Path) -> Path:
    # copies the gene trees and the inputs the method requires to <work_dir>/dataset and adds derived inputs
    # there; other files (the truth of a simulation, a published network) never reach the method
    work_ds = work_dir / "dataset"
    if work_ds.exists():
        shutil.rmtree(work_ds)
    work_ds.mkdir(parents=True)
    shutil.copy(dataset.path / "gene_trees.nwk", work_ds)
    for req, name in (("ploidy", "ploidy.tsv"), ("species_tree", "species_tree.nwk")):
        if req in spec.requires and (dataset.path / name).is_file():
            shutil.copy(dataset.path / name, work_ds)
    if "alignments" in spec.requires:
        shutil.copytree(dataset.path / "alignments", work_ds / "alignments")
    if "ploidy" in spec.requires and not (work_ds / "ploidy.tsv").is_file():
        trees = (work_ds / "gene_trees.nwk").read_text().splitlines()
        write_ploidy_tsv(infer_ploidy([t for t in trees if t.strip()]), work_ds / "ploidy.tsv")
    if "species_tree" in spec.requires and not (work_ds / "species_tree.nwk").is_file():
        if cfg["tools"].get("astral") is None:
            raise ValueError(f"{dataset.id}: species_tree.nwk is missing and tools.astral is not set in the config")
        run_astral(work_ds / "gene_trees.nwk", work_ds / "species_tree.nwk", cfg["tools"]["astral"],
                   spec.resources["threads"])
    return work_ds


def method_prepare(spec: MethodSpec, work_ds: Path, work_dir: Path, params: dict) -> None:
    hook = spec.path / "prepare.py"
    if not hook.is_file():
        return
    mod_spec = importlib.util.spec_from_file_location(f"prepare_{spec.name}", hook)
    module = importlib.util.module_from_spec(mod_spec)
    mod_spec.loader.exec_module(module)
    module.prepare(work_ds, work_dir, params)


def run_one(spec: MethodSpec, dataset: Dataset, out_dir, params: dict, label: str, timeout_s: int, cfg: dict) -> dict:
    t0 = time.monotonic()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in ("network.nwk", "status.json", "log.txt"):
        (out_dir / stale).unlink(missing_ok=True)
    params = {**spec.params, **params}
    (out_dir / "params.json").write_text(json.dumps(params))
    work_dir = out_dir / "work"
    if work_dir.exists():
        shutil.rmtree(work_dir)  # tools such as GRANDMA resume from state left by an earlier run
    runtime = 0.0
    try:
        work_ds = prepare_inputs(spec, dataset, cfg, work_dir)
        method_prepare(spec, work_ds, work_dir, params)
        args = [work_ds.as_posix(), out_dir.as_posix(), (out_dir / "params.json").as_posix(), str(spec.resources["threads"])]
        cmd = run_command(spec, cfg, args)
        code, runtime = run_with_timeout(cmd, tool_env(cfg), out_dir / "log.txt", time_budget(timeout_s, time.monotonic() - t0))
    except Exception:
        with open(out_dir / "log.txt", "a") as log:
            log.write(traceback.format_exc())
        state = "failed"
    else:
        network = out_dir / "network.nwk"
        if code is None:
            state = "timeout"
        elif code == 0 and network.is_file() and network.stat().st_size > 0:
            state = "completed"
        else:
            state = "failed"
    status = {"status": state, "runtime_s": runtime, "method": spec.name, "version": spec.version, "params": params}
    (out_dir / "status.json").write_text(json.dumps(status, indent=2))
    return status


def raise_if_not_completed(results: dict) -> None:
    # results: dataset id -> status dict, as returned by a local submit
    bad = {i: s["status"] for i, s in results.items() if s["status"] != "completed"}
    for i, st in bad.items():
        print(f"{i}: {st}", file=sys.stderr)
    if bad:
        raise RuntimeError(f"{len(bad)} of {len(results)} datasets did not complete")


def check_inputs(spec: MethodSpec, datasets: list[Dataset], params: dict, cfg: dict) -> None:
    # species_tree.nwk is inferred by the runner when absent, so it is not a validation requirement
    validate_or_raise(datasets, replace(spec, requires=[r for r in spec.requires if r != "species_tree"],
                                        params={**spec.params, **params}))
    if "species_tree" in spec.requires and cfg["tools"].get("astral") is None:
        absent = [d.id for d in datasets if not (d.path / "species_tree.nwk").is_file()]
        if absent:
            raise ValueError(f"species_tree.nwk missing for {absent} and tools.astral is not set in the config")


def submit(spec: MethodSpec, datasets: list[Dataset], label: str, params: dict, cfg: dict, cfg_path,
           dry_run: bool = False):
    check_inputs(spec, datasets, params, cfg)
    wrong = [d.id for d in datasets if Path(d.root) != Path(cfg["data_root"])]
    if wrong:
        raise ValueError(f"datasets {wrong} are not under data_root {cfg['data_root']}")
    label_dir = (Path(cfg["results_root"]) / label).resolve()
    timeout_s = parse_slurm_time(spec.resources["time"])
    if cfg["executor"]["type"] == "local":
        return {d.id: run_one(spec, d, result_dir(cfg["results_root"], label, d), params, label, timeout_s, cfg)
                for d in datasets}

    list_path, params_path, script_path = label_dir / "dataset_list.txt", label_dir / "params.json", label_dir / "submit.sh"
    script = render_slurm(spec, cfg, len(datasets), label_dir, spec.path.resolve(), params_path, list_path,
                          Path(cfg_path).resolve(), label)  # checks label_dir before anything is written
    (label_dir / "logs").mkdir(parents=True, exist_ok=True)
    list_path.write_text("".join(d.id + "\n" for d in datasets))
    params_path.write_text(json.dumps(params))
    script_path.write_text(script)
    if not dry_run:
        subprocess.run(["sbatch", str(script_path)], check=True)
    return script_path


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="polybench.runner")
    ap.add_argument("--one", action="store_true", required=True)
    ap.add_argument("--method", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--dataset-list", required=True)
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--config", required=True)
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    spec = load_method(a.method)
    ids = Path(a.dataset_list).read_text().split()
    dataset = find_datasets(cfg["data_root"], ids[a.index])[0]
    run_one(spec, dataset, result_dir(cfg["results_root"], a.label, dataset), json.loads(Path(a.params).read_text()),
            a.label, parse_slurm_time(spec.resources["time"]), cfg)


if __name__ == "__main__":
    main()
