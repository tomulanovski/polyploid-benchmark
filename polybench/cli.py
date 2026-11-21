import argparse
import subprocess
import sys
from pathlib import Path

import yaml

from polybench.config import load_config
from polybench.evaluate.compare import write_summary
from polybench.layout import find_datasets
from polybench.paper import NAMES, make_all
from polybench.runner import (MethodSpec, check_inputs, load_method, prepare_inputs,
                              method_prepare, raise_if_not_completed, submit)
from polybench.simulate.alisim import simulate_replicate_alignments
from polybench.simulate.nrf import nrf_table
from polybench.simulate.simphy import CONFIGS, simulate_network
from polybench.validate import validate_or_raise

# prep without --method: the inputs every method can ask for
SHARED = MethodSpec(name="shared", version="", env=None, requires=["ploidy", "species_tree"], output="network",
                    partial_match=False, resources={"threads": 1}, params={}, path=Path("."))


def datasets(cfg: dict, patterns: list[str]):
    found = {d.id: d for p in patterns for d in find_datasets(cfg["data_root"], p)}
    return [found[i] for i in sorted(found)]


def read_params(path) -> dict:
    if path is None:
        return {}
    params = yaml.safe_load(Path(path).read_text())
    if not isinstance(params, dict):
        raise ValueError(f"{path}: params file must be a mapping")
    return params


def cmd_validate(a):
    cfg = load_config(a.config)
    ds = datasets(cfg, a.datasets)
    if a.method:
        check_inputs(load_method(a.method), ds, read_params(a.params), cfg)
    else:
        validate_or_raise(ds)
    print(f"{len(ds)} datasets valid")


def cmd_prep(a):
    cfg = load_config(a.config)
    out = Path(a.out).resolve()
    if out.is_relative_to(Path(cfg["data_root"]).resolve()):
        raise ValueError(f"--out {a.out} is inside data_root, prep never writes into the data folder")
    ds = datasets(cfg, a.datasets)
    params = read_params(a.params)
    if a.method:
        spec = load_method(a.method)
        check_inputs(spec, ds, params, cfg)
    else:
        spec = SHARED
        validate_or_raise(ds)
    for d in ds:
        work_ds = prepare_inputs(spec, d, cfg, out / d.id)
        if a.method:
            method_prepare(spec, work_ds, out / d.id, {**spec.params, **params})
        print(f"{d.id}: {out / d.id}")


def cmd_run(a):
    cfg = load_config(a.config)
    if a.executor:
        cfg["executor"]["type"] = a.executor
    spec = load_method(a.method)
    label = a.label or spec.name
    result = submit(spec, datasets(cfg, a.datasets), label, read_params(a.params), cfg, a.config,
                    dry_run=a.dry_run)
    if cfg["executor"]["type"] == "local":
        raise_if_not_completed(result)
        print(f"{len(result)} datasets completed, results in {Path(cfg['results_root']) / label}")
    else:
        print(f"sbatch script: {result}" + (" (not submitted, --dry-run)" if a.dry_run else ""))


def cmd_evaluate(a):
    cfg = load_config(a.config)
    nrf = nrf_table(cfg["data_root"], a.nrf) if a.nrf else None  # before the slow part, so bad configs fail early
    out, errors = write_summary(cfg["results_root"], cfg["data_root"], a.labels, a.workers, a.timeout,
                                a.methods or ())
    if nrf is not None:
        nrf.to_csv(out / "nrf.csv", index=False)
    print(f"summary written to {out}")
    if errors:
        raise RuntimeError(f"{errors} comparisons ended in ERROR, see the messages above")


def cmd_paper(a):
    if a.configs and a.only and "fig5" not in a.only:
        raise ValueError("--configs only applies to fig5, which --only leaves out")
    paths = make_all(a.summary, a.networks, a.out, a.grampa_iter_fallback, a.only, a.configs)
    print(f"{len(paths)} files written to {a.out}")


def cmd_simulate(a):
    cfg = load_config(a.config)
    if a.alignments and a.gtr_params is None:
        raise ValueError("--alignments needs --gtr-params")
    for tool in ("simphy", "iqtree") if a.alignments else ("simphy",):
        if cfg["tools"].get(tool) is None:
            raise ValueError(f"tools.{tool} is not set in the config")
    simulate_network(cfg["data_root"], a.work_dir, a.sim_config, a.network, simphy_bin=cfg["tools"]["simphy"],
                     subst_rates_pkl=a.subst_rates, seed=a.seed, n_trees=a.n_trees, replicates=a.replicates,
                     batches=a.batches)
    if a.alignments:
        simulate_replicate_alignments(cfg["data_root"], a.sim_config, a.network, a.gtr_params, a.seed, a.replicates,
                                      cfg["tools"]["iqtree"])
    print(f"wrote {Path(cfg['data_root']) / 'simulated' / a.sim_config / a.network}")


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="polybench")
    sub = ap.add_subparsers(dest="command", required=True)

    def command(name, func, help, config=True):
        p = sub.add_parser(name, help=help)
        if config:
            p.add_argument("--config", required=True, help="config.yaml")
        p.set_defaults(func=func)
        return p

    def dataset_args(p, method_required):
        p.add_argument("--datasets", nargs="+", required=True,
                       help="glob(s) under data_root, e.g. 'simulated/*/*/rep*' 'empirical/*'")
        p.add_argument("--method", required=method_required, help="method name or method folder")
        p.add_argument("--params", help="YAML file overriding the method's default params")

    p = command("validate", cmd_validate, "check input datasets, for one method when --method is given")
    dataset_args(p, False)

    p = command("prep", cmd_prep, "validate and write the derived inputs to --out, for inspection")
    dataset_args(p, False)
    p.add_argument("--out", required=True, help="output folder, outside data_root")

    p = command("run", cmd_run, "run a method on datasets")
    dataset_args(p, True)
    p.add_argument("--label", help="results folder name (default: the method name)")
    p.add_argument("--executor", choices=["local", "slurm"], help="override executor.type from the config")
    p.add_argument("--dry-run", action="store_true", help="slurm: write the sbatch script without submitting")

    p = command("evaluate", cmd_evaluate, "compare results to the truth and write results_root/summary")
    p.add_argument("--labels", nargs="+",
                   help="evaluate only these labels and replace their rows in an existing summary (default: rebuild from all)")
    p.add_argument("--methods", nargs="+", help="folders of methods run by path, outside the repo's methods/")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--timeout", type=float, default=3600,
                   help="seconds per comparison (default 3600), 0 for no limit; the paper used no limit")
    p.add_argument("--nrf", nargs="+", metavar="CONFIG", help="also write summary/nrf.csv for these ILS-only configs")

    p = command("paper", cmd_paper, "make the paper figures and tables", config=False)
    p.add_argument("--summary", required=True, help="results_root/summary")
    p.add_argument("--networks", required=True, help="folder of model network .tre files")
    p.add_argument("--out", required=True)
    p.add_argument("--no-grampa-iter-fallback", dest="grampa_iter_fallback", action="store_false",
                   help="do not score failed GRAMPAIter runs with the GRAMPA result (default: as the paper, do)")
    p.add_argument("--only", nargs="+", choices=list(NAMES), help="make only these outputs")
    p.add_argument("--configs", nargs="+", help="Fig 5 columns other than the paper's 12 configs")

    p = command("simulate", cmd_simulate, "simulate gene trees for one model network with SimPhy")
    p.add_argument("--network", required=True, help="name of <data_root>/model_networks/<name>.tre")
    p.add_argument("--sim-config", required=True, choices=list(CONFIGS))
    p.add_argument("--work-dir", required=True, help="folder for the raw SimPhy output")
    p.add_argument("--subst-rates", required=True, help="substitution rate pickle")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--n-trees", type=int, default=1000)
    p.add_argument("--replicates", type=int, default=5)
    p.add_argument("--batches", type=int, default=1,
                   help="run SimPhy this many times per replicate with n-trees/batches trees each, for networks "
                        "where one SimPhy run of all trees fails")
    p.add_argument("--alignments", action="store_true", help="also simulate alignments with AliSim (tools.iqtree)")
    p.add_argument("--gtr-params", help="GTR parameter pickle for --alignments")
    return ap


def main(argv=None) -> None:
    a = parser().parse_args(argv)
    try:
        a.func(a)
    except (ValueError, FileNotFoundError, RuntimeError, subprocess.CalledProcessError) as e:
        print(f"polybench {a.command}: error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
