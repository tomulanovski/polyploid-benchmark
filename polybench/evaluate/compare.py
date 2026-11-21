import io
import json
import multiprocessing
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from polybench.evaluate.aggregate import aggregate
from polybench.layout import find_datasets, result_dir, true_network_path
from polybench.measures.compare_reticulations import pairwise_compare
from polybench.measures.reticulate_tree import ReticulateTree
from polybench.runner import load_method

INVENTORY_COLS = ["kind", "config", "network", "replicate", "label", "status", "runtime_s"]
COMPARISON_COLS = ["kind", "config", "network", "replicate", "label", "metric", "value", "status"]
PAIRWISE_COLS = ["dataset", "label_a", "label_b", "metric", "value"]
PUBLISHED = "published"


def load_network(path) -> ReticulateTree:
    # utf-8 because taxa may carry accents. refold=True unfolds an extended-Newick network and refolds it with
    # HOLM, so it is compared at the same resolution as a MUL-tree output; MUL-tree input is unaffected.
    return ReticulateTree(Path(path).read_text(encoding="utf-8").strip(), refold=True)


def flatten(metrics: dict) -> dict:
    # {'ploidy_diff': {'dist': .1, 'FP': 2}} -> {'ploidy_diff.dist': .1, 'ploidy_diff.FP': 2}, as legacy
    flat = {}
    for key, value in metrics.items():
        if isinstance(value, dict):
            flat.update({f"{key}.{k}": v for k, v in value.items()})
        else:
            flat[key] = value
    return flat


def compare_files(path_a, path_b, partial_match: bool) -> dict:
    trees = []
    for p in (path_a, path_b):
        try:
            trees.append(load_network(p))
        except Exception as e:
            raise RuntimeError(f"failed to load {p}: {e}") from e
    return flatten(pairwise_compare(*trees, partial_match=partial_match))


def _child(conn, fn, args):
    try:
        conn.send(("SUCCESS", fn(*args)))
    except Exception as e:
        conn.send(("ERROR", f"{type(e).__name__}: {e}"))


def isolated(fn, args, timeout_s: float) -> tuple[str, object]:
    # one process per comparison, as legacy, so a hang is killed at timeout_s (0 = no limit) and an
    # out-of-memory kill costs only this comparison
    ctx = multiprocessing.get_context("spawn")
    recv, send = ctx.Pipe(duplex=False)
    p = ctx.Process(target=_child, args=(send, fn, args))
    p.start()
    send.close()
    try:
        if recv.poll(timeout_s or None):
            return recv.recv()
        return "ERROR", f"timeout after {timeout_s}s"
    except EOFError:
        p.join()
        return "ERROR", f"comparison process died (exit code {p.exitcode})"
    finally:
        p.kill()
        p.join()
        recv.close()


def _run_all(jobs: list[tuple], workers: int, timeout_s: float) -> list[tuple[str, object]]:
    # jobs: (path_a, path_b, partial_match); threads only wait on the per-comparison processes
    with ThreadPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(lambda j: isolated(compare_files, j, timeout_s), jobs))


def _report_errors(errors: list[str], total: int, what: str) -> None:
    for e in errors:
        print(e, file=sys.stderr)
    print(f"{what}: {total - len(errors)} of {total} succeeded, {len(errors)} errors", file=sys.stderr)


def build_inventory(results_root, data_root, labels, method_dirs=()) -> pd.DataFrame:
    # method_dirs: folders of methods outside the repo's methods/, as passed to run --method <path>
    results_root = Path(results_root)
    absent = [lab for lab in labels if not (results_root / lab).is_dir()]
    if absent:
        raise FileNotFoundError(f"no result folder for labels {absent} under {results_root}")
    datasets = []
    for kind in ("simulated", "empirical"):
        if (Path(data_root) / kind).is_dir():
            datasets += find_datasets(data_root, f"{kind}/*" + ("/*/*" if kind == "simulated" else ""))
    if not datasets:
        raise FileNotFoundError(f"no simulated or empirical datasets under {data_root}")
    specs = {s.name: s for s in map(load_method, method_dirs)}
    rows = []
    for ds in datasets:
        for label in labels:
            d = result_dir(results_root, label, ds)
            row = {"kind": ds.kind, "config": ds.config, "network": ds.network, "replicate": ds.replicate,
                   "label": label, "status": "missing", "runtime_s": None, "method": None,
                   "partial_match": None, "dataset": ds.id, "data_path": str(ds.path),
                   "network_path": str(d / "network.nwk")}
            if d.is_dir():
                if not (d / "status.json").is_file():
                    raise FileNotFoundError(f"{d / 'status.json'} not found")
                st = json.loads((d / "status.json").read_text())
                if st["method"] not in specs:
                    try:
                        specs[st["method"]] = load_method(st["method"])
                    except FileNotFoundError as e:
                        raise FileNotFoundError(f"{d}: unknown method '{st['method']}', pass its folder in "
                                                f"method_dirs (evaluate --methods): {e}") from e
                row.update(status=st["status"], runtime_s=st["runtime_s"], method=st["method"],
                           partial_match=specs[st["method"]].partial_match)
                if st["status"] == "completed" and not (d / "network.nwk").is_file():
                    raise FileNotFoundError(f"{d} is completed but has no network.nwk")
            rows.append(row)
    inv = pd.DataFrame(rows)
    inv["replicate"] = inv["replicate"].astype("Int64")
    return inv


def _score(rows: pd.DataFrame, data_root, workers: int, timeout_s: float, ok: str, what: str) -> pd.DataFrame:
    # rows: inventory rows whose network_path is compared with the true network of their dataset
    jobs = []
    for r in rows.itertuples():
        truth = true_network_path(data_root, find_datasets(data_root, r.dataset)[0])
        if not truth.is_file():
            raise FileNotFoundError(f"true network {truth} not found for {r.dataset}")
        jobs.append((truth, r.network_path, bool(r.partial_match)))
    out, errors = [], []
    for r, (status, res) in zip(rows.itertuples(), _run_all(jobs, workers, timeout_s)):
        key = {"kind": r.kind, "config": r.config, "network": r.network, "replicate": r.replicate, "label": r.label}
        if status == "SUCCESS":
            out += [{**key, "metric": m, "value": v, "status": ok} for m, v in res.items()]
        else:
            errors.append(f"{r.dataset} {r.label}: {res}")
            out.append({**key, "metric": "ERROR", "value": float("nan"), "status": "ERROR"})
    _report_errors(errors, len(jobs), what)
    df = pd.DataFrame(out, columns=COMPARISON_COLS)
    df["replicate"] = df["replicate"].astype("Int64")
    return df


def compare_to_truth(inventory: pd.DataFrame, data_root, workers: int, timeout_s: float = 0) -> pd.DataFrame:
    done = inventory[(inventory["kind"] == "simulated") & (inventory["status"] == "completed")]
    return _score(done, data_root, workers, timeout_s, "SUCCESS", "comparisons to the true network")


def compare_fallback(inventory: pd.DataFrame, data_root, workers: int, timeout_s: float = 0) -> pd.DataFrame:
    # a GRAMPAIter run that failed or timed out is scored with the GRAMPA network of the same dataset, with
    # GRAMPAIter's partial_match, as the paper; rows labelled grampa_iter with status FALLBACK
    sim = inventory[inventory["kind"] == "simulated"]
    failed = sim[(sim["label"] == "grampa_iter") & sim["status"].isin(["failed", "timeout"])]
    grampa = sim.loc[(sim["label"] == "grampa") & (sim["status"] == "completed"), ["dataset", "network_path"]]
    rows = failed.drop(columns="network_path").merge(grampa, on="dataset")
    return _score(rows, data_root, workers, timeout_s, "FALLBACK", "GRAMPAIter fallback comparisons")


def compare_pairwise_empirical(inventory: pd.DataFrame, labels, workers: int = 1,
                               timeout_s: float = 0, involving=None) -> pd.DataFrame:
    # every pair of labels that completed on a dataset, plus the published network when the dataset has one;
    # involving: keep only the pairs with one of these labels
    emp = inventory[inventory["kind"] == "empirical"]
    jobs, keys = [], []
    for ds_id, g in emp.groupby("dataset", sort=True):
        done = g[(g["status"] == "completed") & g["label"].isin(labels)].set_index("label")
        entries = [(lab, done.loc[lab, "network_path"], bool(done.loc[lab, "partial_match"]))
                   for lab in labels if lab in done.index]
        published = Path(g["data_path"].iloc[0]) / "published_network.nwk"
        if published.is_file():
            entries.append((PUBLISHED, str(published), False))
        for i, (la, pa, pma) in enumerate(entries):
            for lb, pb, pmb in entries[i + 1:]:
                if involving is not None and la not in involving and lb not in involving:
                    continue
                jobs.append((pa, pb, pma or pmb))
                keys.append((g["network"].iloc[0], la, lb))
    rows, errors = [], []
    for (dataset, la, lb), (status, out) in zip(keys, _run_all(jobs, workers, timeout_s)):
        key = {"dataset": dataset, "label_a": la, "label_b": lb}
        if status == "SUCCESS":
            rows += [{**key, "metric": m, "value": v} for m, v in out.items()]
        else:
            errors.append(f"empirical/{dataset} {la} vs {lb}: {out}")
            rows.append({**key, "metric": "ERROR", "value": float("nan")})
    _report_errors(errors, len(jobs), "empirical pairwise comparisons")
    return pd.DataFrame(rows, columns=PAIRWISE_COLS)


def _merge_csv(path: Path, new: pd.DataFrame, drop) -> None:
    # rows of other labels are kept as written (read as text), drop(old) marks the rows replaced by new
    old = pd.read_csv(path, dtype=str, keep_default_na=False)
    if list(old.columns) != list(new.columns):
        raise ValueError(f"{path}: columns {list(old.columns)}, expected {list(new.columns)}")
    new = pd.read_csv(io.StringIO(new.to_csv(index=False)), dtype=str, keep_default_na=False)
    pd.concat([old[~drop(old)], new]).to_csv(path, index=False)


def _check_not_shrinking(out: Path, inv: pd.DataFrame, labels, merge: bool) -> None:
    # an unextracted Dryad package has no simulated data, rebuilding would drop every simulated row
    if not (out / "inventory.csv").is_file() or (inv["kind"] == "simulated").any():
        return
    old = pd.read_csv(out / "inventory.csv")
    if merge:
        old = old[old["label"].isin(labels)]
    if (old["kind"] == "simulated").any():
        raise ValueError(f"{out / 'inventory.csv'} has simulated rows but no simulated dataset was found; extract "
                         "the archives (data/simulated, results/<label>/simulated) first, or pass --labels to "
                         "update only some labels")


def write_summary(results_root, data_root, labels=None, workers: int = 1, timeout_s: float = 0,
                  method_dirs=()) -> tuple[Path, int]:
    # without labels the summary is rebuilt from every result folder; with labels and an existing summary only
    # those labels are scored and their rows replace the old ones, other labels' rows are kept as they are.
    # Returns the summary folder and the number of comparisons that ended in ERROR.
    results_root = Path(results_root)
    out = results_root / "summary"
    on_disk = sorted(p.name for p in results_root.iterdir() if p.is_dir() and p.name != "summary")
    merge = labels is not None and (out / "inventory.csv").is_file()
    if labels is None:
        labels = on_disk
    inv = build_inventory(results_root, data_root, labels, method_dirs)
    _check_not_shrinking(out, inv, labels, merge)
    cmp = compare_to_truth(inv, data_root, workers, timeout_s)
    out.mkdir(exist_ok=True)
    if not merge:
        fb = compare_fallback(inv, data_root, workers, timeout_s)
        pw = compare_pairwise_empirical(inv, labels, workers, timeout_s)
        inv[INVENTORY_COLS].to_csv(out / "inventory.csv", index=False)
        cmp.to_csv(out / "comparisons.csv", index=False)
        fb.to_csv(out / "fallback.csv", index=False)
        pw.to_csv(out / "empirical_pairwise.csv", index=False)
        aggregate(cmp).to_csv(out / "aggregated.csv", index=False)
        return out, _errors(cmp, fb, pw)
    others = [lab for lab in on_disk if lab not in labels]
    both = pd.concat([inv, build_inventory(results_root, data_root, others, method_dirs)])
    pw = compare_pairwise_empirical(both, sorted(set(labels) | set(others)), workers, timeout_s, involving=labels)
    fb = pd.DataFrame(columns=COMPARISON_COLS)
    if {"grampa", "grampa_iter"} & set(labels) or not (out / "fallback.csv").is_file():
        fb = compare_fallback(both, data_root, workers, timeout_s)
        fb.to_csv(out / "fallback.csv", index=False)
    _merge_csv(out / "inventory.csv", inv[INVENTORY_COLS], lambda d: d["label"].isin(labels))
    _merge_csv(out / "comparisons.csv", cmp, lambda d: d["label"].isin(labels))
    _merge_csv(out / "empirical_pairwise.csv", pw, lambda d: d["label_a"].isin(labels) | d["label_b"].isin(labels))
    aggregate(pd.read_csv(out / "comparisons.csv")).to_csv(out / "aggregated.csv", index=False)
    return out, _errors(cmp, fb, pw)


def _errors(cmp: pd.DataFrame, fb: pd.DataFrame, pw: pd.DataFrame) -> int:
    return int((cmp["status"] == "ERROR").sum() + (fb["status"] == "ERROR").sum() + (pw["metric"] == "ERROR").sum())
