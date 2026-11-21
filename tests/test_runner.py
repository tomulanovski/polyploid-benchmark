import json
import re
import shutil
import warnings
from pathlib import Path

import pytest

from polybench import runner
from polybench.layout import Dataset
from polybench.runner import load_method, parse_slurm_time, run_one, submit

FIX = Path(__file__).parent / "fixtures"
METHODS = FIX / "methods"
GOOD = FIX / "bad_datasets" / "good"


@pytest.fixture
def dataset(tmp_path):
    root = tmp_path / "data"
    d = root / "empirical" / "good"
    shutil.copytree(GOOD, d)
    return Dataset("empirical", None, "good", None, d, root)


@pytest.fixture
def cfg(tmp_path):
    return {"data_root": str(tmp_path / "data"), "results_root": str(tmp_path / "res"),
            "tools": {"astral": "astral_x"},
            "executor": {"type": "local", "slurm_header": ["#SBATCH --partition=P", "#SBATCH --account=A"],
                         "conda_activate": ""}}


def status_of(out):
    return json.loads((out / "status.json").read_text())


def go(name, dataset, cfg, tmp_path, timeout_s=60):
    out = tmp_path / "out"
    run_one(load_method(METHODS / name), dataset, out, {}, "lbl", timeout_s, cfg)
    return out


def test_completed(dataset, cfg, tmp_path):
    out = go("echo_astral", dataset, cfg, tmp_path)
    st = status_of(out)
    assert st["status"] == "completed" and st["runtime_s"] > 0
    assert st["method"] == "echo_astral" and st["version"] == "0.1" and st["params"] == {"k": 1}
    assert (out / "network.nwk").read_text() == (GOOD / "species_tree.nwk").read_text()


def test_params_override(dataset, cfg, tmp_path):
    out = tmp_path / "out"
    run_one(load_method(METHODS / "echo_astral"), dataset, out, {"k": 5}, "lbl", 60, cfg)
    assert status_of(out)["params"] == {"k": 5}
    assert json.loads((out / "params.json").read_text()) == {"k": 5}


def test_crash_is_failed(dataset, cfg, tmp_path):
    out = go("crash", dataset, cfg, tmp_path)
    assert status_of(out)["status"] == "failed"
    assert "boom-on-stderr" in (out / "log.txt").read_text()


def test_no_output_is_failed(dataset, cfg, tmp_path):
    out = go("noout", dataset, cfg, tmp_path)
    assert status_of(out)["status"] == "failed"


def test_timeout(dataset, cfg, tmp_path):
    out = go("sleep", dataset, cfg, tmp_path, timeout_s=1)
    st = status_of(out)
    assert st["status"] == "timeout" and st["runtime_s"] < 15


def test_popen_error_is_failed(dataset, cfg, tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "run_command", lambda spec, cfg, args: [str(tmp_path / "no_such_binary")])
    out = go("echo_astral", dataset, cfg, tmp_path)
    assert status_of(out)["status"] == "failed"
    assert "Traceback" in (out / "log.txt").read_text()


def test_method_sees_only_its_inputs(dataset, cfg, tmp_path):
    for f in ("published_network.nwk", "simphy_species_tree.nwk", "simphy_params.txt"):
        (dataset.path / f).write_text("x")
    work = tmp_path / "work"
    work_ds = runner.prepare_inputs(load_method(METHODS / "echo_astral"), dataset, cfg, work)
    assert sorted(p.name for p in work_ds.iterdir()) == ["gene_trees.nwk", "species_tree.nwk"]
    spec = runner.MethodSpec(name="x", version="", env=None, requires=["alignments", "ploidy"], output="network",
                             partial_match=False, resources={"threads": 1}, params={}, path=Path("."))
    work_ds = runner.prepare_inputs(spec, dataset, cfg, work)
    assert sorted(p.name for p in work_ds.iterdir()) == ["alignments", "gene_trees.nwk", "ploidy.tsv"]


def test_ploidy_derived_into_work_dir(tmp_path, cfg):
    root = tmp_path / "data"
    d = root / "empirical" / "good"
    shutil.copytree(GOOD, d)
    (d / "ploidy.tsv").unlink()
    ds = Dataset("empirical", None, "good", None, d, root)
    spec = load_method(METHODS / "echo_astral")
    spec.requires.append("ploidy")
    out = tmp_path / "out"
    run_one(spec, ds, out, {}, "lbl", 60, cfg)
    assert not (d / "ploidy.tsv").exists()
    assert (out / "work" / "dataset" / "ploidy.tsv").is_file()
    assert (out / "work" / "dataset" / "gene_trees.nwk").is_file()


def test_missing_species_tree_runs_astral(tmp_path, cfg, monkeypatch):
    root = tmp_path / "data"
    d = root / "empirical" / "good"
    shutil.copytree(GOOD, d)
    (d / "species_tree.nwk").unlink()
    ds = Dataset("empirical", None, "good", None, d, root)
    calls = []

    def fake_astral(gt, out, astral_bin, threads):
        calls.append(astral_bin)
        Path(out).write_text("(a,b);")

    monkeypatch.setattr(runner, "run_astral", fake_astral)
    out = tmp_path / "out"
    run_one(load_method(METHODS / "echo_astral"), ds, out, {}, "lbl", 60, cfg)
    assert calls == ["astral_x"] and status_of(out)["status"] == "completed"
    assert not (d / "species_tree.nwk").exists()


def test_method_prepare_hook(dataset, cfg, tmp_path):
    m = tmp_path / "m"
    shutil.copytree(METHODS / "echo_astral", m)
    (m / "prepare.py").write_text(
        "def prepare(dataset_dir, work_dir, params):\n    (dataset_dir / 'species_tree.nwk').write_text('(x,y);')\n")
    out = tmp_path / "out"
    run_one(load_method(m), dataset, out, {}, "lbl", 60, cfg)
    assert (out / "network.nwk").read_text() == "(x,y);"


def test_tool_env_vars(dataset, cfg, tmp_path):
    m = tmp_path / "m"
    shutil.copytree(METHODS / "echo_astral", m)
    (m / "run.sh").write_text('echo "$POLYBENCH_ASTRAL" > "$2/network.nwk"\n')
    out = tmp_path / "out"
    run_one(load_method(m), dataset, out, {}, "lbl", 60, cfg)
    assert (out / "network.nwk").read_text().strip() == "astral_x"


def test_conda_command():
    spec = load_method(METHODS / "echo_astral")
    script = (spec.path / "run.sh").as_posix()
    spec.env = "myenv"
    cfg = {"executor": {"conda_activate": "source c.sh && conda activate"}}
    cmd = runner.run_command(spec, cfg, ["a b", "c"])
    assert cmd[1:3] == ["-c", f"source c.sh && conda activate myenv && bash {script} 'a b' c"]
    spec.env = None
    assert runner.run_command(spec, cfg, ["a"])[1:] == [script, "a"]


def test_parse_slurm_time():
    assert parse_slurm_time("5-00:00:00") == 432000
    assert parse_slurm_time("02:30:00") == 9000


def test_load_method_errors(tmp_path):
    (tmp_path / "method.yaml").write_text("name: x\n")
    with pytest.raises(ValueError, match="version"):
        load_method(tmp_path)
    with pytest.raises(FileNotFoundError):
        load_method("no_such_method")


def datasets3(tmp_path):
    root = tmp_path / "data"
    out = []
    for k in (1, 2, 3):
        d = root / "simulated" / "c" / "Net" / f"rep{k}"
        shutil.copytree(GOOD, d)
        out.append(Dataset("simulated", "c", "Net", k, d, root))
    (root / "model_networks").mkdir()
    (root / "model_networks" / "Net.tre").write_text("((a,b),(c,(c,d)));")
    return out


def test_slurm_script_render(tmp_path, cfg):
    cfg["executor"]["type"] = "slurm"
    spec = load_method(METHODS / "echo_astral")
    script = submit(spec, datasets3(tmp_path), "lbl", {"k": 2}, cfg, tmp_path / "config.yaml", dry_run=True)
    text = script.read_text()
    assert "#SBATCH --partition=P" in text and "#SBATCH --account=A" in text
    assert text.count("--time") == 1 and "#SBATCH --time=5-00:00:00" in text and "#SBATCH --array=0-2" in text
    assert "#SBATCH --mem=1G" in text and "#SBATCH --cpus-per-task=1" in text
    assert "--index $SLURM_ARRAY_TASK_ID" in text and "-m polybench.runner --one" in text
    allowed = (tmp_path.as_posix(), METHODS.as_posix())
    for path in re.findall(r"(?<![\w$-])(?:[A-Za-z]:)?/[\w./:-]+", text.replace("\\", "/")):
        assert path.startswith(allowed) or path == "/usr/bin/env", path
    lst = (tmp_path / "res" / "lbl" / "dataset_list.txt").read_text().split()
    assert lst == ["simulated/c/Net/rep1", "simulated/c/Net/rep2", "simulated/c/Net/rep3"]
    assert json.loads((tmp_path / "res" / "lbl" / "params.json").read_text()) == {"k": 2}


@pytest.mark.filterwarnings("ignore:.*fewer than 10:UserWarning")  # small fixtures on purpose
def test_submit_validation_failure_runs_nothing(tmp_path, cfg):
    ds = datasets3(tmp_path)
    (ds[1].path / "gene_trees.nwk").write_text("(a,b,c);\n")
    spec = load_method(METHODS / "echo_astral")
    with pytest.raises(ValueError, match="rep2"):
        submit(spec, ds, "lbl", {}, cfg, tmp_path / "config.yaml")
    assert not (tmp_path / "res").exists()


def test_submit_local_loop_and_one(tmp_path, cfg):
    ds = datasets3(tmp_path)
    spec = load_method(METHODS / "echo_astral")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = submit(spec, ds, "lbl", {}, cfg, tmp_path / "config.yaml")
    assert set(results) == {d.id for d in ds}
    for d in ds:
        assert status_of(tmp_path / "res" / "lbl" / d.id)["status"] == "completed"
    # the --one entry point resolves the dataset from the list and runs it
    cfgfile = tmp_path / "config.yaml"
    cfgfile.write_text(f"data_root: {cfg['data_root']}\nresults_root: {tmp_path / 'res2'}\ntools: {{astral: a}}\nexecutor: {{type: local}}\n")
    lst = tmp_path / "list.txt"
    lst.write_text("\n".join(d.id for d in ds) + "\n")
    pj = tmp_path / "p.json"
    pj.write_text("{}")
    runner.main(["--one", "--method", str(METHODS / "echo_astral"), "--label", "x", "--params", str(pj),
                 "--dataset-list", str(lst), "--index", "1", "--config", str(cfgfile)])
    assert status_of(tmp_path / "res2" / "x" / ds[1].id)["status"] == "completed"


def test_path_loaded_method_validate_hook(tmp_path, cfg):
    m = tmp_path / "mymethod"
    shutil.copytree(METHODS / "echo_astral", m)
    (m / "validate.py").write_text("def extra_checks(dataset_dir, params):\n    return ['hook ran']\n")
    with pytest.raises(ValueError, match="hook ran"):
        submit(load_method(m), datasets3(tmp_path), "lbl", {}, cfg, tmp_path / "config.yaml")


def test_prep_exception_is_failed_and_loop_continues(tmp_path, cfg):
    ds = datasets3(tmp_path)
    m = tmp_path / "m"
    shutil.copytree(METHODS / "echo_astral", m)
    (m / "prepare.py").write_text(
        "def prepare(dataset_dir, work_dir, params):\n"
        "    if 'rep1' in str(work_dir):\n        raise RuntimeError('prep-blew-up')\n")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = submit(load_method(m), ds, "lbl", {}, cfg, tmp_path / "config.yaml")
    out1 = tmp_path / "res" / "lbl" / ds[0].id
    assert status_of(out1)["status"] == "failed" and "prep-blew-up" in (out1 / "log.txt").read_text()
    assert results[ds[1].id]["status"] == "completed" and results[ds[2].id]["status"] == "completed"
    with pytest.raises(RuntimeError, match="1 of 3"):
        runner.raise_if_not_completed(results)


def test_raise_if_not_completed_ok():
    runner.raise_if_not_completed({"a": {"status": "completed"}})


def test_stale_status_removed(dataset, cfg, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "status.json").write_text('{"status": "completed"}')
    spec = load_method(METHODS / "echo_astral")
    spec.requires = ["species_tree"]
    cfg["tools"]["astral"] = None
    (dataset.path / "species_tree.nwk").unlink()
    run_one(spec, dataset, out, {}, "lbl", 60, cfg)
    assert status_of(out)["status"] == "failed" and "tools.astral" in (out / "log.txt").read_text()


def test_time_budget():
    from polybench.executors import time_budget
    assert time_budget(432000, 100) == 432000 - 100 - 600
    assert time_budget(3600, 0) == 3420
    assert time_budget(1, 5) == 1.0


def test_quoted_paths_with_space(tmp_path, cfg):
    m = tmp_path / "my method"
    shutil.copytree(METHODS / "echo_astral", m)
    spec = load_method(m)
    spec.env = "e"
    cmd = runner.run_command(spec, {"executor": {"conda_activate": "act"}}, ["x y"])
    assert f"bash '{m.as_posix()}/run.sh' 'x y'" in cmd[-1]
    cfg["executor"]["type"] = "slurm"
    text = submit(spec, datasets3(tmp_path), "lbl", {}, cfg, tmp_path / "my config.yaml", dry_run=True).read_text()
    assert f"--method '{m.as_posix()}'" in text and f"--config '{(tmp_path / 'my config.yaml').as_posix()}'" in text


def test_slurm_rejects_whitespace_results_root(tmp_path, cfg):
    cfg["executor"]["type"] = "slurm"
    cfg["results_root"] = str(tmp_path / "my res")
    with pytest.raises(ValueError, match="whitespace"):
        submit(load_method(METHODS / "echo_astral"), datasets3(tmp_path), "lbl", {}, cfg, tmp_path / "c.yaml", dry_run=True)
    assert not (tmp_path / "my res").exists()


def test_root_mismatch_rejected(tmp_path, cfg):
    cfg["data_root"] = str(tmp_path / "elsewhere")
    with pytest.raises(ValueError, match="data_root"):
        submit(load_method(METHODS / "echo_astral"), datasets3(tmp_path), "lbl", {}, cfg, tmp_path / "c.yaml")


def test_rerun_starts_with_clean_work_dir(dataset, cfg, tmp_path):
    out = go("sentinel", dataset, cfg, tmp_path)
    assert status_of(out)["status"] == "completed"
    assert (out / "work" / "sentinel").is_file()
    out = go("sentinel", dataset, cfg, tmp_path)
    assert status_of(out)["status"] == "completed", (out / "log.txt").read_text()


def test_bool_flags_key(tmp_path):
    assert load_method(METHODS / "echo_astral").bool_flags == "value"
    assert load_method("grampa_iter").bool_flags == "switch"
    assert load_method("polyphest").bool_flags == "value"
    d = tmp_path / "m"
    shutil.copytree(METHODS / "echo_astral", d)
    with open(d / "method.yaml", "a") as f:
        f.write("bool_flags: maybe\n")
    with pytest.raises(ValueError, match="bool_flags"):
        load_method(d)
