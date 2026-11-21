import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
TOY = ROOT / "examples" / "toy"
METHODS = ROOT / "tests" / "fixtures" / "methods"
ALL = ["simulated/*/*/rep*", "empirical/*"]
SIMULATED = {("Shahrestani_2015", 1), ("Wisecaver_2023", 1)}


def cli(cwd, *args):
    return subprocess.run([sys.executable, "-m", "polybench.cli", *map(str, args)], cwd=cwd,
                          capture_output=True, text=True)


def ok(r):
    assert r.returncode == 0, r.stdout + r.stderr
    return r


@pytest.fixture
def toy(tmp_path):
    # commands run from tmp_path/work, not the toy folder: config paths resolve against the config file
    shutil.copytree(TOY, tmp_path / "toy", ignore=shutil.ignore_patterns("results"))
    (tmp_path / "work").mkdir()
    return tmp_path / "toy"


def polybench(toy, *args):
    return cli(toy.parent / "work", *args)


def test_toy_pipeline(toy):
    cfg, work = toy / "config.yaml", toy.parent / "work"
    ok(polybench(toy, "validate", "--config", cfg, "--datasets", *ALL))
    ok(polybench(toy, "prep", "--config", cfg, "--datasets", *ALL, "--method", METHODS / "echo_astral", "--out", "prep"))
    assert (work / "prep" / "empirical" / "toy_empirical" / "dataset" / "species_tree.nwk").is_file()
    assert (work / "prep" / "simulated" / "toy" / "Wisecaver_2023" / "rep1" / "dataset" / "gene_trees.nwk").is_file()
    ok(polybench(toy, "run", "--config", cfg, "--method", METHODS / "echo_astral", "--datasets", *ALL))
    r = polybench(toy, "evaluate", "--config", cfg, "--methods", METHODS / "echo_astral", "--nrf", "toy")
    assert r.returncode == 1 and "unknown config toy" in r.stderr  # checked before any comparison is written
    summary = toy / "results" / "summary"
    assert not summary.exists()
    ok(polybench(toy, "evaluate", "--config", cfg, "--methods", METHODS / "echo_astral"))
    cmp = pd.read_csv(summary / "comparisons.csv")
    assert set(zip(cmp["network"], cmp["replicate"])) == SIMULATED
    assert (cmp["status"] == "SUCCESS").all() and "mu_distance" in set(cmp["metric"])
    pw = pd.read_csv(summary / "empirical_pairwise.csv")
    assert set(pw["dataset"]) == {"toy_empirical"} and set(pw["label_b"]) == {"published"}
    inv = pd.read_csv(summary / "inventory.csv")
    assert len(inv) == 3 and (inv["status"] == "completed").all()
    ok(polybench(toy, "paper", "--summary", summary, "--networks", toy / "data" / "model_networks", "--out", "paper",
                 "--only", "fig5", "table_s2", "--configs", "toy"))
    assert (work / "paper" / "fig5.pdf").stat().st_size > 0
    f5 = pd.read_csv(work / "paper" / "fig5_completion.csv")
    assert list(f5["label"]) == ["echo_astral"] and f5.loc[0, "toy"] == 100
    assert len(pd.read_csv(work / "paper" / "table_s2.csv")) == 2

    (toy / "results" / "echo_astral" / "simulated" / "toy" / "Wisecaver_2023" / "rep1" / "network.nwk").write_text("(")
    r = polybench(toy, "evaluate", "--config", cfg, "--methods", METHODS / "echo_astral")
    assert r.returncode == 1 and "1 comparisons ended in ERROR" in r.stderr


def test_slurm_script_uses_absolute_paths(toy):
    ok(polybench(toy, "run", "--config", "../toy/config.yaml", "--method", METHODS / "echo_astral", "--datasets", *ALL,
                 "--executor", "slurm", "--dry-run"))
    label_dir = (toy / "results" / "echo_astral").resolve().as_posix()
    script = (toy / "results" / "echo_astral" / "submit.sh").read_text()
    assert f"--output={label_dir}/logs/" in script
    assert f"--dataset-list {label_dir}/dataset_list.txt" in script
    assert f"--config {(toy / 'config.yaml').resolve().as_posix()}" in script


def test_run_fails_when_a_dataset_fails(toy):
    r = polybench(toy, "run", "--config", toy / "config.yaml", "--method", METHODS / "crash", "--datasets", *ALL)
    assert r.returncode != 0
    assert "simulated/toy/Shahrestani_2015/rep1: failed" in r.stderr
    assert "3 of 3 datasets did not complete" in r.stderr


def test_validate_lists_problems(toy):
    (toy / "data" / "empirical" / "toy_empirical" / "gene_trees.nwk").write_text("((a,b);\n")
    (toy / "data" / "simulated" / "toy" / "Wisecaver_2023" / "rep1" / "species_tree.nwk").unlink()
    r = polybench(toy, "validate", "--config", toy / "config.yaml", "--datasets", *ALL, "--method", "grampa")
    assert r.returncode != 0
    assert "empirical/toy_empirical:" in r.stderr and "not valid Newick" in r.stderr
    assert "simulated/toy/Shahrestani_2015" not in r.stderr
    # grampa needs a species tree, the runner infers it with ASTRAL, which the toy config does not set
    r = polybench(toy, "run", "--config", toy / "config.yaml", "--method", "grampa", "--datasets",
                  "simulated/toy/Wisecaver_2023/rep1")
    assert r.returncode != 0 and "tools.astral is not set" in r.stderr


@pytest.mark.needs_astral
def test_toy_run_infers_species_tree_with_astral(toy):
    astral = os.environ.get("POLYBENCH_ASTRAL")
    if not astral:
        pytest.skip("tools.astral is unset: set POLYBENCH_ASTRAL to an ASTRAL binary to run this test")
    cfg = toy / "config.yaml"
    cfg.write_text(cfg.read_text().replace("tools: {}", f"tools: {{astral: '{astral}'}}"))
    ds = Path("simulated") / "toy" / "Shahrestani_2015" / "rep1"
    (toy / "data" / ds / "species_tree.nwk").unlink()
    ok(polybench(toy, "run", "--config", cfg, "--method", METHODS / "echo_astral", "--datasets", ds.as_posix()))
    assert (toy / "results" / "echo_astral" / ds / "network.nwk").stat().st_size > 0
