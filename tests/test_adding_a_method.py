# Follows docs/adding_a_method.md in a copy of the repository: a new method astral_only (returns the species
# tree) made from methods/_template, run on the toy, evaluated and shown in Fig. 5 and Table 4 next to a
# stand-in for a published method.
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ALL = ["simulated/*/*/rep*", "empirical/*"]
ILS_ONLY = ["conf_ils_low_10M", "conf_ils_medium_10M", "conf_ils_high_10M"]


def repo_copy(tmp_path) -> Path:
    repo = tmp_path / "repo"
    shutil.copytree(ROOT / "polybench", repo / "polybench", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(ROOT / "methods" / "_template", repo / "methods" / "_template")
    shutil.copytree(ROOT / "examples" / "toy", repo / "examples" / "toy", ignore=shutil.ignore_patterns("results"))
    return repo


def polybench(repo, *args):
    # the copy's polybench, so methods/ resolves to the copy's methods folder
    env = dict(os.environ, PYTHONPATH=str(repo))
    r = subprocess.run([sys.executable, "-m", "polybench.cli", *map(str, args)], cwd=repo, env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    return r


def test_adding_a_method(tmp_path):
    repo = repo_copy(tmp_path)
    method = repo / "methods" / "astral_only"
    shutil.copytree(repo / "methods" / "_template", method)
    yaml = method / "method.yaml"
    yaml.write_text(yaml.read_text().replace("name: my_method", "name: astral_only"))
    before = {p: p.read_bytes() for p in repo.rglob("*") if p.is_file() and method not in p.parents}

    cfg = "examples/toy/config.yaml"
    polybench(repo, "validate", "--config", cfg, "--datasets", *ALL, "--method", "astral_only")
    polybench(repo, "run", "--config", cfg, "--method", "astral_only", "--datasets", *ALL)
    polybench(repo, "evaluate", "--config", cfg, "--labels", "astral_only")
    summary = repo / "examples" / "toy" / "results" / "summary"
    polybench(repo, "paper", "--summary", summary, "--networks", "examples/toy/data/model_networks",
              "--out", "toy_paper", "--only", "fig5", "--configs", "toy")
    f5 = pd.read_csv(repo / "toy_paper" / "fig5_completion.csv")
    assert list(f5["label"]) == ["astral_only"] and f5.loc[0, "toy"] == 100

    # Table 4 needs the three ILS-only configs: the toy replicates under those names stand in for the Dryad data,
    # and a run of the template under the label polyphest, with its summary, for the published results
    dryad = tmp_path / "dryad"
    shutil.copytree(repo / "examples" / "toy" / "data", dryad / "data", ignore=shutil.ignore_patterns("simulated"))
    for config in ILS_ONLY:
        shutil.copytree(repo / "examples" / "toy" / "data" / "simulated" / "toy", dryad / "data" / "simulated" / config)
    (dryad / "config.yaml").write_text("data_root: data\nresults_root: results\ntools: {}\nexecutor:\n  type: local\n")
    polybench(repo, "run", "--config", dryad / "config.yaml", "--method", "methods/_template", "--label", "polyphest",
              "--datasets", *ALL)
    polybench(repo, "evaluate", "--config", dryad / "config.yaml", "--methods", "methods/_template")
    summary = dryad / "results" / "summary"
    published = {f: (summary / f).read_text().splitlines() for f in ("inventory.csv", "comparisons.csv")}

    polybench(repo, "run", "--config", dryad / "config.yaml", "--method", "astral_only", "--datasets", *ALL)
    # the stand-in label needs --methods here too: the new method is compared with it on the empirical data
    polybench(repo, "evaluate", "--config", dryad / "config.yaml", "--labels", "astral_only",
              "--methods", "methods/_template")
    for f, lines in published.items():
        assert (summary / f).read_text().splitlines()[:len(lines)] == lines  # published rows kept as they were
    pw = pd.read_csv(summary / "empirical_pairwise.csv")
    assert set(zip(pw["label_a"], pw["label_b"])) == {("astral_only", "polyphest"), ("astral_only", "published"),
                                                      ("polyphest", "published")}
    polybench(repo, "paper", "--summary", dryad / "results" / "summary", "--networks", dryad / "data" / "model_networks",
              "--out", "paper", "--only", "table4")
    t4 = pd.read_csv(repo / "paper" / "table4.csv")
    assert list(t4["label"]) == ["polyphest", "astral_only"]
    row = t4.set_index("label").loc["astral_only"]
    assert row["completed"] == row["total"] == 6 and 0 <= row["mu_distance"] <= 1

    # only the new method's folder (and the results and outputs the commands write) changed
    after = {p: p.read_bytes() for p in before if p.is_file()}
    assert after == before
