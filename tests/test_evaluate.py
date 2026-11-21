import json
import os
import shutil
import time

import pandas as pd
import pytest

from polybench.evaluate import compare
from polybench.evaluate.aggregate import aggregate
from polybench.evaluate.compare import (build_inventory, compare_pairwise_empirical, compare_to_truth,
                                        write_summary)

# Table 3 of the manuscript: truth A, inferred B and D
A = "(((q,p),(r,p)),(s,(t,u)));"
B = "(((q,p),(s,p)),(r,(t,u)));"
D = "(((q,p),(r,p)),((s,u),(t,u)));"
A_ENEWICK = "(((q,(p)#H1),(r,#H1)),(s,(t,u)));"
SIM = "simulated/cfg/A/rep1"
EMP = "empirical/emp"


def put_result(results, label, ds_id, method, status="completed", newick=None, runtime=12.5):
    d = results / label / ds_id
    d.mkdir(parents=True)
    if newick is not None:
        (d / "network.nwk").write_text(newick + "\n")
    (d / "status.json").write_text(json.dumps({"status": status, "runtime_s": runtime, "method": method,
                                               "version": "x", "params": {}}))


@pytest.fixture
def layout(tmp_path):
    data, results = tmp_path / "data", tmp_path / "results"
    (data / "model_networks").mkdir(parents=True)
    (data / "model_networks" / "A.tre").write_text(A + "\n")
    for ds in (SIM, EMP):
        (data / ds).mkdir(parents=True)
        (data / ds / "gene_trees.nwk").write_text(A + "\n")
    put_result(results, "polyphest_p50", SIM, "polyphest", newick=B)
    put_result(results, "padre", SIM, "padre", newick=D)
    put_result(results, "grampa", SIM, "grampa", status="failed", runtime=3.0)
    put_result(results, "mpsugar", SIM, "mpsugar", newick=A_ENEWICK)
    put_result(results, "polyphest_p50", EMP, "polyphest", newick=B)
    put_result(results, "padre", EMP, "padre", newick=D)
    return data, results


def metric(df, label, name):
    rows = df[(df["label"] == label) & (df["metric"] == name)]
    assert len(rows) == 1, rows
    return rows["value"].iloc[0]


def test_inventory_statuses(layout):
    data, results = layout
    inv = build_inventory(results, data, ["polyphest_p50", "padre", "grampa"])
    sim = inv[inv["kind"] == "simulated"].set_index("label")
    assert len(sim) == 3
    assert sim.loc["polyphest_p50", "status"] == "completed"
    assert sim.loc["padre", "status"] == "completed"
    assert sim.loc["grampa", "status"] == "failed"
    assert sim.loc["grampa", "runtime_s"] == 3.0
    assert sim.loc["polyphest_p50", "method"] == "polyphest"
    emp = inv[inv["kind"] == "empirical"].set_index("label")
    assert emp.loc["grampa", "status"] == "missing"
    assert pd.isna(emp.loc["grampa", "runtime_s"])
    assert pd.isna(emp.loc["padre", "replicate"]) and pd.isna(emp.loc["padre", "config"])


def test_inventory_fails_on_unknown_method(layout):
    data, results = layout
    put_result(results, "mystery", SIM, "no_such_method", newick=B)
    with pytest.raises(FileNotFoundError, match="no_such_method"):
        build_inventory(results, data, ["mystery"])


def test_inventory_fails_on_unknown_label(layout):
    data, results = layout
    with pytest.raises(FileNotFoundError, match="typo"):
        build_inventory(results, data, ["typo"])


def test_inventory_fails_on_folder_without_status(layout):
    data, results = layout
    (results / "padre" / EMP / "status.json").unlink()
    with pytest.raises(FileNotFoundError, match="status.json"):
        build_inventory(results, data, ["padre"])


def test_compare_to_truth_table3(layout):
    data, results = layout
    inv = build_inventory(results, data, ["polyphest_p50", "padre", "grampa"])
    cmp = compare_to_truth(inv, data, workers=2)
    assert set(cmp["status"]) == {"SUCCESS"}
    assert set(cmp["label"]) == {"polyphest_p50", "padre"}  # failed runs and empirical data are not compared
    assert metric(cmp, "polyphest_p50", "ret_leaf_jaccard.dist") == 0.0
    assert metric(cmp, "padre", "ret_leaf_jaccard.dist") == 0.5
    assert metric(cmp, "padre", "num_rets_bias") == 1
    row = cmp.iloc[0]
    assert (row["kind"], row["config"], row["network"], row["replicate"]) == ("simulated", "cfg", "A", 1)
    assert list(cmp.columns) == ["kind", "config", "network", "replicate", "label", "metric", "value", "status"]


def test_enewick_result_is_refolded(layout):
    data, results = layout
    cmp = compare_to_truth(build_inventory(results, data, ["mpsugar"]), data, workers=1)
    assert metric(cmp, "mpsugar", "mu_distance") == 0.0
    assert metric(cmp, "mpsugar", "ret_leaf_jaccard.dist") == 0.0


def test_partial_match_comes_from_the_method(layout):
    data, results = layout
    put_result(results, "grampa_x", SIM, "grampa", newick=D)
    inv = build_inventory(results, data, ["grampa_x"])
    assert inv[inv["kind"] == "simulated"]["partial_match"].tolist() == [True]


def test_comparison_error_is_recorded(layout):
    data, results = layout
    (results / "padre" / SIM / "network.nwk").write_text("this is not newick\n")
    cmp = compare_to_truth(build_inventory(results, data, ["padre"]), data, workers=1)
    assert cmp["status"].tolist() == ["ERROR"]
    assert cmp["metric"].tolist() == ["ERROR"]


def test_missing_truth_fails_loud(layout):
    data, results = layout
    (data / "model_networks" / "A.tre").unlink()
    with pytest.raises(FileNotFoundError, match="A.tre"):
        compare_to_truth(build_inventory(results, data, ["padre"]), data, workers=1)


def test_empirical_pairwise_never_needs_truth(layout):
    data, results = layout
    (data / "model_networks" / "A.tre").unlink()
    (data / EMP / "published_network.nwk").write_text(A + "\n")
    labels = ["polyphest_p50", "padre", "grampa"]
    pw = compare_pairwise_empirical(build_inventory(results, data, labels), labels)
    assert list(pw.columns) == ["dataset", "label_a", "label_b", "metric", "value"]
    assert set(pw["dataset"]) == {"emp"}
    pairs = set(zip(pw["label_a"], pw["label_b"]))
    assert pairs == {("polyphest_p50", "padre"), ("polyphest_p50", "published"), ("padre", "published")}
    vs_pub = pw[(pw["label_a"] == "padre") & (pw["label_b"] == "published")]
    assert vs_pub.set_index("metric").loc["ret_leaf_jaccard.dist", "value"] == 0.5


def test_write_summary(layout):
    data, results = layout
    out, errors = write_summary(results, data, workers=1)
    assert errors == 0
    inv = pd.read_csv(out / "inventory.csv")
    assert list(inv.columns) == ["kind", "config", "network", "replicate", "label", "status", "runtime_s"]
    assert len(inv) == 8  # 2 datasets x 4 labels
    assert inv["replicate"].dropna().tolist() == [1, 1, 1, 1]
    cmp = pd.read_csv(out / "comparisons.csv")
    assert set(cmp["label"]) == {"polyphest_p50", "padre", "mpsugar"}
    assert set(pd.read_csv(out / "empirical_pairwise.csv")["dataset"]) == {"emp"}
    agg = pd.read_csv(out / "aggregated.csv")
    assert list(agg.columns) == ["network", "config", "label", "metric", "mean", "std", "min", "max", "n_valid"]
    assert list(pd.read_csv(out / "fallback.csv").columns) == list(cmp.columns)


def test_write_summary_counts_errors(layout):
    data, results = layout
    (results / "padre" / SIM / "network.nwk").write_text("this is not newick\n")
    _, errors = write_summary(results, data, workers=1)
    assert errors == 1


def test_grampa_iter_fallback_uses_its_own_partial_match(layout):
    # D scores 0 against A with grampa's partial_match and 0.5 without it, as grampa_iter is scored
    data, results = layout
    (data / "simulated/cfg/A/rep2").mkdir()
    (data / "simulated/cfg/A/rep2/gene_trees.nwk").write_text(A + "\n")
    shutil.rmtree(results / "grampa")
    put_result(results, "grampa", SIM, "grampa", newick=D)
    put_result(results, "grampa", "simulated/cfg/A/rep2", "grampa", newick=D)  # grampa_iter missing: no fallback
    put_result(results, "grampa_iter", SIM, "grampa_iter", status="timeout")
    out, _ = write_summary(results, data, workers=1)
    cmp, fb = pd.read_csv(out / "comparisons.csv"), pd.read_csv(out / "fallback.csv")
    grampa = cmp[(cmp["label"] == "grampa") & (cmp["replicate"] == 1)]
    assert metric(grampa, "grampa", "ret_leaf_jaccard.dist") == 0.0
    assert metric(fb, "grampa_iter", "ret_leaf_jaccard.dist") == 0.5
    assert set(fb["status"]) == {"FALLBACK"} and set(fb["replicate"]) == {1}
    # the --labels merge path recomputes it
    (results / "grampa" / SIM / "network.nwk").write_text(A + "\n")
    write_summary(results, data, labels=["grampa"], workers=1)
    fb = pd.read_csv(out / "fallback.csv")
    assert metric(fb, "grampa_iter", "ret_leaf_jaccard.dist") == 0.0


def test_rebuild_refuses_to_drop_simulated_rows(layout):
    data, results = layout
    write_summary(results, data, workers=1)
    shutil.rmtree(data / "simulated")
    with pytest.raises(ValueError, match="extract"):
        write_summary(results, data, workers=1)


def test_aggregate_over_replicates():
    rows = [("A", "cfg", "x", r, "mu_distance", v, "SUCCESS") for r, v in ((1, 0.2), (2, 0.4))]
    rows += [("A", "cfg", "x", 3, "ERROR", float("nan"), "ERROR"),
             ("B", "cfg", "x", 1, "mu_distance", 0.1, "SUCCESS")]
    df = pd.DataFrame(rows, columns=["network", "config", "label", "replicate", "metric", "value", "status"])
    agg = aggregate(df).set_index("network")
    assert agg.loc["A", "mean"] == pytest.approx(0.3)
    assert agg.loc["A", "n_valid"] == 2
    assert agg.loc["B", "std"] == 0  # one replicate
    assert set(agg["metric"]) == {"mu_distance"}


def test_isolated_timeout_and_crash():
    t0 = time.monotonic()
    status, msg = compare.isolated(time.sleep, (30,), timeout_s=2)
    assert status == "ERROR" and "timeout" in msg
    assert time.monotonic() - t0 < 20
    status, msg = compare.isolated(os._exit, (3,), timeout_s=0)
    assert status == "ERROR" and "exit code 3" in msg


def test_write_summary_merges_given_labels(layout):
    data, results = layout
    (data / EMP / "published_network.nwk").write_text(A + "\n")
    out, _ = write_summary(results, data, labels=["polyphest_p50", "padre", "grampa"], workers=1)
    before = {f: (out / f).read_text().splitlines() for f in ("inventory.csv", "comparisons.csv",
                                                             "empirical_pairwise.csv")}
    # padre is rerun with another network, mpsugar is new on the empirical dataset
    (results / "padre" / SIM / "network.nwk").write_text(A + "\n")
    put_result(results, "mpsugar", EMP, "mpsugar", newick=A_ENEWICK)
    write_summary(results, data, labels=["padre", "mpsugar"], workers=1)

    for f, old in before.items():
        new = (out / f).read_text().splitlines()
        col = 4 if f != "empirical_pairwise.csv" else None
        def touched(line):
            cells = line.split(",")
            return (cells[col] in ("padre", "mpsugar")) if col is not None else bool({"padre", "mpsugar"} & set(cells[1:3]))
        kept = [l for l in old[1:] if not touched(l)]
        assert new[0] == old[0] and new[1:1 + len(kept)] == kept, f  # other labels' rows byte-identical
    cmp = pd.read_csv(out / "comparisons.csv")
    assert set(cmp["label"]) == {"polyphest_p50", "padre", "mpsugar"}
    assert metric(cmp, "padre", "mu_distance") == 0  # replaced by the rerun
    inv = pd.read_csv(out / "inventory.csv")
    assert len(inv) == 8 and (inv.groupby("label").size() == 2).all()
    pw = pd.read_csv(out / "empirical_pairwise.csv")
    pairs = set(zip(pw["label_a"], pw["label_b"]))
    assert pairs == {("padre", "polyphest_p50"), ("padre", "published"), ("polyphest_p50", "published"),
                     ("mpsugar", "padre"), ("mpsugar", "polyphest_p50"), ("mpsugar", "published")}
    assert len(pw) == len(pw.drop_duplicates(["dataset", "label_a", "label_b", "metric"]))
    agg = pd.read_csv(out / "aggregated.csv")
    assert set(agg["label"]) == {"polyphest_p50", "padre", "mpsugar"}
