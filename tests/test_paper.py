import io
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from polybench.paper import make_all
from polybench.paper import common, fig5, fig6_fig7, stats, table4
from polybench.paper.fig3_table3 import table3
from polybench.paper.table_s2 import network_stats

MODEL_NETWORKS = Path(__file__).parent / "data" / "model_networks"
NETWORKS = ["Ding_2023", "Wu_2015", "Morales-Briones_2021"]  # Ding, Wu tetraploid; Morales never fractionated
LABELS = ["polyphest_p50", "polyphest_p70", "polyphest_p90", "grampa", "grampa_iter", "mpsugar", "alloppnet",
          "padre", "newmethod"]
BASE = {"polyphest_p50": 0.1, "polyphest_p70": 0.5, "polyphest_p90": 0.6, "grampa": 0.5, "grampa_iter": 0.45,
        "mpsugar": 0.35, "alloppnet": 0.8, "padre": 0.1, "newmethod": 0.3}
METRICS = ["ret_leaf_jaccard.dist", "ret_sisters_jaccard.dist", "mu_distance", "num_rets_bias"]
CONFIGS = common.BASELINE + [c for _, c in common.FRACTIONATION[1:]]
P50_FAILED = ("conf_ils_low_10M", "Ding_2023", 1)    # polyphest p50 failed, p70 completed
ITER_FAILED = ("conf_ils_low_10M", "Wu_2015", 2)     # grampa_iter failed, grampa completed
EMP_NETS = {
    "polyphest_p50": "(((a,b),(c,b)),(d,(e,f)));",
    "polyphest_p70": "(((a,b),(c,b)),((d,e),(f,e)));",
    "grampa_iter": "(((a,b),(c,b)),((d,a),(e,f)));",
    "padre": "((a,b),(c,(d,(e,f))));",
    "newmethod": "(((a,c),(b,c)),(d,(e,f)));",
}
EMP_DATASETS = ["Popp_2005", "Wu_2015", "Hori_2014"]
FALLBACK_VALUE = 0.99  # grampa_iter's partial_match differs from grampa's, so fallback.csv has its own values


def status_of(label, config, network, rep):
    if (config, network, rep) == P50_FAILED and label == "polyphest_p50":
        return "failed"
    if (config, network, rep) == ITER_FAILED and label == "grampa_iter":
        return "failed"
    if label == "alloppnet" and network not in common.TETRAPLOID_NETWORKS:
        return "missing"
    if label == "padre" and rep == 2:
        return "timeout"
    return "completed"


def value(rng, label, metric):
    if metric == "num_rets_bias":
        return float(rng.integers(-3, 4))
    return BASE[label] + rng.uniform(0, 0.2)


def build_results(root: Path) -> Path:
    rng = np.random.default_rng(0)
    inv, cmp = [], []
    for config in CONFIGS:
        nets = [n for n in NETWORKS if not (config.endswith(("075", "050", "025")) and n in common.FRACTIONATION_EXCLUDE)]
        for network in nets:
            for rep in (1, 2):
                for label in LABELS:
                    st = status_of(label, config, network, rep)
                    inv.append(("simulated", config, network, rep, label, st, 10.0))
                    if st == "completed":
                        cmp += [("simulated", config, network, rep, label, m, value(rng, label, m), "SUCCESS")
                                for m in METRICS]
    pw = []
    for ds in EMP_DATASETS:
        done = []
        for label in EMP_NETS:
            st = "failed" if (label == "polyphest_p50" and ds == "Hori_2014") else "completed"
            inv.append(("empirical", None, ds, None, label, st, 10.0))
            if st == "completed":
                d = root / label / "empirical" / ds
                d.mkdir(parents=True)
                (d / "network.nwk").write_text(EMP_NETS[label] + "\n")
                done.append(label)
        for i, a in enumerate(done):
            for b in done[i + 1:]:
                pw += [(ds, a, b, m, rng.uniform(0, 1)) for m in METRICS[:3]]
                pw += [(ds, a, b, "mu_scored", 1.0)]
    summary = root / "summary"
    summary.mkdir(parents=True)
    pd.DataFrame(inv, columns=["kind", "config", "network", "replicate", "label", "status", "runtime_s"]) \
        .to_csv(summary / "inventory.csv", index=False)
    cols = ["kind", "config", "network", "replicate", "label", "metric", "value", "status"]
    pd.DataFrame(cmp, columns=cols).to_csv(summary / "comparisons.csv", index=False)
    fb = [("simulated", *ITER_FAILED, "grampa_iter", m, FALLBACK_VALUE, "FALLBACK") for m in METRICS]
    pd.DataFrame(fb, columns=cols).to_csv(summary / "fallback.csv", index=False)
    pd.DataFrame(pw, columns=["dataset", "label_a", "label_b", "metric", "value"]) \
        .to_csv(summary / "empirical_pairwise.csv", index=False)
    nrf = [(c, n, r, rng.uniform(0, 0.5)) for c in common.ILS_ONLY for n in NETWORKS for r in (1, 2)]
    pd.DataFrame(nrf, columns=["config", "network", "replicate", "mean_nrf"]).to_csv(summary / "nrf.csv", index=False)
    return summary


@pytest.fixture(scope="module")
def paper(tmp_path_factory):
    root = tmp_path_factory.mktemp("paper")
    summary = build_results(root / "results")
    nets = root / "model_networks"
    nets.mkdir()
    for n in ("Ding_2023", "Shahrestani_2015"):
        shutil.copy(MODEL_NETWORKS / f"{n}.tre", nets)
    out = root / "out"
    return summary, nets, out, make_all(summary, nets, out)


EXPECTED = ["fig3.pdf", "fig3.png", "table3.csv", "table2.csv", "fig5.pdf", "fig5.png", "fig5_completion.csv",
            "fig6.pdf", "fig6.png", "fig6.csv", "fig7.pdf", "fig7.png", "fig7.csv", "table4.csv",
            "simulation_tests.csv", "simulation_tests.txt", "accuracy.csv", "section33.txt",
            "section33_per_dataset.csv", "table_s2.csv"]


def test_make_all_writes_every_file(paper):
    _, _, out, paths = paper
    assert sorted(p.name for p in paths) == sorted(EXPECTED)
    for name in EXPECTED:
        assert (out / name).is_file() and (out / name).stat().st_size > 0, name


def test_newmethod_in_table4_and_fig5(paper):
    _, _, out, _ = paper
    t4 = pd.read_csv(out / "table4.csv")
    assert "newmethod" in set(t4["label"])
    assert "grampa" not in set(t4["label"])  # Table 4 leaves GRAMPA out, as the manuscript
    f5 = pd.read_csv(out / "fig5_completion.csv")
    assert "newmethod" in set(f5["label"])
    acc = pd.read_csv(out / "accuracy.csv")
    assert "newmethod" in set(acc["label"])
    tests = pd.read_csv(out / "simulation_tests.csv")
    assert tests["block"].str.contains("newmethod").any()


def test_method_order_starts_with_paper_order(paper):
    _, _, out, _ = paper
    f5 = pd.read_csv(out / "fig5_completion.csv")
    assert list(f5["label"]) == ["polyphest", "grampa", "grampa_iter", "mpsugar", "alloppnet", "padre", "newmethod"]
    t4 = pd.read_csv(out / "table4.csv")
    assert list(t4["label"]) == ["polyphest", "mpsugar", "grampa_iter", "alloppnet", "padre", "newmethod"]
    assert list(t4["method"][:5]) == ["Polyphest", "MPAllopp", "GRAMPAIter", "AlloppNET", "PADRE"]
    assert common.order(["zeta", "padre", "alpha", "polyphest"]) == ["polyphest", "padre", "alpha", "zeta"]
    c = common.colors(["polyphest", "grampa_iter", "newmethod"])
    assert c["polyphest"] == "#DE8F05" and c["grampa_iter"] == "#56B4E9" and c["newmethod"] == common.EXTRA_COLORS[0]


def test_fig5_completion_values(paper):
    _, _, out, _ = paper
    f5 = pd.read_csv(out / "fig5_completion.csv").set_index("label")
    # p50 failed but p70 completed, so polyphest completed every run of that config
    assert f5.loc["polyphest", "conf_ils_low_10M"] == 100
    # one of six grampa_iter runs failed; padre timed out on every rep2; alloppnet never run on Morales
    assert f5.loc["grampa_iter", "conf_ils_low_10M"] == pytest.approx(500 / 6)
    assert f5.loc["padre", "conf_ils_low_10M"] == 50
    assert f5.loc["alloppnet", "conf_ils_low_10M"] == pytest.approx(400 / 6)


def test_polyphest_lowest_completed_percentile(paper):
    summary, _, _, _ = paper
    inv, cmp = common.load_simulated(summary)
    assert not inv["label"].str.startswith("polyphest_").any()
    raw = pd.read_csv(summary / "comparisons.csv")
    key = dict(zip(["config", "network", "replicate"], P50_FAILED))

    def pick(df, label):
        m = (df["label"] == label) & (df["metric"] == "mu_distance")
        for k, v in key.items():
            m &= df[k] == v
        return df.loc[m, "value"].tolist()

    assert pick(cmp, "polyphest") == pick(raw, "polyphest_p70")
    assert inv.set_index(["label", "config", "network", "replicate"]).loc[("polyphest", *P50_FAILED), "status"] \
        == "completed"
    other = ("conf_ils_low_10M", "Ding_2023", 2)
    key = dict(zip(["config", "network", "replicate"], other))
    assert pick(cmp, "polyphest") == pick(raw, "polyphest_p50")


def test_polyphest_none_completed_counts_as_failed():
    inv = pd.DataFrame({"config": ["c"] * 3, "network": ["n"] * 3, "replicate": [1] * 3,
                        "label": ["polyphest_p90", "polyphest_p50", "polyphest_p70"],
                        "status": ["timeout", "failed", "timeout"]})
    merged = common.merge_polyphest_inventory(inv)
    assert merged[["label", "status"]].values.tolist() == [["polyphest", "failed"]]


def test_grampa_iter_fallback(paper):
    summary, _, _, _ = paper
    key = dict(zip(["config", "network", "replicate"], ITER_FAILED))

    def rows(df, label):
        m = (df["label"] == label) & (df["metric"] == "mu_distance")
        for k, v in key.items():
            m &= df[k] == v
        return df.loc[m, "value"].tolist()

    inv, cmp = common.load_simulated(summary, grampa_iter_fallback=False)
    assert rows(cmp, "grampa_iter") == []
    inv_fb, cmp_fb = common.load_simulated(summary)  # default, as the paper
    assert rows(cmp_fb, "grampa_iter") == [FALLBACK_VALUE]  # scored by evaluate with grampa_iter's partial_match
    assert rows(cmp_fb, "grampa") not in ([], [FALLBACK_VALUE])
    for df in (inv, inv_fb):  # completion still counts the failed run
        st = df.set_index(["label", "config", "network", "replicate"]).loc[("grampa_iter", *ITER_FAILED), "status"]
        assert st == "failed"


def test_legacy_label_names(tmp_path):
    summary = tmp_path / "summary"
    summary.mkdir()
    pd.DataFrame({"kind": ["simulated"], "config": ["c"], "network": ["n"], "replicate": [1],
                  "label": ["grandma_split"], "status": ["completed"], "runtime_s": [1.0]}) \
        .to_csv(summary / "inventory.csv", index=False)
    assert common.read(summary, "inventory.csv")["label"].tolist() == ["grampa_iter"]


def test_table3_matches_manuscript():
    t = table3().set_index("measure")
    manuscript = {"Reticulation descendants": (1.0, 0, 0.5), "Reticulation sisters": (0, 0.5, 0.5),
                  "Difference in reticulation count (signed)": (0, 0, 1), "mu-distance": (1.0, 0.43, 0.38)}
    for measure, vals in manuscript.items():
        got = t.loc[measure, ["Wrong polyploid", "Wrong parental lineage", "Extra reticulation"]].tolist()
        assert [round(v + 1e-9, 2) for v in got] == list(vals), measure


# legacy run_reticulation_stats.py run with the current legacy reticulate_tree.py (identical to the port). The
# committed legacy mul_tree_final_stats.csv predates later folding fixes and differs on 15 networks. Morales-Briones
# matches the manuscript: 82 species, 49 polyploidization events.
LEGACY_S2 = """Filename,Num_Species,Num_Polyploids,Max_Copies,H_Strict,H_Relaxed,H_Diff,Num_Autopolyploidization_Events,Total_WGD
Diaz-Perez_2018.tre,18,9,3,10,13,3,3,13
Ding_2023.tre,12,7,2,0,1,1,1,1
Morales-Briones_2021.tre,82,49,3,49,45,-4,0,49
Popp_2005.tre,13,3,3,2,4,2,2,4
Wu_2015.tre,18,9,2,8,8,0,0,8
"""


def test_table_s2_matches_legacy(tmp_path):
    legacy = pd.read_csv(io.StringIO(LEGACY_S2))
    for f in legacy["Filename"]:
        shutil.copy(MODEL_NETWORKS / f, tmp_path)
    got = network_stats(tmp_path).set_index("Filename")
    for _, row in legacy.iterrows():
        for col in legacy.columns[1:]:
            assert got.loc[row["Filename"], col] == row[col], (row["Filename"], col)


def test_missing_config_fails_loud(paper, tmp_path):
    summary, nets, _, _ = paper
    s = tmp_path / "summary"
    s.mkdir()
    for f in ("inventory.csv", "comparisons.csv", "fallback.csv"):
        df = pd.read_csv(summary / f)
        df[df["config"] != "conf_dup_loss_high_10M_ne2M"].to_csv(s / f, index=False)
    for module in (fig5, fig6_fig7, stats):
        with pytest.raises(ValueError, match="conf_dup_loss_high_10M_ne2M"):
            module.make(s, nets, tmp_path / "out")


def test_missing_columns_fail_loud(tmp_path):
    (tmp_path / "inventory.csv").write_text("kind,config,network\nsimulated,c,n\n")
    with pytest.raises(ValueError, match="label"):
        common.read(tmp_path, "inventory.csv")
    with pytest.raises(FileNotFoundError):
        common.read(tmp_path, "comparisons.csv")


def test_fig6_methods_override(paper, tmp_path):
    summary, nets, _, _ = paper
    fig6_fig7.make(summary, nets, tmp_path, methods=["polyphest", "grampa_iter", "newmethod"])
    assert set(pd.read_csv(tmp_path / "fig6.csv")["label"]) == {"polyphest", "grampa_iter", "newmethod"}
    f7 = pd.read_csv(tmp_path / "fig7.csv")
    assert "Morales-Briones_2021" not in f7.columns
    assert (f7["n_networks"] == 2).all()  # Morales dropped from the whole series, 1.00 baseline included


def test_holm():
    assert stats.holm([0.01, 0.04, 0.03]).tolist() == pytest.approx([0.03, 0.06, 0.06])
    assert stats.holm([0.5, 0.9]).tolist() == pytest.approx([1.0, 1.0])


def test_table4_counts(paper):
    _, _, out, _ = paper
    t4 = pd.read_csv(out / "table4.csv").set_index("label")
    # 2 tetraploid networks x 2 reps x 3 ILS configs in the synthetic data
    assert t4.loc["polyphest", "completed"] == 12 and t4.loc["polyphest", "total"] == 12
    assert t4.loc["padre", "completed"] == 6
    assert t4.loc["grampa_iter", "completed"] == 11


def test_fallback_flag_changes_outputs(paper, tmp_path):
    summary, nets, out, _ = paper
    table4.make(summary, nets, tmp_path / "off", grampa_iter_fallback=False)
    on = pd.read_csv(out / "table4.csv").set_index("label").loc["grampa_iter"]  # default is on
    off = pd.read_csv(tmp_path / "off" / "table4.csv").set_index("label").loc["grampa_iter"]
    assert on["completed"] == off["completed"] == 11  # completion still counts the failed run
    assert on["mu_distance"] != off["mu_distance"]
    make_all(summary, nets, tmp_path / "all_off", grampa_iter_fallback=False)
    assert pd.read_csv(tmp_path / "all_off" / "table4.csv").set_index("label").loc["grampa_iter", "mu_distance"] \
        == off["mu_distance"]


def test_tables_use_fallback_csv(paper, tmp_path):
    summary, nets, out, _ = paper
    copy = tmp_path / "summary"
    shutil.copytree(summary, copy)
    fb = pd.read_csv(copy / "fallback.csv")
    fb.loc[fb["metric"] == "mu_distance", "value"] = 0.0
    fb.to_csv(copy / "fallback.csv", index=False)
    table4.make(copy, nets, tmp_path / "t4")
    before = pd.read_csv(out / "table4.csv").set_index("label").loc["grampa_iter", "mu_distance"]
    after = pd.read_csv(tmp_path / "t4" / "table4.csv").set_index("label").loc["grampa_iter", "mu_distance"]
    assert after < before
    (copy / "fallback.csv").unlink()
    with pytest.raises(FileNotFoundError, match="fallback.csv"):
        common.load_simulated(copy)
    common.load_simulated(copy, grampa_iter_fallback=False)


def test_table2_nrf_only_on_ils_only_rows(paper):
    _, _, out, _ = paper
    t2 = pd.read_csv(out / "table2.csv")
    no_dup = t2["dup_loss"] == "none"
    assert t2.loc[no_dup, "nrf_mean"].notna().all() and len(t2) == 12
    assert t2.loc[~no_dup, ["nrf_mean", "nrf_min", "nrf_max", "n_networks"]].isna().all().all()


def test_accuracy_bias_is_a_count(paper):
    _, _, out, _ = paper
    acc = pd.read_csv(out / "accuracy.csv")
    assert "num_rets_bias_count_mean" in acc.columns and "num_rets_bias_mean" not in acc.columns


def test_section33_fails_loud_without_pairs(paper, tmp_path):
    from polybench.paper import section33
    summary, nets, _, _ = paper
    root = tmp_path / "results"
    shutil.copytree(summary.parent, root)
    pw = pd.read_csv(root / "summary" / "empirical_pairwise.csv")
    drop = (pw["label_a"] == "grampa_iter") | (pw["label_b"] == "grampa_iter")
    pw[~drop].to_csv(root / "summary" / "empirical_pairwise.csv", index=False)
    with pytest.raises(ValueError, match="no Polyphest vs GRAMPAIter"):
        section33.make(root / "summary", nets, tmp_path / "out")
    pd.read_csv(summary / "empirical_pairwise.csv").to_csv(root / "summary" / "empirical_pairwise.csv", index=False)
    inv = pd.read_csv(root / "summary" / "inventory.csv")
    emp_iter = (inv["kind"] == "empirical") & (inv["label"] == "grampa_iter")
    inv.loc[emp_iter, "status"] = "failed"
    inv.to_csv(root / "summary" / "inventory.csv", index=False)
    with pytest.raises(ValueError, match="no completed empirical runs"):
        section33.make(root / "summary", nets, tmp_path / "out")
