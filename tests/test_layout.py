import pytest

from polybench.config import load_config
from polybench.layout import Dataset, find_datasets, result_dir, true_network_path


@pytest.fixture
def root(tmp_path):
    for rep in (1, 2):
        d = tmp_path / "simulated" / "conf_a" / "Ding_2023" / f"rep{rep}"
        d.mkdir(parents=True)
        (d / "gene_trees.nwk").write_text("(a,b);\n")
    e = tmp_path / "empirical" / "Wheat"
    e.mkdir(parents=True)
    (e / "gene_trees.nwk").write_text("(a,b);\n")
    (tmp_path / "model_networks").mkdir()
    (tmp_path / "model_networks" / "Ding_2023.tre").write_text("(a,b);")
    return tmp_path


def test_find_simulated(root):
    ds = find_datasets(root, "simulated/*/*/*")
    assert [d.id for d in ds] == [
        "simulated/conf_a/Ding_2023/rep1",
        "simulated/conf_a/Ding_2023/rep2",
    ]
    assert ds[0] == Dataset("simulated", "conf_a", "Ding_2023", 1, root / "simulated/conf_a/Ding_2023/rep1", root)


def test_find_empirical(root):
    ds = find_datasets(root, "empirical/*")
    assert [d.id for d in ds] == ["empirical/Wheat"]
    assert ds[0].kind == "empirical"


def test_true_network_path(root):
    sim = find_datasets(root, "simulated/*/*/*")[0]
    assert true_network_path(root, sim) == root / "model_networks" / "Ding_2023.tre"
    assert true_network_path(root, find_datasets(root, "empirical/*")[0]) is None


def test_missing_gene_trees_is_left_to_validate(root):
    (root / "empirical" / "Wheat" / "gene_trees.nwk").unlink()
    assert [d.id for d in find_datasets(root, "empirical/*")] == ["empirical/Wheat"]


def test_bad_layout_raises(root):
    (root / "empirical" / "Wheat" / "extra").mkdir()
    with pytest.raises(ValueError, match="is not <data_root>"):
        find_datasets(root, "empirical/*/*")


def test_no_match_raises(root):
    with pytest.raises(FileNotFoundError, match="no datasets"):
        find_datasets(root, "simulated/nope/*/*")


def test_result_dir(root, tmp_path):
    ds = find_datasets(root, "empirical/*")[0]
    assert result_dir(tmp_path, "padre", ds) == tmp_path / "padre" / "empirical" / "Wheat"


def test_load_config(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("data_root: /d\nresults_root: /r\ntools: {astral: /a.jar}\nexecutor: {type: local}\n")
    cfg = load_config(p)
    assert cfg["tools"]["astral"] == "/a.jar"
    assert cfg["executor"] == {"type": "local", "slurm_header": [], "conda_activate": ""}


def test_load_config_resolves_relative_paths(tmp_path, monkeypatch):
    (tmp_path / "conf").mkdir()
    p = tmp_path / "conf" / "c.yaml"
    p.write_text("data_root: data\nresults_root: ../res\ntools: {astral: bin/astral, simphy: simphy, iqtree: null}\n"
                 "executor: {type: local}\n")
    monkeypatch.chdir(tmp_path)
    cfg = load_config("conf/c.yaml")
    assert cfg["data_root"] == str(tmp_path / "conf" / "data")
    assert cfg["results_root"] == str(tmp_path / "res")
    assert cfg["tools"] == {"astral": str(tmp_path / "conf" / "bin" / "astral"), "simphy": "simphy", "iqtree": None}


def test_load_config_bad_executor(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("data_root: /d\nresults_root: /r\ntools: {}\nexecutor: {type: pbs}\n")
    with pytest.raises(ValueError, match="executor.type"):
        load_config(p)


def test_load_config_missing_key(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("data_root: /d\n")
    with pytest.raises(ValueError, match="results_root"):
        load_config(p)


def test_load_config_accepts_new_tool_keys(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("data_root: d\nresults_root: r\ntools: {my_tool2: bin/mt}\nexecutor: {type: local}\n")
    assert load_config(p)["tools"]["my_tool2"] == str(tmp_path / "bin" / "mt")
    p.write_text("data_root: d\nresults_root: r\ntools: {my-tool: x}\nexecutor: {type: local}\n")
    with pytest.raises(ValueError, match="my-tool"):
        load_config(p)
