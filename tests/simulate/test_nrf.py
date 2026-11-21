import pandas as pd
import pytest

from polybench.simulate.nrf import nrf, nrf_table

T = "((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"
OTHER = "((A:1,C:1):1,(B:1,(D:1,E:1):1):1);"


def test_nrf_identical_zero():
    assert nrf(T, T) == 0.0


def test_nrf_ignores_child_order():
    assert nrf("((E:1,D:1):1,(C:1,(B:1,A:1):1):1);", T) == 0.0


def test_nrf_different_topology_positive():
    assert 0 < nrf(OTHER, T) <= 1


def test_nrf_leaf_mismatch_fails():
    with pytest.raises(ValueError):
        nrf("((A,B),(C,D));", T)


def _rep(root, config, net, k, genes):
    d = root / "simulated" / config / net / f"rep{k}"
    d.mkdir(parents=True)
    (d / "simphy_gene_trees.nwk").write_text("\n".join(genes) + "\n")
    (d / "simphy_species_tree.nwk").write_text(T + "\n")
    (d / "gene_trees.nwk").write_text("x;\n")


def test_nrf_table(tmp_path):
    _rep(tmp_path, "conf_ils_low_10M", "N1", 1, [T, T])
    _rep(tmp_path, "conf_ils_low_10M", "N1", 2, [T, OTHER])
    out = tmp_path / "summary" / "nrf.csv"
    df = nrf_table(tmp_path, ["conf_ils_low_10M"], out_csv=out)
    assert list(df.columns) == ["config", "network", "replicate", "mean_nrf"]
    assert df.loc[df.replicate == 1, "mean_nrf"].item() == 0.0
    assert df.loc[df.replicate == 2, "mean_nrf"].item() == pytest.approx(nrf(OTHER, T) / 2)
    assert pd.read_csv(out).shape == (2, 4)


def test_nrf_table_rejects_dup_loss_config(tmp_path):
    with pytest.raises(ValueError):
        nrf_table(tmp_path, ["conf_dup_loss_high_10M"])
