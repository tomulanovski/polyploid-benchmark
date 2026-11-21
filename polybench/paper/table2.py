import pandas as pd

from polybench.paper.common import FAMILIES, ILS_LEVELS, ILS_ONLY, read, require_configs, write_csv
from polybench.simulate.simphy import CONFIGS


def table2(nrf: pd.DataFrame) -> pd.DataFrame:
    # the 12 baseline configs; nRF (mean over replicates per network, then mean and range across networks, as
    # legacy aggregate_ad_reusable) only on the ILS-only rows, the configs it is computed for
    require_configs(nrf, ILS_ONLY, "Table 2")
    rows = []
    for ils, ils_config, family in zip(ILS_LEVELS, ILS_ONLY, list(FAMILIES.values())[1:]):
        per_net = nrf[nrf["config"] == ils_config].groupby("network")["mean_nrf"].mean()
        dup_configs = {"none": ils_config, **family}
        for dup, config in dup_configs.items():
            c = CONFIGS[config]
            nrf_cols = ({"nrf_mean": per_net.mean(), "nrf_min": per_net.min(), "nrf_max": per_net.max(),
                         "n_networks": len(per_net)} if dup == "none" else {})
            rows.append({"config": config, "ils": ils, "ne": c["ne"], "dup_loss": dup, "dup_rate": c["dup"],
                         "loss_rate": c["loss"], **nrf_cols})
    return pd.DataFrame(rows, columns=["config", "ils", "ne", "dup_loss", "dup_rate", "loss_rate", "nrf_mean",
                                       "nrf_min", "nrf_max", "n_networks"])


def make(summary_dir, networks_dir, out_dir) -> list:
    return [write_csv(table2(read(summary_dir, "nrf.csv")), out_dir, "table2.csv")]
