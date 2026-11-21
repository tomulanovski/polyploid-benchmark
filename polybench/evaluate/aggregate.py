import pandas as pd

KEYS = ["network", "config", "label", "metric"]


def aggregate(comparisons: pd.DataFrame) -> pd.DataFrame:
    # mean/std/min/max over replicates of every successful comparison, as legacy aggregated_metrics.csv
    ok = comparisons[comparisons["status"] == "SUCCESS"]
    out = ok.groupby(KEYS)["value"].agg(mean="mean", std="std", min="min", max="max",
                                        n_valid="count").reset_index()
    out["std"] = out["std"].fillna(0)  # one replicate
    return out[KEYS + ["mean", "std", "min", "max", "n_valid"]]
