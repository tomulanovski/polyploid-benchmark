import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from polybench.evaluate.compare import COMPARISON_COLS, INVENTORY_COLS, PAIRWISE_COLS  # noqa: E402

# Fig 5 column groups: family -> {level: config}
FAMILIES = {
    "ILS only": {"Low": "conf_ils_low_10M", "Medium": "conf_ils_medium_10M", "High": "conf_ils_high_10M"},
    "Dup/Loss\n(low ILS)": {"Low": "conf_dup_loss_low_10M", "Medium": "conf_dup_loss_medium_10M",
                            "High": "conf_dup_loss_high_10M"},
    "Dup/Loss\n(med ILS)": {"Low": "conf_dup_loss_low_10M_ne1M", "Medium": "conf_dup_loss_medium_10M_ne1M",
                            "High": "conf_dup_loss_high_10M_ne1M"},
    "Dup/Loss\n(high ILS)": {"Low": "conf_dup_loss_low_10M_ne2M", "Medium": "conf_dup_loss_medium_10M_ne2M",
                             "High": "conf_dup_loss_high_10M_ne2M"},
}
BASELINE = [c for fam in FAMILIES.values() for c in fam.values()]
ILS_ONLY = list(FAMILIES["ILS only"].values())
ILS_LEVELS = ["Low", "Medium", "High"]
# Fig 6 columns: dup/loss level -> {ILS level: config}
DUP_PANELS = {
    "No Dup/Loss": dict(zip(ILS_LEVELS, ILS_ONLY)),
    **{f"{lvl} Dup/Loss": {ils: FAMILIES[fam][lvl] for ils, fam in zip(ILS_LEVELS, list(FAMILIES)[1:])}
       for lvl in ILS_LEVELS},
}
# (retention, config); the 1.00 baseline is the medium ILS, medium dup/loss config
FRACTIONATION = [("1.00", "conf_dup_loss_medium_10M_ne1M"),
                 ("0.75", "conf_dup_loss_medium_10M_ne1M_fix075"),
                 ("0.50", "conf_dup_loss_medium_10M_ne1M_fix050"),
                 ("0.25", "conf_dup_loss_medium_10M_ne1M_fix025")]
FRACTIONATION_EXCLUDE = ["Morales-Briones_2021"]  # never simulated under fractionation
TETRAPLOID_NETWORKS = ["Bendiksby_2011", "Ding_2023", "Koenen_2020", "Liu_2023", "Shahrestani_2015",
                       "Wisecaver_2023", "Wu_2015", "Zhao_2021"]

TRIO = {"ret_leaf_jaccard.dist": "Ret. Descendants Measure",
        "ret_sisters_jaccard.dist": "Ret. Sisters Measure",
        "mu_distance": r"$\mu$-distance"}

METHOD_ORDER = ["polyphest", "grampa", "grampa_iter", "mpsugar", "alloppnet", "padre"]
DISPLAY = {"polyphest": "Polyphest", "grampa": "GRAMPA", "grampa_iter": "GRAMPAIter", "mpsugar": "MPAllopp",
           "alloppnet": "AlloppNET", "padre": "PADRE"}
PLOT_DISPLAY = {**DISPLAY, "grampa_iter": r"GRAMPA$^{Iter}$"}
COLORS = {"polyphest": "#DE8F05", "grampa": "#0173B2", "grampa_iter": "#56B4E9", "mpsugar": "#8B4513",
          "alloppnet": "#DC143C", "padre": "#ECE133"}
EXTRA_COLORS = ["#029E73", "#CC78BC", "#949494", "#FBAFE4", "#CA9161", "#000000"]
LEGACY_LABELS = {"grandma_split": "grampa_iter"}
POLYPHEST_VARIANT = re.compile(r"polyphest_p(\d+)$")

COLUMNS = {
    "inventory.csv": INVENTORY_COLS,
    "comparisons.csv": COMPARISON_COLS,
    "fallback.csv": COMPARISON_COLS,
    "empirical_pairwise.csv": PAIRWISE_COLS,
    "nrf.csv": ["config", "network", "replicate", "mean_nrf"],
}
SIM_KEYS = ["config", "network", "replicate"]


def order(labels) -> list[str]:
    present = set(labels)
    return [m for m in METHOD_ORDER if m in present] + sorted(present - set(METHOD_ORDER))


def colors(labels) -> dict[str, str]:
    extras = [lab for lab in order(labels) if lab not in COLORS]
    return {**{lab: COLORS[lab] for lab in labels if lab in COLORS},
            **{lab: EXTRA_COLORS[i % len(EXTRA_COLORS)] for i, lab in enumerate(extras)}}


def name(label: str, plot: bool = False) -> str:
    return (PLOT_DISPLAY if plot else DISPLAY).get(label, label)


def read(summary_dir, filename: str) -> pd.DataFrame:
    path = Path(summary_dir) / filename
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found")
    df = pd.read_csv(path)
    missing = [c for c in COLUMNS[filename] if c not in df.columns]
    if missing:
        raise ValueError(f"{path} lacks columns {missing}")
    for col in ("label", "label_a", "label_b"):
        if col in df.columns:
            df[col] = df[col].replace(LEGACY_LABELS)
    return df


def _variants(labels) -> list[str]:
    labels = set(labels)
    variants = sorted((lab for lab in labels if POLYPHEST_VARIANT.match(lab)),
                      key=lambda lab: int(POLYPHEST_VARIANT.match(lab)[1]))
    if variants and "polyphest" in labels:
        raise ValueError(f"both 'polyphest' and {variants} labels present, cannot merge")
    return variants


def merge_polyphest_inventory(inv: pd.DataFrame, keys=SIM_KEYS) -> pd.DataFrame:
    # per dataset the lowest percentile that completed; if none did, the lowest percentile row (a failure)
    variants = _variants(inv["label"])
    if not variants:
        return inv
    poly = inv[inv["label"].isin(variants)]
    rank = poly["label"].map(variants.index) + np.where(poly["status"] == "completed", 0, len(variants))
    best = poly.assign(_rank=rank).sort_values("_rank").drop_duplicates(keys).drop(columns="_rank")
    return pd.concat([inv[~inv["label"].isin(variants)], best.assign(label="polyphest")], ignore_index=True)


def merge_polyphest_comparisons(cmp: pd.DataFrame, keys=SIM_KEYS) -> pd.DataFrame:
    # per dataset the rows of the lowest percentile with a successful comparison, as legacy
    variants = _variants(cmp["label"])
    if not variants:
        return cmp
    poly = cmp[cmp["label"].isin(variants)].assign(_rank=lambda d: d["label"].map(variants.index))
    best = poly[poly["status"] == "SUCCESS"].groupby(keys)["_rank"].min().rename("_best").reset_index()
    poly = poly.merge(best, on=keys)
    poly = poly[poly["_rank"] == poly["_best"]].drop(columns=["_rank", "_best"])
    return pd.concat([cmp[~cmp["label"].isin(variants)], poly.assign(label="polyphest")], ignore_index=True)


def apply_grampa_iter_fallback(cmp: pd.DataFrame, fallback: pd.DataFrame) -> pd.DataFrame:
    # fallback.csv: GRAMPAIter runs that failed or timed out, scored by evaluate with the GRAMPA network of the
    # same dataset; completion is untouched
    if set(fallback["label"]) - {"grampa_iter"}:
        raise ValueError(f"fallback.csv has labels other than grampa_iter: {sorted(set(fallback['label']))}")
    fb = fallback.assign(status=fallback["status"].replace({"FALLBACK": "SUCCESS"}))
    return pd.concat([cmp, fb], ignore_index=True)


def load_simulated(summary_dir, grampa_iter_fallback: bool = True):
    # fallback on by default, as the paper: a failed GRAMPAIter run is scored with its GRAMPA result
    inv = read(summary_dir, "inventory.csv")
    inv = merge_polyphest_inventory(inv[inv["kind"] == "simulated"])
    cmp = read(summary_dir, "comparisons.csv")
    cmp = merge_polyphest_comparisons(cmp[cmp["kind"] == "simulated"])
    if grampa_iter_fallback:
        cmp = apply_grampa_iter_fallback(cmp, read(summary_dir, "fallback.csv"))
    return inv, cmp


def require_configs(df: pd.DataFrame, configs, what: str) -> None:
    missing = [c for c in configs if c not in set(df["config"])]
    if missing:
        raise ValueError(f"{what} needs configs {missing}, not found in the summary")


def require_labels(df: pd.DataFrame, labels, what: str) -> None:
    missing = [lab for lab in labels if lab not in set(df["label"])]
    if missing:
        raise ValueError(f"{what} needs results for {missing}, not found in the summary")


def require_metrics(cmp: pd.DataFrame, metrics, what: str) -> None:
    missing = [m for m in metrics if m not in set(cmp.loc[cmp["status"] == "SUCCESS", "metric"])]
    if missing:
        raise ValueError(f"{what} needs metrics {missing}, not found in the summary")


def success(cmp: pd.DataFrame, metric: str) -> pd.DataFrame:
    return cmp[(cmp["metric"] == metric) & (cmp["status"] == "SUCCESS")]


def per_network(df: pd.DataFrame) -> pd.Series:
    # replicates of one network are not independent, so every mean and SEM is taken across networks
    return df.groupby("network")["value"].mean()


def mean_sem(vals: pd.Series) -> tuple[float, float]:
    if len(vals) == 0:
        return np.nan, 0.0
    return vals.mean(), (vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0.0)


def save(fig, out_dir, stem: str) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = [out_dir / f"{stem}.pdf", out_dir / f"{stem}.png"]
    fig.savefig(paths[0], bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight", dpi=300)
    plt.close(fig)
    return paths


def write_csv(df: pd.DataFrame, out_dir, filename: str) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / filename, index=False)
    return out_dir / filename
