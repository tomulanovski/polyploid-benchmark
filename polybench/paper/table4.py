import pandas as pd

from polybench.paper.common import (ILS_ONLY, TETRAPLOID_NETWORKS, TRIO, load_simulated, name, order, per_network,
                                    require_configs, success, write_csv)

TABLE4_ORDER = ["polyphest", "mpsugar", "grampa_iter", "alloppnet", "padre"]  # manuscript Table 4 rows


def table4(inv: pd.DataFrame, cmp: pd.DataFrame) -> pd.DataFrame:
    # 8 tetraploid networks x 5 replicates x 3 ILS-only configs; GRAMPA is not in the table (legacy TETRA_METHODS)
    require_configs(inv, ILS_ONLY, "Table 4")
    inv = inv[inv["config"].isin(ILS_ONLY) & inv["network"].isin(TETRAPLOID_NETWORKS)]
    cmp = cmp[cmp["config"].isin(ILS_ONLY) & cmp["network"].isin(TETRAPLOID_NETWORKS)]
    if inv.empty:
        raise ValueError(f"Table 4: no datasets of the tetraploid networks {TETRAPLOID_NETWORKS}")
    labels = [lab for lab in order(inv["label"]) if lab != "grampa"]
    labels = [lab for lab in TABLE4_ORDER if lab in labels] + [lab for lab in labels if lab not in TABLE4_ORDER]
    rows = []
    for label in labels:
        runs = inv[inv["label"] == label]
        done = int((runs["status"] == "completed").sum())
        row = {"method": name(label), "label": label, "completed": done, "total": len(runs),
               "completed_pct": 100 * done / len(runs)}
        for metric in TRIO:
            # replicate values pooled per network across the three configs, then the mean over networks
            row[metric] = per_network(success(cmp, metric).query("label == @label")).mean()
        rows.append(row)
    return pd.DataFrame(rows)


def make(summary_dir, networks_dir, out_dir, grampa_iter_fallback: bool = True) -> list:
    return [write_csv(table4(*load_simulated(summary_dir, grampa_iter_fallback)), out_dir, "table4.csv")]
