from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from polybench.evaluate.compare import PUBLISHED, load_network
from polybench.paper.common import _variants, merge_polyphest_inventory, name, read, write_csv

P, G = "polyphest", "grampa_iter"
TRIO = {"mu_distance": "mu-distance", "ret_leaf_jaccard.dist": "ret. descendants",
        "ret_sisters_jaccard.dist": "ret. sister"}
POPP, EPHEDRA = "Popp_2005", "Wu_2015"


def empirical(summary_dir):
    # polyphest is the lowest percentile that completed per dataset; pairs with the other percentiles are dropped
    inv = read(summary_dir, "inventory.csv")
    inv = inv[inv["kind"] == "empirical"]
    variants = _variants(inv["label"])
    chosen = inv[inv["label"].isin(variants) & (inv["status"] == "completed")]
    chosen = chosen.assign(_r=chosen["label"].map(variants.index)).sort_values("_r").drop_duplicates("network")
    chosen = dict(zip(chosen["network"], chosen["label"]))
    pw = read(summary_dir, "empirical_pairwise.csv")
    keep = pd.Series(True, index=pw.index)
    for col in ("label_a", "label_b"):
        is_var = pw[col].isin(variants)
        keep &= ~is_var | (pw[col] == pw["dataset"].map(chosen))
        pw.loc[is_var, col] = P
    inv = merge_polyphest_inventory(inv, keys=["network"])
    inv["path_label"] = inv.apply(lambda r: chosen.get(r["network"], r["label"]) if r["label"] == P else r["label"],
                                  axis=1)
    return inv, pw[keep & (pw["label_a"] != PUBLISHED) & (pw["label_b"] != PUBLISHED)]


def pair(pw, a, b, metric) -> pd.Series:
    rows = pw[(((pw["label_a"] == a) & (pw["label_b"] == b)) | ((pw["label_a"] == b) & (pw["label_b"] == a)))
              & (pw["metric"] == metric)]
    return rows.set_index("dataset")["value"]


def counts(inv, results_root) -> pd.DataFrame:
    # reticulations, polyploid species and mean descendant species per reticulation, per label and dataset
    rows = []
    for r in inv[inv["label"].isin([P, G]) & (inv["status"] == "completed")].itertuples():
        path = Path(results_root) / r.path_label / "empirical" / r.network / "network.nwk"
        if not path.is_file():
            raise FileNotFoundError(f"{path} not found for a completed run")
        net = load_network(path)
        leaves = net.get_reticulation_leaves()
        rows.append({"dataset": r.network, "label": r.label, "reticulations": net.get_reticulation_count(),
                     "polyploids": sum(v > 1 for v in net.get_leaf_counts().values()),
                     "mean_descendants": np.mean([len(set(v)) for v in leaves.values()]) if leaves else np.nan})
    return pd.DataFrame(rows, columns=["dataset", "label", "reticulations", "polyploids", "mean_descendants"])


def paired_report(tab, what, w):
    both = tab[[P, G]].dropna() if {P, G} <= set(tab.columns) else pd.DataFrame(columns=[P, G])
    diff = both[G] - both[P]
    w(f"--- {what} (n = {len(both)} datasets) ---")
    for m in (P, G):
        w(f"  {name(m):<11} mean {both[m].mean():7.2f}   median {both[m].median():7.2f}")
    w(f"  {name(G)} higher on {(diff > 0).sum()}, lower on {(diff < 0).sum()}, equal on {(diff == 0).sum()}")
    if (diff != 0).sum() < 2:
        w("  Wilcoxon: not computable (fewer than two non-zero differences)")
    else:
        stat, p = wilcoxon(both[P], both[G], alternative="two-sided")
        w(f"  Wilcoxon signed-rank, two-sided: W = {stat:.1f}, p = {p:.4f}")
    w("")


def make(summary_dir, networks_dir, out_dir) -> list:
    inv, pw = empirical(summary_dir)
    out = []
    w = out.append

    w(f"--- Popp: Polyphest vs PADRE on {POPP} ---")
    for metric in ("ret_sisters_jaccard.dist", "mu_distance"):
        v = pair(pw, P, "padre", metric)
        w(f"  {TRIO[metric]:<18}: " + (f"{v[POPP]:.2f}" if POPP in v.index else "not available"))
    w("")

    w(f"--- {EPHEDRA}: every method pair ---")
    ep = pw[pw["dataset"] == EPHEDRA]
    tab = ep[ep["metric"].isin(TRIO)].pivot_table(index=["label_a", "label_b"], columns="metric", values="value")
    if tab.empty:
        w("  not available")
    else:
        for (a, b), r in tab.iterrows():
            w(f"  {name(a) + ' vs ' + name(b):<32}" + "".join(f"{TRIO[m]} {r.get(m, np.nan):.3f}   " for m in TRIO))
        mu = tab["mu_distance"].dropna() if "mu_distance" in tab else pd.Series(dtype=float)
        if len(mu):
            w(f"  mu-distance range {mu.min():.2f} to {mu.max():.2f} over {len(mu)} pairs")
        for m in TRIO:
            if m in tab and tab[m].notna().any():
                a, b = tab[m].idxmin()
                w(f"  closest pair on {TRIO[m]:<18}: {name(a)} vs {name(b)} ({tab[m].min():.3f})")
    w("")

    w(f"--- {name(P)} vs {name(G)} across datasets ---")
    per_ds = pd.DataFrame({m: pair(pw, P, G, m) for m in [*TRIO, "mu_scored"]})
    missing = [m for m in [*TRIO, "mu_scored"] if per_ds.get(m, pd.Series(dtype=float)).dropna().empty]
    if missing:
        raise ValueError(f"Section 3.3: no {name(P)} vs {name(G)} empirical pairs with {missing}")
    for m in TRIO:
        s = per_ds[m].dropna() if m in per_ds else pd.Series(dtype=float)
        w(f"  {TRIO[m]:<18} min {s.min():.2f}  max {s.max():.2f}  mean {s.mean():.2f}  median {s.median():.2f}"
          f"  sd {s.std():.2f}  (n = {len(s)})")
        if len(s):
            w(f"  {'':<18} lowest on {s.idxmin()}, highest on {s.idxmax()}")
    scored = per_ds["mu_scored"].dropna()
    w(f"  mu-distance scored on {int((scored == 1).sum())} of {len(scored)} datasets")
    w("")

    cnt = counts(inv, Path(summary_dir).parent)
    absent = [lab for lab in (P, G) if lab not in set(cnt["label"])]
    if absent:
        raise ValueError(f"Section 3.3: no completed empirical runs for {absent}")
    for col, what in (("polyploids", "Number of polyploid species"), ("reticulations", "Number of reticulation events"),
                      ("mean_descendants", "Mean descendant species per reticulation (larger = deeper)")):
        paired_report(cnt.pivot_table(index="dataset", columns="label", values=col, aggfunc="first"), what, w)

    w("--- Spearman correlations across datasets ---")
    for a, b in (("mu_distance", "ret_leaf_jaccard.dist"), ("mu_distance", "ret_sisters_jaccard.dist"),
                 ("ret_leaf_jaccard.dist", "ret_sisters_jaccard.dist")):
        both = per_ds[[a, b]].dropna()
        if len(both) < 3:
            w(f"  {TRIO[a]:<18} vs {TRIO[b]:<18} not computable (n = {len(both)})")
            continue
        rho, p = spearmanr(both[a], both[b])
        w(f"  {TRIO[a]:<18} vs {TRIO[b]:<18} rho = {rho:6.3f}, p = {p:.4f}, n = {len(both)}")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "section33.txt").write_text("\n".join(out) + "\n")
    table = per_ds.rename_axis("dataset").reset_index()
    return [out_dir / "section33.txt", write_csv(table, out_dir, "section33_per_dataset.csv")]
