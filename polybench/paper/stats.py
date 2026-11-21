from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from polybench.paper.common import (BASELINE, FRACTIONATION, FRACTIONATION_EXCLUDE, ILS_ONLY, METHOD_ORDER,
                                    TETRAPLOID_NETWORKS, TRIO, FAMILIES, load_simulated, name, order,
                                    per_network, require_configs, require_metrics, write_csv)

P, G, M = "polyphest", "grampa_iter", "mpsugar"
MEASURE = {"mu_distance": "mu-distance", "ret_leaf_jaccard.dist": "ret. descendants",
           "ret_sisters_jaccard.dist": "ret. sister"}
# slope series, (x, config): dup/loss none..high at each ILS level, and x = 1 - retention for fractionation
DUPLOSS_SERIES = {f"dup/loss, {ils.lower()} ILS": list(enumerate([ils_cfg, *fam.values()]))
                  for ils, ils_cfg, fam in zip(["Low", "Medium", "High"], ILS_ONLY, list(FAMILIES.values())[1:])}
RETENTION_SERIES = [(1 - float(r), c) for r, c in FRACTIONATION]


def holm(pvals) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    adj = np.empty(len(p))
    running = 0.0
    for rank, i in enumerate(np.argsort(p)):
        running = max(running, (len(p) - rank) * p[i])
        adj[i] = min(running, 1.0)
    return adj


def network_means(cmp, configs, metric) -> pd.DataFrame:
    # network x label table of means over replicates (and over the configs given)
    d = cmp[cmp["config"].isin(configs) & (cmp["metric"] == metric)]
    return d.groupby(["network", "label"])["value"].mean().unstack("label")


def paired(a: pd.Series, b: pd.Series, what: str, strict: bool = True) -> dict:
    # two-sided Wilcoxon signed-rank paired by network; a paper test that cannot run is an error, a test of
    # an added method is reported as not computable
    both = pd.concat([a, b], axis=1, keys=["a", "b"]).dropna()
    diff = both["b"] - both["a"]
    out = {"n": len(both), "mean_1": both["a"].mean(), "mean_2": both["b"].mean(),
           "first_lower_on": int((diff > 0).sum()), "W": np.nan, "p": np.nan, "note": ""}
    if (diff != 0).sum() < 2:
        if strict:
            raise ValueError(f"{what}: fewer than two non-zero differences")
        return {**out, "note": "not computable: fewer than two non-zero differences"}
    stat, p = wilcoxon(both["a"], both["b"], alternative="two-sided")
    return {**out, "W": stat, "p": p}


def column(tab, label, what):
    if label not in tab.columns:
        raise ValueError(f"{what}: no results for {label}")
    return tab[label]


def per_config(cmp, configs, m1, m2, strict=True) -> list:
    rows = []
    for cfg in configs:
        for metric, measure in MEASURE.items():
            tab = network_means(cmp, [cfg], metric)
            what = f"{cfg} {metric}"
            if not strict and (m1 not in tab.columns or m2 not in tab.columns):
                rows.append({"test": cfg, "measure": measure, "n": 0, "note": "no results"})
                continue
            rows.append({"test": cfg, "measure": measure,
                         **paired(column(tab, m1, what), column(tab, m2, what), what, strict)})
    return rows


def slopes(cmp, series, metric, label) -> pd.Series:
    # per-network least-squares slope against x, networks with a value at every level
    cols = {x: column(network_means(cmp, [cfg], metric), label, f"{cfg} {metric}") for x, cfg in series}
    tab = pd.DataFrame(cols).dropna()
    xs = np.array(list(cols), dtype=float)
    return tab.apply(lambda row: np.polyfit(xs, row.values, 1)[0], axis=1)


def slope_block(cmp, series_map, m1, m2) -> list:
    return [{"test": f"slope, {sname}", "measure": measure,
             **paired(slopes(cmp, series, metric, m1), slopes(cmp, series, metric, m2), f"{sname} {metric}")}
            for sname, series in series_map.items() for metric, measure in MEASURE.items()]


def tetraploid_block(cmp, m1, m2) -> list:
    rows = []
    for metric, measure in MEASURE.items():
        tab = network_means(cmp, ILS_ONLY, metric)
        tab = tab[tab.index.isin(TETRAPLOID_NETWORKS)]
        rows.append({"test": "tetraploid subset, ILS-only configs", "measure": measure,
                     **paired(column(tab, m1, "tetraploid"), column(tab, m2, "tetraploid"), f"tetraploid {metric}")})
    return rows


def block(rows, title, m1, m2) -> pd.DataFrame:
    # Holm within the block, separately for each measure
    res = pd.DataFrame(rows).assign(block=title, method_1=name(m1), method_2=name(m2))
    res["p_holm"] = np.nan
    for _, idx in res.groupby("measure").groups.items():
        ok = res.loc[idx, "p"].dropna()
        res.loc[ok.index, "p_holm"] = holm(ok.values)
    return res


def report(res: pd.DataFrame) -> str:
    out = ['Means are over networks. "1 lower on" counts networks where the first method has the lower value.',
           "For slope tests the means are mean slopes per step of x.", ""]
    for title, r in res.groupby("block", sort=False):
        n1, n2 = r["method_1"].iloc[0], r["method_2"].iloc[0]
        out.append(f"=== {title} ({len(r)} tests, Holm within block, separately per measure) ===")
        out.append(f'  {"test":<46}{"measure":<18}{"n":>3}{n1:>12}{n2:>12}{"1 lower on":>11}{"W":>7}'
                   f'{"p":>10}{"p_Holm":>10}')
        for _, x in r.iterrows():
            if pd.isna(x["p"]):
                out.append(f'  {x["test"]:<46}{x["measure"]:<18}{x["n"]:>3}  {x["note"]}')
                continue
            out.append(f'  {x["test"]:<46}{x["measure"]:<18}{x["n"]:>3}{x["mean_1"]:>12.3f}{x["mean_2"]:>12.3f}'
                       f'{f"{x.first_lower_on}/{x.n}":>11}{x["W"]:>7.1f}{x["p"]:>10.2e}{x["p_holm"]:>10.2e}')
        tested = r["p_holm"].dropna()
        out.append(f"  significant after Holm: {(tested < 0.05).sum()}/{len(tested)}"
                   + (f", max p_Holm {tested.max():.2e}" if len(tested) else ""))
        out.append("")
    return "\n".join(out)


def accuracy(inv, cmp, configs) -> pd.DataFrame:
    # per label and config: completion and the per-network mean of each measure, for every label in the summary
    rows = []
    for label in order(inv["label"]):
        for cfg in configs:
            runs = inv[(inv["label"] == label) & (inv["config"] == cfg)]
            row = {"method": name(label), "label": label, "config": cfg,
                   "completion_pct": 100 * (runs["status"] == "completed").mean() if len(runs) else np.nan}
            for metric in [*TRIO, "num_rets_bias"]:
                vals = network_means(cmp, [cfg], metric)
                vals = vals[label].dropna() if label in vals.columns else pd.Series(dtype=float)
                col = "num_rets_bias_count" if metric == "num_rets_bias" else metric  # signed raw count, not %
                row[f"{col}_mean"] = vals.mean()
                row[f"{col}_median"] = vals.median()
                row[f"{col}_n_networks"] = len(vals)
            rows.append(row)
    return pd.DataFrame(rows)


def make(summary_dir, networks_dir, out_dir, grampa_iter_fallback: bool = True) -> list:
    inv, cmp = load_simulated(summary_dir, grampa_iter_fallback)
    frac = [c for _, c in FRACTIONATION]
    require_configs(cmp, BASELINE + frac[1:], "simulation tests")
    require_metrics(cmp, TRIO, "simulation tests")
    cmp = cmp[cmp["status"] == "SUCCESS"]
    excluded = cmp["network"].isin(FRACTIONATION_EXCLUDE)
    if (excluded & cmp["config"].isin(frac[1:])).any():
        raise ValueError(f"{FRACTIONATION_EXCLUDE} have fractionation results but are listed as never simulated")
    cmp_frac = cmp[~excluded]
    blocks = [
        block(per_config(cmp, BASELINE, P, G), "A. Baseline configurations", P, G),
        block(per_config(cmp_frac, frac, P, G), "B. Fractionation levels (excluded networks dropped)", P, G),
        block(slope_block(cmp, DUPLOSS_SERIES, P, G), "C1. Sensitivity to duplication and loss (x = level 0..3)",
              P, G),
        block(slope_block(cmp_frac, {"retention 1.0 -> 0.25": RETENTION_SERIES}, P, G),
              "C2. Sensitivity to fractionation (x = 1 - retention)", P, G),
        block(tetraploid_block(cmp, M, G), "D. Tetraploid subset, MPAllopp vs GRAMPAIter", M, G),
    ]
    for extra in [lab for lab in order(cmp["label"]) if lab not in METHOD_ORDER]:
        for ref in (P, G):
            blocks.append(block(per_config(cmp, BASELINE, extra, ref, strict=False),
                                f"E. {extra} vs {name(ref)}, baseline configurations", extra, ref))
    res = pd.concat(blocks, ignore_index=True)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "simulation_tests.txt").write_text(report(res) + "\n")
    return [write_csv(res, out_dir, "simulation_tests.csv"), out_dir / "simulation_tests.txt",
            write_csv(accuracy(inv, cmp, BASELINE + frac[1:]), out_dir, "accuracy.csv")]
