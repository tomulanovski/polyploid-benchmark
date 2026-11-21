import matplotlib.pyplot as plt
import pandas as pd

from polybench.paper.common import BASELINE, FAMILIES, load_simulated, name, order, require_configs, save, write_csv


def completion(inv: pd.DataFrame, configs=BASELINE) -> pd.DataFrame:
    # % of datasets in the inventory that completed; a label never run on a dataset ("missing") counts as not
    # completed, as legacy where every method had a row for all 21 networks
    require_configs(inv, configs, "Fig 5")
    inv = inv[inv["config"].isin(configs)]
    pct = (inv["status"] == "completed").groupby([inv["label"], inv["config"]]).mean().unstack(fill_value=0) * 100
    return pct.reindex(index=order(pct.index), columns=configs, fill_value=0).rename_axis(index="label", columns=None)


def make(summary_dir, networks_dir, out_dir, grampa_iter_fallback: bool = True, configs=None) -> list:
    # configs: columns other than the paper's 12 (e.g. the toy example), drawn without the family grouping
    inv, _ = load_simulated(summary_dir, grampa_iter_fallback)
    pct = completion(inv, BASELINE if configs is None else list(configs))

    fig, ax = plt.subplots(figsize=(14, max(5, len(pct) * 0.7)))
    im = ax.imshow(pct.values, cmap="RdYlGn", aspect="auto", vmin=0, vmax=100)
    for i in range(pct.shape[0]):
        for j in range(pct.shape[1]):
            val = pct.values[i, j]
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center", fontsize=10,
                    fontweight="bold", color="white" if val < 40 or val > 85 else "black")
    ax.set_xticks(range(pct.shape[1]))
    ax.set_yticks(range(len(pct)))
    ax.set_yticklabels([name(m, plot=True) for m in pct.index], fontsize=11)
    if configs is None:
        ax.set_xticklabels([lvl for fam in FAMILIES.values() for lvl in fam], fontsize=10)
        for k, fam in enumerate(FAMILIES):
            start = 3 * k
            if k < len(FAMILIES) - 1:
                ax.axvline(x=start + 2.5, color="black", linewidth=2)
            ax.text(start + 1, 1.02, fam, ha="center", va="bottom", fontsize=10, fontweight="bold",
                    transform=ax.get_xaxis_transform())
    else:
        ax.set_xticklabels(pct.columns, fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02).set_label("Completion Rate (%)", fontsize=12)
    fig.tight_layout()

    table = pct.reset_index()
    table.insert(1, "method", [name(m) for m in pct.index])
    return save(fig, out_dir, "fig5") + [write_csv(table, out_dir, "fig5_completion.csv")]
