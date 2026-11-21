import matplotlib.pyplot as plt
import pandas as pd

from polybench.paper.common import (DUP_PANELS, FRACTIONATION, FRACTIONATION_EXCLUDE, ILS_LEVELS, TRIO, colors,
                                    load_simulated, mean_sem, name, per_network, require_configs, require_labels,
                                    require_metrics, save, success, write_csv)

METHODS = ["polyphest", "grampa_iter"]


def check(cmp, configs, methods, what):
    require_configs(cmp, configs, what)
    require_labels(cmp, methods, what)
    require_metrics(cmp, TRIO, what)


def point(cmp, metric, label, config) -> dict:
    vals = per_network(success(cmp, metric).query("label == @label and config == @config"))
    mean, sem = mean_sem(vals)
    return {"label": label, "config": config, "metric": metric, "mean": mean, "sem": sem, "n_networks": len(vals)}


def fig6(cmp, methods, out_dir) -> list:
    # rows: descendants, sisters, mu; columns: dup/loss none..high; x: ILS level (legacy 06b_degradation_lines)
    check(cmp, [c for p in DUP_PANELS.values() for c in p.values()], methods, "Fig 6")
    col = colors(methods)
    fig, axes = plt.subplots(len(TRIO), len(DUP_PANELS), figsize=(5 * len(DUP_PANELS), 4 * len(TRIO)),
                             squeeze=False, sharey="row")
    rows = []
    for r, (metric, ylabel) in enumerate(TRIO.items()):
        for c, (panel, configs) in enumerate(DUP_PANELS.items()):
            ax = axes[r, c]
            for label in methods:
                pts = [{**point(cmp, metric, label, configs[lvl]), "panel": panel, "level": lvl} for lvl in ILS_LEVELS]
                rows += pts
                ax.errorbar(range(len(ILS_LEVELS)), [p["mean"] for p in pts], yerr=[p["sem"] for p in pts],
                            marker="o", markersize=8, capsize=4, color=col[label], linewidth=2.5,
                            label=name(label, plot=True) if r == 0 and c == 0 else None)
            ax.set_xticks(range(len(ILS_LEVELS)))
            ax.set_xticklabels(ILS_LEVELS)
            ax.set_xlim(-0.3, len(ILS_LEVELS) - 0.7)
            ax.set_ylim(0, 1)
            ax.grid(True, alpha=0.25, linestyle="--")
            if r == 0:
                ax.set_title(panel, fontsize=13, fontweight="bold", pad=10)
            if r == len(TRIO) - 1:
                ax.set_xlabel("ILS Level", fontsize=11, fontweight="bold")
            if c == 0:
                ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(methods), fontsize=12, framealpha=0.9,
               bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return save(fig, out_dir, "fig6") + [write_csv(pd.DataFrame(rows), out_dir, "fig6.csv")]


def fig7(cmp, methods, out_dir) -> list:
    # retention 1.00 -> 0.25 on medium ILS, medium dup/loss; excluded networks dropped from every level
    configs = [c for _, c in FRACTIONATION]
    check(cmp, configs, methods, "Fig 7")
    excluded = cmp["network"].isin(FRACTIONATION_EXCLUDE)
    if (excluded & cmp["config"].isin(configs[1:]) & (cmp["status"] == "SUCCESS")).any():
        raise ValueError(f"{FRACTIONATION_EXCLUDE} have fractionation results but are listed as never simulated")
    cmp = cmp[~excluded]
    col = colors(methods)
    retention = [r for r, _ in FRACTIONATION]
    fig, axes = plt.subplots(1, len(TRIO), figsize=(5.2 * len(TRIO), 4.6), squeeze=False)
    rows = []
    for ax, (metric, ylabel) in zip(axes[0], TRIO.items()):
        for label in methods:
            pts = [{**point(cmp, metric, label, config), "retention": r} for r, config in FRACTIONATION]
            rows += pts
            ax.errorbar(retention, [p["mean"] for p in pts], yerr=[p["sem"] for p in pts], marker="o",
                        markersize=8, capsize=4, color=col[label], linewidth=2.5, label=name(label, plot=True))
        ax.grid(True, alpha=0.25, linestyle="--")
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
        ax.set_xlabel("Post-WGD retention rate", fontsize=11, fontweight="bold")
        ax.set_ylim(0, 1)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(methods), fontsize=12, framealpha=0.9,
               bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return save(fig, out_dir, "fig7") + [write_csv(pd.DataFrame(rows), out_dir, "fig7.csv")]


def make(summary_dir, networks_dir, out_dir, methods=METHODS, grampa_iter_fallback: bool = True) -> list:
    _, cmp = load_simulated(summary_dir, grampa_iter_fallback)
    return fig6(cmp, list(methods), out_dir) + fig7(cmp, list(methods), out_dir)
