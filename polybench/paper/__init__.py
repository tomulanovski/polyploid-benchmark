from polybench.paper import fig3_table3, fig5, fig6_fig7, section33, stats, table2, table4, table_s2

MODULES = [fig3_table3, table2, fig5, fig6_fig7, table4, stats, section33, table_s2]
SIMULATED = {fig5, fig6_fig7, table4, stats}  # modules whose make takes grampa_iter_fallback
NAMES = {m.__name__.rsplit(".", 1)[1]: m for m in MODULES}


def make_all(summary_dir, networks_dir, out_dir, grampa_iter_fallback: bool = True, only=None,
             fig5_configs=None) -> list:
    # only: module names from NAMES to make (default all); fig5_configs: Fig 5 columns other than the paper's 12
    unknown = sorted(set(only or []) - set(NAMES))
    if unknown:
        raise ValueError(f"unknown paper outputs {unknown}, expected some of {list(NAMES)}")
    paths = []
    for name, m in NAMES.items():
        if only is not None and name not in only:
            continue
        extra = {"grampa_iter_fallback": grampa_iter_fallback} if m in SIMULATED else {}
        if m is fig5 and fig5_configs is not None:
            extra["configs"] = fig5_configs
        paths += m.make(summary_dir, networks_dir, out_dir, **extra)
    return paths
