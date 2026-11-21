import contextlib
import io
from pathlib import Path

import pandas as pd

from polybench.measures.reticulate_tree import ReticulateTree
from polybench.paper.common import write_csv

RELAXED_THRESHOLD = 0.2  # Polyphest folding threshold for H_Relaxed


def stats_row(path: Path) -> dict:
    # legacy run_reticulation_stats.analyze_mul_trees, one model network
    newick = path.read_text(encoding="utf-8").strip()
    strict = ReticulateTree(newick, is_multree=True)
    counts = strict.measure(printout=False)["leaf_counts"]
    polyploids = sorted(n for n, c in counts.items() if c > 1)
    relaxed = ReticulateTree(newick, is_multree=True, threshold=RELAXED_THRESHOLD, normalize=True)
    auto, h_strict = strict.get_auto_event_count(), strict.get_allo_event_count()
    h_relaxed = relaxed.measure(printout=False)["reticulation_count"]
    return {"Filename": path.name, "Num_Species": len(counts), "Num_Polyploids": len(polyploids),
            "Max_Copies": max(counts.values()), "H_Strict": h_strict, "H_Relaxed": h_relaxed,
            "H_Diff": h_relaxed - h_strict, "Num_Autopolyploidization_Events": auto, "Total_WGD": auto + h_strict,
            "Polyploid_Names": ", ".join(polyploids)}


def network_stats(networks_dir) -> pd.DataFrame:
    paths = sorted(Path(networks_dir).glob("*.tre"))
    if not paths:
        raise FileNotFoundError(f"no .tre model networks in {networks_dir}")
    rows = []
    for p in paths:
        try:
            with contextlib.redirect_stdout(io.StringIO()):  # the measures print a note on Polyphest folding
                rows.append(stats_row(p))
        except Exception as e:
            raise RuntimeError(f"{p}: {e}") from e
    return pd.DataFrame(rows)


def make(summary_dir, networks_dir, out_dir) -> list:
    return [write_csv(network_stats(networks_dir), out_dir, "table_s2.csv")]
