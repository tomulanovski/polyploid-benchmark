from collections import Counter, defaultdict
from pathlib import Path

from ete3 import Tree


def _tip_counts(newick: str) -> Counter:
    return Counter(leaf.name for leaf in Tree(newick, format=1).iter_leaves())


def copy_number_distribution(gene_trees: list[str], taxon: str) -> dict[int, int]:
    dist = Counter()
    for newick in gene_trees:
        n = _tip_counts(newick)[taxon]
        if n:
            dist[n] += 1
    return dict(dist)


def smooth_copy_number(counts: dict[int, int], width: int = 2) -> int:
    ordered = sorted(counts.items())
    if len(ordered) <= 2:
        # most frequent, ties go to the smaller copy number
        return max(ordered, key=lambda kv: (kv[1], -kv[0]))[0]
    best, best_score = None, 0.0
    for center, _ in ordered:
        score = sum(
            f * (1 - abs(center - c) / (width + 1))
            for c, f in ordered
            if abs(center - c) <= width
        )
        if score > best_score:  # strict: ties go to the smaller copy number
            best, best_score = center, score
    return best


def infer_ploidy(gene_trees: list[str]) -> dict[str, int]:
    dists = defaultdict(Counter)
    for newick in gene_trees:
        for taxon, n in _tip_counts(newick).items():
            dists[taxon][n] += 1
    return {t: smooth_copy_number(dict(d)) for t, d in sorted(dists.items())}


def write_ploidy_tsv(ploidy: dict[str, int], path) -> None:
    lines = [f"{t}\t{k}\n" for t, k in sorted(ploidy.items())]
    Path(path).write_text("".join(lines))


def read_ploidy_tsv(path) -> dict[str, int]:
    ploidy = {}
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) != 2:
            raise ValueError(f"{path} line {n}: expected 2 tab-separated columns, got {len(cols)}")
        taxon, copies = cols
        if not copies.strip().isdigit() or int(copies) < 1:
            raise ValueError(f"{path} line {n}: copies must be an integer >= 1, got '{copies}'")
        ploidy[taxon] = int(copies)
    return ploidy


def write_polyphest_multiset(ploidy: dict[str, int], path) -> None:
    lines = [f"{t}\n" * k for t, k in sorted(ploidy.items())]
    Path(path).write_text("".join(lines))
