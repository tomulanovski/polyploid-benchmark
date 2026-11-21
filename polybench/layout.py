from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Dataset:
    kind: str  # "simulated" or "empirical"
    config: str | None
    network: str  # model network name, or the dataset name for empirical data
    replicate: int | None
    path: Path
    root: Path  # data root

    @property
    def id(self) -> str:
        if self.kind == "empirical":
            return f"empirical/{self.network}"
        return f"simulated/{self.config}/{self.network}/rep{self.replicate}"


def _parse(data_root: Path, path: Path) -> Dataset:
    parts = path.relative_to(data_root).parts
    if parts[0] == "empirical" and len(parts) == 2:
        return Dataset("empirical", None, parts[1], None, path, data_root)
    if parts[0] == "simulated" and len(parts) == 4 and parts[3].startswith("rep") and parts[3][3:].isdigit():
        return Dataset("simulated", parts[1], parts[2], int(parts[3][3:]), path, data_root)
    raise ValueError(f"{path} is not <data_root>/simulated/<config>/<network>/rep<k> or <data_root>/empirical/<dataset>")


def find_datasets(data_root, pattern: str) -> list[Dataset]:
    data_root = Path(data_root)
    paths = sorted(p for p in data_root.glob(pattern) if p.is_dir())
    if not paths:
        raise FileNotFoundError(f"no datasets match '{pattern}' under {data_root}")
    # a missing gene_trees.nwk is reported by validate_dataset, not here
    return sorted((_parse(data_root, p) for p in paths), key=lambda d: d.id)


def result_dir(results_root, label: str, dataset: Dataset) -> Path:
    return Path(results_root) / label / dataset.id


def true_network_path(data_root, dataset: Dataset) -> Path | None:
    if dataset.kind == "empirical":
        return None
    return Path(data_root) / "model_networks" / f"{dataset.network}.tre"


def read_taxon_map(alignments_dir) -> tuple[dict[str, str], dict[str, tuple[str, str]]]:
    # alignments/taxon_map.tsv: seq_id<TAB>taxon[<TAB>individual:genome]. The optional third column fixes the
    # AlloppNET homeolog pairing of a tetraploid's sequences (genome A or B of one individual).
    map_file = Path(alignments_dir) / "taxon_map.tsv"
    taxon_of: dict[str, str] = {}
    homeolog_of: dict[str, tuple[str, str]] = {}
    if not map_file.is_file():
        return taxon_of, homeolog_of
    for n, line in enumerate(map_file.read_text().splitlines(), 1):
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) not in (2, 3) or not all(c.strip() for c in cols):
            raise ValueError(f"{map_file} line {n}: expected 'seq_id<TAB>taxon' or "
                             "'seq_id<TAB>taxon<TAB>individual:genome'")
        if taxon_of.setdefault(cols[0], cols[1]) != cols[1]:
            raise ValueError(f"{map_file} line {n}: '{cols[0]}' mapped to two taxa")
        if len(cols) == 3:
            ind, sep, genome = cols[2].rpartition(":")
            if not sep or not ind or genome not in ("A", "B"):
                raise ValueError(f"{map_file} line {n}: homeolog '{cols[2]}' is not 'individual:A' or 'individual:B'")
            homeolog_of[cols[0]] = (ind, genome)
    return taxon_of, homeolog_of


def sequence_taxa(alignments_dir) -> dict[str, dict[str, str]]:
    # file name -> {sequence id: taxon}; ids missing from taxon_map.tsv (or no map) are their own taxon
    alignments_dir = Path(alignments_dir)
    taxon_of, _ = read_taxon_map(alignments_dir)
    files = sorted(alignments_dir.glob("*.fasta"))
    if not files:
        raise ValueError(f"{alignments_dir}: no .fasta files")
    out, problems = {}, []
    for f in files:
        ids = []
        for n, line in enumerate(f.read_text().splitlines(), 1):
            if line.startswith(">"):
                fields = line[1:].split()
                if fields:
                    ids.append(fields[0])
                else:
                    problems.append(f"{f.name} line {n}: empty sequence id")
            elif line.strip() and not ids:
                problems.append(f"{f.name} line {n}: sequence data before the first '>' header")
        if not ids:
            problems.append(f"{f.name}: no sequences")
        dups = sorted({i for i in ids if ids.count(i) > 1})
        if dups:
            problems.append(f"{f.name}: duplicate sequence ids {dups}")
        out[f.name] = {i: taxon_of.get(i, i) for i in ids}
    if problems:
        raise ValueError("; ".join(problems))
    return out
