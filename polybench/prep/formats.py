import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from ete3 import Tree

from polybench.layout import read_taxon_map, sequence_taxa
from polybench.prep.ploidy import read_ploidy_tsv
from polybench.prep.trees import clean_gene_trees, make_unique_labels, map_leaves, relabel_leaves


def _numbered(newick: str, always: bool, taxon_labels: dict[str, set]) -> str:
    # copy numbers in leaf order; with always=False only repeated taxa are numbered
    tree = Tree(newick, format=1)  # format 1 also reads support values on internal nodes
    totals = Counter(leaf.name for leaf in tree)
    seen: Counter = Counter()
    for leaf in tree:
        taxon = leaf.name
        seen[taxon] += 1
        if always or totals[taxon] > 1:
            leaf.name = f"{taxon}_{seen[taxon]}"
        taxon_labels[taxon].add(leaf.name)
    return tree.write(format=9)


def to_grampa(gene_trees: list[str], mapping: dict[str, str]) -> list[str]:
    out = []
    for newick in gene_trees:
        seen: Counter = Counter()

        def number(name):
            seen[name] += 1
            return f"{seen[name]}_{name}"

        out.append(map_leaves(relabel_leaves(newick, mapping), number))
    return out


def clean_grampa_tree(newick: str) -> str:
    # drops internal labels such as <5> or <P*>, then the * and + copy marks on leaves
    newick = re.sub(r"<[^>]*>", "", newick)
    return re.sub(r"(?<=[A-Za-z0-9_])[*+]+", "", newick)


def write_grampa_inputs(dataset_dir, inputs_dir) -> None:
    # grampa_trees.tre (<copy>_<species> labels), species_tree.tre and taxa_map.tsv (original<TAB>grampa label,
    # one line per taxon) after the substring fix
    dataset_dir, inputs_dir = Path(dataset_dir), Path(inputs_dir)
    lines = [l for l in (dataset_dir / "gene_trees.nwk").read_text().splitlines() if l.strip()]
    trees = clean_gene_trees(lines)
    taxa = sorted({leaf.name for t in trees for leaf in Tree(t, format=1)})
    mapping = make_unique_labels(taxa)
    species = relabel_leaves((dataset_dir / "species_tree.nwk").read_text().strip(), mapping)
    inputs_dir.mkdir(parents=True, exist_ok=True)
    (inputs_dir / "grampa_trees.tre").write_text("".join(t + "\n" for t in to_grampa(trees, mapping)))
    (inputs_dir / "species_tree.tre").write_text(species + "\n")
    (inputs_dir / "taxa_map.tsv").write_text("".join(f"{t}\t{mapping.get(t, t)}\n" for t in taxa))


def to_mpsugar_nexus(gene_trees: list[str]) -> tuple[str, str]:
    labels: dict[str, set] = defaultdict(set)
    trees = [_numbered(t, True, labels) for t in gene_trees]
    lines = ["#NEXUS", "BEGIN TREES;"]
    lines += [f"Tree geneTree{i} = {t}" for i, t in enumerate(trees, 1)]
    lines.append("END;")
    # copies in numeric order (A_2 before A_10)
    taxon_map = {taxon: sorted(found, key=lambda l: int(l.rsplit("_", 1)[1]))
                 for taxon, found in sorted(labels.items())}
    return "\n".join(lines) + "\n", json.dumps(taxon_map, indent=4, sort_keys=True)


PHYLONET_COMMAND = "InferNetwork_MP_Allopp (all) 1 -x 1 -pl 12 -di"


def to_phylonet_nexus(gene_trees: list[str], command: str = PHYLONET_COMMAND) -> str:
    # the taxa map (-a) is built from the copy labels and appended to the command
    labels: dict[str, set] = defaultdict(set)
    trees = [_numbered(t, False, labels) for t in gene_trees]
    taxa_map = "; ".join(f"{t}:{','.join(sorted(found))}" for t, found in sorted(labels.items()))
    lines = ["#NEXUS", "BEGIN TREES;"]
    lines += [f"Tree geneTree{i} = {t}" for i, t in enumerate(trees, 1)]
    lines += ["END;", "", "BEGIN PHYLONET;", f"{command.strip().rstrip(';')} -a <{taxa_map}>;", "END;"]
    return "\n".join(lines) + "\n"


def strip_copy_suffixes(newick: str) -> str:
    newick = re.sub(r"\[[^\]]*\]", "", newick)
    tree = Tree(newick, format=5 if ":" in newick else 9)
    for leaf in tree:
        leaf.name = re.sub(r"_\d+$", "", leaf.name)
    return tree.write(format=9)


def strip_alloppnet_suffixes(newick: str, ploidy: dict[str, int]) -> str:
    # AlloppNET appends 0/1 (no separator) to tetraploid species for homeologs
    polyploids = {t for t, p in ploidy.items() if p == 4}
    newick = re.sub(r"\[[^\]]*\]", "", newick)
    tree = Tree(newick, format=5 if ":" in newick else 9)
    for leaf in tree:
        if leaf.name[-1:] in ("0", "1") and leaf.name[:-1] in polyploids:
            leaf.name = leaf.name[:-1]
    return tree.write(format=9)


def read_fasta(path) -> dict[str, str]:
    # sequence id = first word of the header, as in layout.sequence_taxa
    seqs: dict[str, list[str]] = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in seqs:
                raise ValueError(f"duplicate sequence {name} in {path}")
            seqs[name] = []
        elif line.strip():
            seqs[name].append(line.strip())
    return {n: "".join(parts) for n, parts in seqs.items()}


def alignment_to_nexus(seqs: dict[str, str]) -> str:
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError("alignment sequences differ in length")
    rows = "".join(f"    {name}  {seq}\n" for name, seq in seqs.items())
    return (
        "#NEXUS\nBEGIN DATA;\n"
        f"  DIMENSIONS NTAX={len(seqs)} NCHAR={lengths.pop()};\n"
        "  FORMAT DATATYPE=DNA MISSING=N GAP=-;\n"
        f"  MATRIX\n{rows}  ;\nEND;\n"
    )


def homeolog_problems(alignments: dict[str, dict[str, str]], copies: dict[str, int],
                      homeologs: dict[str, tuple[str, str]]) -> list[str]:
    # the third taxon_map.tsv column: only tetraploids, all or none of a taxon's sequences, an individual in one
    # taxon, and per gene file at most one sequence per (individual, genome)
    problems, taxon_of_ind, with_col, without_col = [], {}, set(), set()
    for gene in sorted(alignments):
        seen = set()
        for seq_id, taxon in alignments[gene].items():
            if seq_id not in homeologs:
                without_col.add(taxon)
                continue
            with_col.add(taxon)
            ind, genome = homeologs[seq_id]
            if copies.get(taxon) != 2:
                problems.append(f"{gene}: sequence {seq_id} has a homeolog but taxon {taxon} is not a tetraploid")
            if taxon_of_ind.setdefault(ind, taxon) != taxon:
                problems.append(f"individual {ind} is used for taxa {taxon_of_ind[ind]} and {taxon}")
            if (ind, genome) in seen:
                problems.append(f"{gene}: individual {ind} has two sequences for genome {genome}")
            seen.add((ind, genome))
    for taxon in sorted(with_col & without_col):
        problems.append(f"taxon {taxon}: homeolog column given for some sequences only")
    return sorted(set(problems))


def taxa_table(alignments: dict[str, dict[str, str]], copies: dict[str, int],
               homeologs: dict[str, tuple[str, str]] | None = None) -> str:
    # AlloppNET taxa table, one row per sequence id over all genes. alignments: gene -> {sequence id: taxon}.
    # Diploid (1 copy) sequences are their own individual. Tetraploid (2 copies) sequences take the individual
    # and genome given in homeologs (taxon_map.tsv third column) when present. Otherwise they are paired into
    # individuals with genomes A and B, gene by gene in sorted id order, and an id keeps the individual it got in
    # the first gene it appears in. An individual with one genome gets a placeholder <id>_miss for the other.
    homeologs = homeologs or {}
    problems = homeolog_problems(alignments, copies, homeologs)
    if problems:
        raise ValueError("; ".join(problems))
    individual: dict[str, tuple[str, str, str]] = {}  # id -> (taxon, individual, genome)
    count: Counter = Counter()
    for gene in sorted(alignments):
        by_taxon = defaultdict(list)
        for seq_id, taxon in alignments[gene].items():
            if taxon not in copies:
                raise ValueError(f"{gene}: taxon {taxon} of sequence {seq_id} is not in ploidy.tsv")
            if copies[taxon] == 1:
                individual[seq_id] = (taxon, f"{seq_id}_ind", "A")
            elif seq_id in homeologs:
                individual[seq_id] = (taxon, *homeologs[seq_id])
            elif seq_id not in individual:
                by_taxon[taxon].append(seq_id)
        for taxon, ids in by_taxon.items():
            ids.sort()
            for i in range(0, len(ids), 2):
                count[taxon] += 1
                ind = f"{taxon}_I{count[taxon]}"
                pair = ids[i:i + 2]
                if len(pair) == 1:
                    pair.append(f"{pair[0]}_miss")
                individual[pair[0]] = (taxon, ind, "A")
                individual[pair[1]] = (taxon, ind, "B")
    genomes = defaultdict(dict)  # individual -> {genome: id}, for the given homeologs
    for seq_id in homeologs:
        if seq_id in individual:
            taxon, ind, genome = individual[seq_id]
            genomes[(taxon, ind)][genome] = seq_id
    for (taxon, ind), present in genomes.items():
        if len(present) == 1:
            (genome, seq_id), = present.items()
            individual[f"{seq_id}_miss"] = (taxon, ind, "B" if genome == "A" else "A")
    rows = sorted(individual.items(), key=lambda kv: (kv[1][0], kv[0]))
    return "\n".join(["ID species individual genome"] + [f"{i} {t} {ind} {g}" for i, (t, ind, g) in rows]) + "\n"


def write_alloppnet_inputs(dataset_dir, inputs_dir) -> None:
    # <gene>.nex per alignments/<gene>.fasta and taxa_table.txt; ploidy.tsv copies 1 = diploid, 2 = tetraploid
    dataset_dir, inputs_dir = Path(dataset_dir), Path(inputs_dir)
    copies = read_ploidy_tsv(dataset_dir / "ploidy.tsv")
    if not any(k == 2 for k in copies.values()):
        raise ValueError(f"{dataset_dir}: AlloppNET needs at least one tetraploid (2 copies in ploidy.tsv)")
    seq_taxa = sequence_taxa(dataset_dir / "alignments")
    inputs_dir.mkdir(parents=True, exist_ok=True)
    genes = {}
    for fname, mapping in seq_taxa.items():
        gene = Path(fname).stem
        (inputs_dir / f"{gene}.nex").write_text(alignment_to_nexus(read_fasta(dataset_dir / "alignments" / fname)))
        genes[gene] = mapping
    _, homeologs = read_taxon_map(dataset_dir / "alignments")
    (inputs_dir / "taxa_table.txt").write_text(taxa_table(genes, copies, homeologs))
