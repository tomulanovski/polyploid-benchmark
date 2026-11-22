#!/usr/bin/env python3
# Packs the legacy gene2net folders into the Dryad layout of the README (data/, results/<label>/, distributions/).
# Run on the cluster from the repo root with the polybench env active, see dryad/README_template.md.
import argparse
import json
import math
import re
import shutil
import sys
import tarfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import yaml
from ete3 import Tree

from polybench.evaluate.compare import build_inventory
from polybench.layout import find_datasets, read_taxon_map, sequence_taxa
from polybench.prep.formats import clean_grampa_tree, taxa_table
from polybench.prep.ploidy import infer_ploidy
from polybench.prep.trees import split_newick
from polybench.runner import load_method
from polybench.simulate.fractionation import species_tree_leaf
from polybench.simulate.simphy import CONFIGS, natural_key, restore_simphy_labels
from polybench.validate import validate_or_raise

HERE = Path(__file__).resolve().parent
POLYPHEST_EXAMPLE = HERE.parent / "methods" / "polyphest" / "examples" / "p50.yaml"

# results label -> (folder under <net>/results/<config>/, result file written by postprocess_results.py)
SIM_RESULTS = {
    "padre": ("padre", "padre_result.tre"),
    "polyphest_p50": ("polyphest_p50", "polyphest_result.tre"),
    "polyphest_p70": ("polyphest_p70", "polyphest_result.tre"),
    "polyphest_p90": ("polyphest_p90", "polyphest_result.tre"),
    "grampa": ("grampa", "grampa_result.tre"),
    "grampa_iter": ("grandma_split", "grandma_result.tre"),
    "mpsugar": ("mpsugar", "mpsugar_result.tre"),
    "alloppnet": ("alloppnet", "alloppnet_result.tre"),
}
POLYPHEST_LABEL = re.compile(r"polyphest_p(\d+)$")
LABEL_METHOD = {"padre": "padre", "grampa": "grampa", "grampa_iter": "grampa_iter", "mpsugar": "mpsugar",
                "alloppnet": "alloppnet", "mpallopp_phylonet": "mpallopp_phylonet"}
PICKLES = ["substitution_rates.pkl", "gtr_parameters_all.pkl"]


@dataclass
class SimRep:
    config: str
    network: str
    replicate: int
    rep_dir: Path
    gene_files: list
    species_file: Path
    sim_config_file: Path
    alloppnet_input: Path | None


@dataclass
class Result:
    label: str
    dataset_id: str
    status: str
    source: Path | None  # network file to copy, None for failed runs


@dataclass
class Plan:
    networks: list = field(default_factory=list)  # model network files
    sim: list = field(default_factory=list)
    emp: dict = field(default_factory=dict)  # dataset -> {kind: source path}
    results: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    copied_grampa_iter: list = field(default_factory=list)


def grampa_iter_copied(run_dir: Path) -> bool:
    # GRANDMA's own output final_multree.tre is absent: the grandma_result.tre next to it is the first-iteration
    # GRAMPA result copied in by hand, so the run is packed as failed and the pipeline's fallback rescores it
    return (run_dir / "grandma_result.tre").is_file() and not nonempty(run_dir / "final_multree.tre")


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0 and path.read_text(encoding="utf-8").strip() != ""


def topology(newick: str) -> str:
    # leaf labels and topology only, children sorted, after the GRAMPA label cleanup
    def walk(node):
        if node.is_leaf():
            return node.name
        return "(" + ",".join(sorted(walk(c) for c in node.children)) + ")"
    return walk(Tree(clean_grampa_tree(newick.strip()), format=1))


def check_copies(plan: "Plan") -> None:
    # every run packed as a copied GRAMPA result should hold the same tree as that replicate's grampa_result.tre
    match, mismatch, no_grampa = [], [], []
    for dataset_id, run_dir in plan.copied_grampa_iter:
        grampa = run_dir.parents[1] / "grampa" / run_dir.name / "grampa_result.tre"
        if not nonempty(grampa):
            no_grampa.append(dataset_id)
        elif topology((run_dir / "grandma_result.tre").read_text()) == topology(grampa.read_text()):
            match.append(dataset_id)
        else:
            mismatch.append(dataset_id)
    print(f"copied GRAMPAIter runs vs grampa_result.tre: {len(match)} match, {len(mismatch)} mismatch, "
          f"{len(no_grampa)} without a GRAMPA result")
    for what, ids in (("mismatch", mismatch), ("no GRAMPA result", no_grampa)):
        for i in ids:
            print(f"  {what}: {i}")


def sim_result(net_dir: Path, legacy: str, label: str, k: int, plan: Plan, dataset_id: str) -> None:
    # legacy: the config's folder name in gene2net
    folder, fname = SIM_RESULTS[label]
    run_dir = net_dir / "results" / legacy / folder / f"replicate_{k}"
    # PADRE wrote its output into processed/, postprocess_results.py only made the results folder when it found one
    ran = (net_dir / "processed" / legacy / "padre_input" / f"replicate_{k}").is_dir() if label == "padre" \
        else run_dir.is_dir()
    if not ran:
        return
    result = run_dir / fname
    if label == "grampa_iter":
        if grampa_iter_copied(run_dir):
            plan.copied_grampa_iter.append((dataset_id, run_dir))
            plan.results.append(Result(label, dataset_id, "failed", None))
            return
        if nonempty(run_dir / "final_multree.tre") and not nonempty(result):
            raise FileNotFoundError(f"{run_dir}: final_multree.tre exists but {fname} is missing or empty, "
                                    f"run postprocess_results.py {legacy} --methods grandma_split first")
    ok = nonempty(result)
    plan.results.append(Result(label, dataset_id, "completed" if ok else "failed", result if ok else None))


def species_file_for(rep_dir: Path, config: str, net_dir: Path, k: int, legacy_names: dict) -> Path:
    found = sorted(rep_dir.rglob("s_tree.trees"), key=natural_key)
    if found:
        return found[0]
    base = CONFIGS[config].get("base")
    if base:  # fractionated configs hold only the pruned gene trees, the species tree is the base config's
        base_dir = net_dir / "data" / legacy_names.get(base, base) / f"replicate_{k}"
        found = sorted(base_dir.rglob("s_tree.trees"), key=natural_key)
        if found:
            return found[0]
    raise FileNotFoundError(f"no s_tree.trees under {rep_dir}" + (f" or the {base} replicate" if base else ""))


def plan_simulated(root: Path, configs, n_trees: int, plan: Plan, legacy_names: dict | None = None) -> None:
    # legacy_names: config -> its folder name in gene2net, when it differs
    sims = root / "simulations" / "simulations"
    networks = sorted(p for p in (root / "simulations" / "networks").glob("*.tre"))
    if not networks:
        raise FileNotFoundError(f"no model networks in {root / 'simulations' / 'networks'}")
    plan.networks = networks
    legacy_names = legacy_names or {}
    for config in configs:
        if config not in CONFIGS:
            raise ValueError(f"unknown config {config}, expected one of {list(CONFIGS)}")
        legacy = legacy_names.get(config, config)
        base = CONFIGS[config].get("base", config)
        cfg_dir_name = legacy_names.get(base, base)
        found = 0
        for net in networks:
            net_dir = sims / net.stem
            if not net_dir.is_dir():
                raise FileNotFoundError(f"{net_dir} not found for model network {net.name}")
            data = net_dir / "data" / legacy
            if not data.is_dir():
                plan.notes.append(f"no data: {data}")
                continue
            sim_config_file = net_dir / "data" / cfg_dir_name / "simulation_config.txt"
            reps = sorted((p for p in data.glob("replicate_*") if p.name[10:].isdigit()),
                          key=lambda p: int(p.name[10:]))
            if not reps:
                raise FileNotFoundError(f"{data}: no replicate_<k> folders")
            for rep_dir in reps:
                k = int(rep_dir.name[10:])
                # legacy prep order: shell glob over batch_*/ then g_*.trees, i.e. lexicographic per folder
                genes = sorted(rep_dir.rglob("g_trees*.trees"), key=lambda f: f.relative_to(rep_dir).parts)
                if len(genes) != n_trees:
                    raise ValueError(f"{rep_dir}: expected {n_trees} g_trees*.trees files, found {len(genes)}")
                alloppnet_input = net_dir / "processed" / legacy / "alloppnet_input" / f"replicate_{k}"
                species = species_file_for(rep_dir, config, net_dir, k, legacy_names)
                plan.sim.append(SimRep(config, net.stem, k, rep_dir, genes, species, sim_config_file,
                                       alloppnet_input if alloppnet_input.is_dir() else None))
                dataset_id = f"simulated/{config}/{net.stem}/rep{k}"
                for label in SIM_RESULTS:
                    sim_result(net_dir, legacy, label, k, plan, dataset_id)
            found += len(reps)
            if len(reps) != 5:
                plan.notes.append(f"{len(reps)} replicates (not 5): {data}")
        if not found:
            raise FileNotFoundError(f"config {config}: no replicates under {sims}/<network>/data/{legacy}"
                                    + ("" if legacy != config else ", map a differently named folder with --config-map"))


def plan_empirical(root: Path, sources_file: Path, plan: Plan) -> None:
    src = yaml.safe_load(sources_file.read_text())
    fill = [f"{k}" for k in src["results"] if "FILL" in k] + [f"results.{k}" for k, v in src["results"].items()
                                                              if "FILL" in str(v)]
    for ds, entry in src["datasets"].items():
        fill += [f"datasets.{ds}.{k}" for k, v in (entry or {}).items() if "FILL" in str(v)]
    if fill:
        raise ValueError(f"{sources_file}: fill in the FILL entries first: {fill}")
    for label in src["results"]:
        if label not in LABEL_METHOD and not POLYPHEST_LABEL.match(label):
            raise ValueError(f"{sources_file}: unknown results label {label}")
    papers = root / "papers"
    for ds, entry in src["datasets"].items():
        d = papers / ds
        if not d.is_dir():
            raise FileNotFoundError(f"{d} not found (dataset {ds} in {sources_file})")
        if "gene_trees" not in (entry or {}):
            raise ValueError(f"{sources_file}: dataset {ds} has no gene_trees entry")
        files = {}
        for kind, rel in entry.items():
            if kind not in ("gene_trees", "ploidy", "species_tree", "published_network", "alignments"):
                raise ValueError(f"{sources_file}: dataset {ds}: unknown key {kind}")
            if rel is None:
                continue
            p = d / rel
            if not (p.is_dir() if kind == "alignments" else p.is_file()):
                raise FileNotFoundError(f"{p} not found ({ds}.{kind} in {sources_file})")
            files[kind] = p
        plan.emp[ds] = files
        for label, rel in src["results"].items():
            result = d / "networks" / rel
            if not result.parent.is_dir():
                continue  # method never run on this dataset
            ok = nonempty(result)
            plan.results.append(Result(label, f"empirical/{ds}", "completed" if ok else "failed",
                                       result if ok else None))


def method_info(label: str) -> dict:
    m = POLYPHEST_LABEL.match(label)
    spec = load_method("polyphest" if m else LABEL_METHOD[label])
    params = dict(spec.params)
    if m:
        params.update(yaml.safe_load(POLYPHEST_EXAMPLE.read_text()), percentile=int(m[1]))
    return {"method": spec.name, "version": spec.version, "params": params}


def read_sim_config(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found")
    text = path.read_text()

    def grab(pattern):
        m = re.search(pattern, text, re.M)
        if not m:
            raise ValueError(f"{path}: no line matching '{pattern}'")
        return m[1].strip()
    return {"ne": int(float(grab(r"^Effective Population Size \(Ne\):\s*(\S+)"))),
            "dup": float(grab(r"^Duplication rate:.*\((\S+) per gen\)")),
            "loss": float(grab(r"^Loss rate:.*\((\S+) per gen\)")),
            "height": grab(r"^Tree Height:\s*(\S+)"),
            "subst_rate": grab(r"^Substitution rate \(sampled\):\s*(\S+)"),
            "seed": grab(r"^Base seed:\s*(\S+)")}


def simphy_params(r: SimRep, n_trees: int) -> str:
    c = CONFIGS[r.config]
    legacy = read_sim_config(r.sim_config_file)
    if legacy["height"] != "10M":
        raise ValueError(f"{r.sim_config_file}: tree height {legacy['height']}, expected 10M")
    for key in ("ne", "dup", "loss"):
        if not math.isclose(legacy[key], c[key], rel_tol=1e-6):
            raise ValueError(f"{r.sim_config_file}: {key} {legacy[key]} does not match {r.config} ({c[key]})")
    info = {"ne": c["ne"], "dup": c["dup"], "loss": c["loss"], "height": 10_000_000,
            "subst_rate": f"f:{legacy['subst_rate']}"}
    if "base" in c:
        summary = json.loads((r.rep_dir / "diploidization_summary.json").read_text())
        if summary["dist"] != "fixed" or not math.isclose(summary["dist_params"]["value"], c["retention"]):
            raise ValueError(f"{r.rep_dir}: diploidization {summary['dist']} {summary['dist_params']} does not match "
                             f"{r.config} (fixed {c['retention']})")
        info["base_config"] = c["base"]
        info["retention_seed"] = summary["seed"]
    info.update(config=r.config, network=r.network, replicate=r.replicate, n_trees=n_trees,
                retention=c.get("retention"), legacy_base_seed=legacy["seed"], source="gene2net")
    return "".join(f"{k}={v}\n" for k, v in info.items())


def species_mapping(stree: str, model_net: Path) -> dict:
    # species tree leaf -> taxon (legacy labels: taxon, or taxon_<copy> for the copies of a polyploid)
    leaves = Tree(stree, format=1).get_leaf_names()
    mapping = {leaf: leaf.split("_", 1)[0] for leaf in leaves}
    truth = Counter(Tree(model_net.read_text().strip(), format=1).get_leaf_names())
    if Counter(mapping.values()) != truth:
        raise ValueError(f"SimPhy species tree leaves {sorted(leaves)} do not match model network {model_net.name}")
    return mapping


def read_phylip(path: Path) -> dict[str, str]:
    seqs = {}
    for line in path.read_text().splitlines()[1:]:
        if line.strip():
            name, seq = line.split(None, 1)
            if name in seqs:
                raise ValueError(f"{path}: duplicate sequence {name}")
            seqs[name] = seq.replace(" ", "")
    if not seqs:
        raise ValueError(f"{path}: no sequences")
    return seqs


def legacy_homeologs(alignment_ids: set[str], alloppnet_input: Path) -> tuple[dict, dict[str, int]]:
    # tetraploid sequence id -> (individual, genome) from the legacy taxa_table.txt, and the legacy copy numbers
    ploidy = json.loads((alloppnet_input / "ploidy_level.json").read_text())
    copies = {t: 1 if p == 2 else 2 for t, p in ploidy.items()}
    homeologs = {}
    for line in (alloppnet_input / "taxa_table.txt").read_text().strip().splitlines()[1:]:
        seq_id, taxon, ind, genome = line.split()
        if seq_id in alignment_ids and copies[taxon] == 2:
            homeologs[seq_id] = (ind, genome)
    return homeologs, copies


def check_pairing(aln_out: Path, copies: dict[str, int], alloppnet_input: Path) -> None:
    # the AlloppNET prepare rule on the packed alignments must give the legacy taxa_table.txt the paper runs used
    _, homeologs = read_taxon_map(aln_out)
    genes = {Path(f).stem: m for f, m in sequence_taxa(aln_out).items()}
    legacy = alloppnet_input / "taxa_table.txt"
    new = sorted(taxa_table(genes, copies, homeologs).strip().splitlines()[1:])
    old = sorted(legacy.read_text().strip().splitlines()[1:])
    if new != old:
        raise ValueError(f"{legacy}: the packed alignments do not reproduce the paper taxa table "
                         f"(only legacy: {sorted(set(old) - set(new))[:5]}, only new: {sorted(set(new) - set(old))[:5]})")


def write_alignments(r: SimRep, mapping: dict, out: Path) -> dict[str, int]:
    aln_dir = r.rep_dir / "1" / "alignments"
    files = sorted(aln_dir.glob("alignment_*.phy"), key=str)  # legacy glob order
    if not files:
        raise FileNotFoundError(f"{aln_dir}: no alignment_*.phy, but {r.alloppnet_input} says AlloppNET ran")
    out.mkdir(parents=True, exist_ok=True)
    taxon_of = {}
    for i, f in enumerate(files, 1):
        seqs = read_phylip(f)
        (out / f"gene_{i:04d}.fasta").write_text("".join(f">{s}\n{q}\n" for s, q in seqs.items()))
        for s in seqs:
            leaf = species_tree_leaf(s)
            if leaf not in mapping:
                raise ValueError(f"{f}: sequence {s} is not from a species tree leaf")
            taxon_of[s] = mapping[leaf]
    homeologs, copies = legacy_homeologs(set(taxon_of), r.alloppnet_input)
    (out / "taxon_map.tsv").write_text("".join(
        f"{s}\t{t}" + (f"\t{homeologs[s][0]}:{homeologs[s][1]}" if s in homeologs else "") + "\n"
        for s, t in sorted(taxon_of.items())))
    check_pairing(out, copies, r.alloppnet_input)
    return copies


def write_sim_rep(r: SimRep, out_data: Path, n_trees: int) -> None:
    out = out_data / "simulated" / r.config / r.network / f"rep{r.replicate}"
    out.mkdir(parents=True, exist_ok=True)
    stree = r.species_file.read_text().strip()
    mapping = species_mapping(stree, out_data / "model_networks" / f"{r.network}.tre")
    genes = []
    for f in r.gene_files:
        lines = [x for x in f.read_text().splitlines() if x.strip()]
        if len(lines) != 1:
            raise ValueError(f"{f}: expected one tree, found {len(lines)}")
        genes.append(lines[0].strip())
    (out / "gene_trees.nwk").write_text("".join(restore_simphy_labels(g, mapping) + "\n" for g in genes))
    (out / "simphy_gene_trees.nwk").write_text("".join(g + "\n" for g in genes))
    (out / "simphy_species_tree.nwk").write_text(stree + "\n")
    (out / "simphy_params.txt").write_text(simphy_params(r, n_trees))
    if r.alloppnet_input:
        legacy = write_alignments(r, mapping, out / "alignments")
        inferred = infer_ploidy([restore_simphy_labels(g, mapping) for g in genes])
        if any(min(inferred[t], 2) != k for t, k in legacy.items()) or any(inferred[t] > 2 for t in legacy):
            # the rerun infers ploidy from the gene trees; AlloppNET then refuses or treats a taxon differently
            print(f"note: {out}: inferred copies {inferred} differ from the AlloppNET ploidy of the paper run {legacy}")


def split_trees(path: Path) -> list[str]:
    trees = split_newick(path.read_text(encoding="utf-8"))
    for t in trees:
        Tree(t, format=1)  # raises with the parse error
    if not trees:
        raise ValueError(f"{path}: no trees")
    return [" ".join(t.split()) for t in trees]


def read_ploidy_source(path: Path) -> dict[str, int]:
    # taxa_copies.json ({taxon: copies}) or a Polyphest multiset file (one taxon per line, repeated per copy)
    if path.suffix == ".json":
        return {t: int(k) for t, k in json.loads(path.read_text()).items()}
    return dict(Counter(line.strip() for line in path.read_text().splitlines() if line.strip()))


def write_empirical(ds: str, files: dict, out_data: Path) -> None:
    out = out_data / "empirical" / ds
    out.mkdir(parents=True, exist_ok=True)
    (out / "gene_trees.nwk").write_text("".join(t + "\n" for t in split_trees(files["gene_trees"])), encoding="utf-8")
    if "ploidy" in files:
        ploidy = read_ploidy_source(files["ploidy"])
        (out / "ploidy.tsv").write_text("".join(f"{t}\t{k}\n" for t, k in sorted(ploidy.items())), encoding="utf-8")
    for kind, name in (("species_tree", "species_tree.nwk"), ("published_network", "published_network.nwk")):
        if kind in files:
            trees = split_trees(files[kind])
            if len(trees) != 1:
                raise ValueError(f"{files[kind]}: expected one tree, found {len(trees)}")
            (out / name).write_text(trees[0] + "\n", encoding="utf-8")
    if "alignments" in files:
        fastas = sorted(files["alignments"].glob("*.fasta"))
        if not fastas:
            raise FileNotFoundError(f"{files['alignments']}: no .fasta files")
        (out / "alignments").mkdir(exist_ok=True)
        for f in fastas + [p for p in [files["alignments"] / "taxon_map.tsv"] if p.is_file()]:
            shutil.copyfile(f, out / "alignments" / f.name)


def write_result(res: Result, out: Path, root: Path, infos: dict) -> None:
    d = out / "results" / res.label / res.dataset_id
    d.mkdir(parents=True, exist_ok=True)
    if res.source is not None:
        (d / "network.nwk").write_text(res.source.read_text(encoding="utf-8").strip() + "\n", encoding="utf-8")
    status = {"status": res.status, "runtime_s": None, **infos[res.label],
              "source": "gene2net legacy run" + (f": {res.source.relative_to(root).as_posix()}" if res.source else "")}
    (d / "status.json").write_text(json.dumps(status, indent=2))


def summary(plan: Plan) -> None:
    counts = defaultdict(Counter)
    for r in plan.results:
        counts[r.label][r.status] += 1
    for label in sorted(counts):
        print(f"  {label}: " + ", ".join(f"{n} {s}" for s, n in sorted(counts[label].items())))
    print(f"grampa_iter: {len(plan.copied_grampa_iter)} runs had grandma_result.tre without GRANDMA's "
          "final_multree.tre (the copied GRAMPA result), packed as failed:")
    by_net = Counter(i.split("/")[2] for i, _ in plan.copied_grampa_iter)
    for net, n in by_net.most_common():
        print(f"  {net}: {n}")
    for note in plan.notes:
        print(f"note: {note}")


def make_archives(out: Path, configs) -> None:
    arch = out / "archives"
    arch.mkdir(exist_ok=True)
    for config in configs:
        data = out / "data" / "simulated" / config
        if not data.is_dir():
            continue
        targets = {f"simulated_{config}.tar.gz": [data],
                   f"results_{config}.tar.gz": sorted(out.glob(f"results/*/simulated/{config}"))}
        for name, dirs in targets.items():
            with tarfile.open(arch / name, "w:gz") as tar:
                for d in dirs:
                    tar.add(d, arcname=d.relative_to(out).as_posix())
            print(f"{arch / name}: {(arch / name).stat().st_size / 1e6:.1f} MB")


def dir_size(paths) -> int:
    return sum(p.stat().st_size for p in paths)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Pack the legacy gene2net folders into the Dryad layout")
    ap.add_argument("--gene2net-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--configs", nargs="+", default=list(CONFIGS))
    ap.add_argument("--config-map", action="append", default=[], metavar="LEGACY=NEW",
                    help="gene2net folder name of a config, e.g. conf_dup_loss_medium_10M_ne1M_dip075="
                         "conf_dup_loss_medium_10M_ne1M_fix075 (repeatable)")
    ap.add_argument("--empirical-sources", type=Path, default=HERE / "empirical_sources.yaml",
                    help="dataset -> legacy files under papers/<dataset>/ and results label -> file under networks/")
    ap.add_argument("--no-empirical", action="store_true",
                    help="pack only the simulations; empirical folders are copied as run (see the Dryad README)")
    ap.add_argument("--n-trees", type=int, default=1000, help="gene trees per replicate")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and the source sizes, write nothing")
    ap.add_argument("--no-archives", action="store_true", help="do not write the per-config .tar.gz files")
    a = ap.parse_args(argv)
    root, out = a.gene2net_root.resolve(), a.out.resolve()

    legacy_names = {}
    for item in a.config_map:
        old, sep, new = item.partition("=")
        if not sep or not old or new not in CONFIGS:
            raise ValueError(f"--config-map {item}: expected LEGACY=NEW with NEW one of {list(CONFIGS)}")
        legacy_names[new] = old
    plan = Plan()
    plan_simulated(root, a.configs, a.n_trees, plan, legacy_names)
    print(f"model networks: {len(plan.networks)}")
    by_config = Counter(r.config for r in plan.sim)
    for config in a.configs:
        reps = [r for r in plan.sim if r.config == config]
        print(f"simulated {config}: {by_config[config]} replicates, "
              f"{dir_size(f for r in reps for f in r.gene_files) / 1e6:.1f} MB gene trees, "
              f"{sum(r.alloppnet_input is not None for r in reps)} with AlloppNET alignments")
    if not a.no_empirical:
        plan_empirical(root, a.empirical_sources, plan)
    dists = [root / "simulations" / "distributions" / p for p in PICKLES]
    for p in dists:
        if not p.is_file():
            raise FileNotFoundError(f"{p} not found")
    print(f"empirical: {len(plan.emp)} datasets ({', '.join(plan.emp)})")
    print("results:")
    summary(plan)
    if a.dry_run:
        check_copies(plan)
        print("dry run, nothing written")
        return

    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"{out} is not empty")
    data = out / "data"
    (data / "model_networks").mkdir(parents=True)
    for net in plan.networks:
        shutil.copyfile(net, data / "model_networks" / net.name)
    for r in plan.sim:
        write_sim_rep(r, data, a.n_trees)
    for ds, files in plan.emp.items():
        write_empirical(ds, files, data)
    infos = {label: method_info(label) for label in {r.label for r in plan.results}}
    for res in plan.results:
        write_result(res, out, root, infos)
    (out / "distributions").mkdir()
    for p in dists:
        shutil.copyfile(p, out / "distributions" / p.name)

    datasets = [d for kind in ("simulated/*/*/rep*", "empirical/*") if any(data.glob(kind))
                for d in find_datasets(data, kind)]
    validate_or_raise(datasets)
    inv = build_inventory(out / "results", data, sorted(infos))
    print(f"validated {len(datasets)} datasets; inventory status counts:")
    print(inv.groupby(["label", "status"]).size().to_string())
    if not a.no_archives:
        make_archives(out, a.configs)

    print(f"""
Next steps:
1. Write {out / 'config.yaml'} with
     data_root: data
     results_root: results
     tools: {{}}
     executor: {{type: local}}
   and run
     polybench evaluate --config {out / 'config.yaml'} --workers 16 --nrf conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M
2. polybench paper --summary {out / 'results' / 'summary'} --networks {data / 'model_networks'} --out {out / 'supplement'} --only table_s2
   and rename supplement/table_s2.csv to Table_S2.csv. Delete config.yaml afterwards.
3. Upload: README.md, data/model_networks, data/empirical, results/<label>/empirical, results/summary,
   distributions, supplement, code, and archives/ (they hold data/simulated and results/<label>/simulated).""")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, FileExistsError) as e:
        print(f"pack_dryad: error: {e}", file=sys.stderr)
        sys.exit(1)
