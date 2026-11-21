import pickle
import random
import re
import subprocess
from collections import Counter
from pathlib import Path

from ete3 import Tree

from polybench.simulate.fractionation import apply_fractionation, species_tree_leaf
from polybench.simulate.rescale import rescale, to_nexus

# rates are per gene per generation, Ne in individuals
CONFIGS = {
    "conf_ils_low_10M": {"ne": 200_000, "dup": 0, "loss": 0},
    "conf_ils_medium_10M": {"ne": 1_000_000, "dup": 0, "loss": 0},
    "conf_ils_high_10M": {"ne": 2_000_000, "dup": 0, "loss": 0},
    "conf_dup_loss_low_10M": {"ne": 200_000, "dup": 1e-9, "loss": 1e-9},
    "conf_dup_loss_low_10M_ne1M": {"ne": 1_000_000, "dup": 1e-9, "loss": 1e-9},
    "conf_dup_loss_low_10M_ne2M": {"ne": 2_000_000, "dup": 1e-9, "loss": 1e-9},
    "conf_dup_loss_medium_10M": {"ne": 200_000, "dup": 1e-8, "loss": 1e-8},
    "conf_dup_loss_medium_10M_ne1M": {"ne": 1_000_000, "dup": 1e-8, "loss": 1e-8},
    "conf_dup_loss_medium_10M_ne2M": {"ne": 2_000_000, "dup": 1e-8, "loss": 1e-8},
    "conf_dup_loss_high_10M": {"ne": 200_000, "dup": 1e-7, "loss": 1e-7},
    "conf_dup_loss_high_10M_ne1M": {"ne": 1_000_000, "dup": 1e-7, "loss": 1e-7},
    "conf_dup_loss_high_10M_ne2M": {"ne": 2_000_000, "dup": 1e-7, "loss": 1e-7},
    # fractionation of the medium/medium gene trees, no new SimPhy run
    **{f"conf_dup_loss_medium_10M_ne1M_fix{tag}": {"ne": 1_000_000, "dup": 1e-8, "loss": 1e-8, "retention": q,
                                                   "base": "conf_dup_loss_medium_10M_ne1M"}
       for tag, q in (("075", 0.75), ("050", 0.50), ("025", 0.25))},
}


def relabel_for_simphy(multree: str):
    # a species appearing k times becomes s_1..s_k; mapping is relabelled leaf -> original name
    tree = Tree(multree, format=1)
    names = [leaf.name for leaf in tree.iter_leaves()]
    bad = sorted({n for n in names if "_" in n})
    if bad:
        raise ValueError(f"taxon names must not contain underscores (SimPhy label parsing): {bad[:5]}")
    total, seen, mapping = Counter(names), Counter(), {}
    for leaf in tree.iter_leaves():
        orig = leaf.name
        if total[orig] > 1:
            seen[orig] += 1
            leaf.name = f"{orig}_{seen[orig]}"
        mapping[leaf.name] = orig
    return tree.write(format=9), mapping


def rename_leaves(newick: str, rename) -> str:
    tree = Tree(newick, format=1)
    for leaf in tree.iter_leaves():
        leaf.name = rename(leaf.name)
    return tree.write(format=5)


def restore_simphy_labels(newick: str, mapping: dict) -> str:
    def rename(tip):
        leaf = species_tree_leaf(tip)
        if leaf not in mapping:
            raise ValueError(f"gene tree tip '{tip}' ('{leaf}') is not a species tree leaf")
        return mapping[leaf]
    return rename_leaves(newick, rename)


def prepare_species_tree(multree: str, height: float = 10_000_000, seed: int = 0):
    relabelled, mapping = relabel_for_simphy(multree)
    return to_nexus(rescale(relabelled, height, seed)), mapping


def _pool(pkl, key=None):
    with open(pkl, "rb") as f:
        data = pickle.load(f)
    pool = [entry[key] if key else entry for entries in data.values() for entry in entries]
    if not pool:
        raise ValueError(f"{pkl}: no entries to sample from")
    return pool


def sample_subst_rate(subst_rates_pkl, seed: int) -> float:
    return random.Random(seed).choice(_pool(subst_rates_pkl, "rate"))


def _f(x) -> str:
    return f"f:{x}" if isinstance(x, int) else f"f:{x:g}"


def simphy_command(species_tree, out_dir, ne, dup_rate, loss_rate, n_trees: int = 1000, replicates: int = 5, *,
                   seed: int, subst_rates_pkl, rate_seed: int | None = None, simphy_bin: str = "simphy") -> list[str]:
    rate = sample_subst_rate(subst_rates_pkl, seed if rate_seed is None else rate_seed)
    return [str(simphy_bin), "-rs", str(replicates), "-rl", _f(n_trees), "-rg", "1",
            "-SR", str(species_tree), "-lb", _f(dup_rate), "-ld", _f(loss_rate), "-sp", _f(ne),
            "-su", f"f:{rate:.15e}", "-sg", "f:1", "-lt", "f:0", "-lg", "f:0",
            "-cs", str(seed), "-o", str(out_dir), "-v", "1", "-od", "1", "-op", "1", "-oc", "1"]


def run_simphy(cmd: list[str]) -> None:
    Path(cmd[cmd.index("-o") + 1]).parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True)


def replicate_seed(seed: int, k: int) -> int:
    return seed + k * 10_000


def batch_seed(seed: int, k: int, b: int) -> int:
    # unique per (replicate, batch) as long as batches < 1_000_000
    return seed + k * 1_000_000 + b


def natural_key(path: Path):
    return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", str(path))]


def _read_replicate(rep_raw: Path, n_trees: int):
    # batch_<b>/ folders, if any, are read in natural order (batch_2 before batch_10)
    stree_files = sorted(rep_raw.rglob("s_tree.trees"), key=natural_key)
    if not stree_files:
        raise FileNotFoundError(f"no s_tree.trees under {rep_raw}")
    genes = []
    for f in sorted(rep_raw.rglob("g_trees*.trees"), key=natural_key):
        lines = [x for x in f.read_text().splitlines() if x.strip()]
        if len(lines) != 1:
            raise ValueError(f"{f}: expected one tree, found {len(lines)}")
        genes.append(lines[0])
    if len(genes) != n_trees:
        raise ValueError(f"{rep_raw}: expected {n_trees} gene trees, found {len(genes)}")
    return stree_files[0].read_text().strip(), genes


def write_replicates(raw_dir, mapping: dict, out_root, config: str, network: str, params: dict, replicates: int,
                     n_trees: int, retention: float | None = None, seed: int = 0, batches: int = 1) -> None:
    # SimPhy output -> <out_root>/simulated/<config>/<network>/rep<k>/, labels restored to taxon names
    raw_dir = Path(raw_dir)
    for k in range(1, replicates + 1):
        rep_raw = raw_dir / str(k)
        if not rep_raw.is_dir():
            raise FileNotFoundError(f"missing SimPhy replicate {rep_raw}")
        stree, genes = _read_replicate(rep_raw, n_trees)
        if retention is not None:
            genes = apply_fractionation(genes, stree, retention, replicate_seed(seed, k))
        out = Path(out_root) / "simulated" / config / network / f"rep{k}"
        out.mkdir(parents=True, exist_ok=True)
        restored = [restore_simphy_labels(g, mapping) for g in genes]
        (out / "gene_trees.nwk").write_text("\n".join(restored) + "\n")
        (out / "simphy_gene_trees.nwk").write_text("\n".join(genes) + "\n")
        (out / "simphy_species_tree.nwk").write_text(stree + "\n")
        info = {**params, "config": config, "network": network, "replicate": k, "n_trees": n_trees,
                "retention": retention, "seed": replicate_seed(seed, k)}
        if batches > 1:
            info["batches"] = batches
            info["batch_seeds"] = ",".join(str(batch_seed(seed, k, b)) for b in range(1, batches + 1))
        (out / "simphy_params.txt").write_text("".join(f"{key}={v}\n" for key, v in info.items()))


def simulate_network(data_root, work_dir, config: str, network: str, *, simphy_bin: str, subst_rates_pkl, seed: int,
                     n_trees: int = 1000, replicates: int = 5, height: float = 10_000_000, batches: int = 1) -> None:
    # batches > 1: SimPhy runs batches times per replicate with n_trees / batches trees each, for networks where
    # one run of all trees fails
    if batches < 1 or n_trees % batches:
        raise ValueError(f"--batches {batches} must be >= 1 and divide --n-trees {n_trees}")
    c = CONFIGS[config]
    used = [replicate_seed(seed, k) if batches == 1 else batch_seed(seed, k, b)
            for k in range(1, replicates + 1) for b in range(1, batches + 1)]
    if len(set(used)) != len(used):
        raise ValueError(f"SimPhy seeds collide for {network} (seed={seed}, replicates={replicates}, batches={batches})")
    multree = (Path(data_root) / "model_networks" / f"{network}.tre").read_text().strip()
    nex, mapping = prepare_species_tree(multree, height, seed)
    raw = Path(work_dir) / c.get("base", config) / network
    params = {"ne": c["ne"], "dup": c["dup"], "loss": c["loss"], "height": height}
    if "base" in c:
        if batches > 1:
            raise ValueError(f"{config} reuses the SimPhy output of {c['base']}, --batches does not apply")
        if not raw.is_dir():
            raise FileNotFoundError(f"{config} needs the SimPhy output of {c['base']} at {raw}")
        params["base_config"] = c["base"]
        base_params = Path(data_root) / "simulated" / c["base"] / network / "rep1" / "simphy_params.txt"
        rates = [x.split("=", 1)[1] for x in base_params.read_text().splitlines() if x.startswith("subst_rate=")]
        if len(rates) != 1:
            raise ValueError(f"{base_params}: expected one subst_rate line, found {len(rates)}")
        params["subst_rate"] = rates[0]
    else:
        raw.mkdir(parents=True, exist_ok=True)
        species_file = raw / "species_tree.nex"
        species_file.write_text(nex)
        for k in range(1, replicates + 1):
            # one SimPhy call per replicate (or batch) with its own seed, one substitution rate shared by all
            if batches == 1:
                runs = [(raw / str(k), n_trees, replicate_seed(seed, k))]
            else:
                runs = [(raw / str(k) / f"batch_{b}", n_trees // batches, batch_seed(seed, k, b))
                        for b in range(1, batches + 1)]
            for out, n, run_seed in runs:
                cmd = simphy_command(species_file, out, c["ne"], c["dup"], c["loss"], n, 1, seed=run_seed,
                                     rate_seed=seed, subst_rates_pkl=subst_rates_pkl, simphy_bin=simphy_bin)
                params["subst_rate"] = cmd[cmd.index("-su") + 1]
                run_simphy(cmd)
    write_replicates(raw, mapping, data_root, config, network, params, replicates, n_trees,
                     c.get("retention"), seed, batches)
