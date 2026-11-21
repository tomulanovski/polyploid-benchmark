import random
import re

from ete3 import Tree

_LEAF = re.compile(r"(?<=[(,])[^(),:;\[\]]+(?=[:),;\[])")
_INTERNAL = re.compile(r"\)[^,:;()\[\]]+")
_ROOT_LENGTH = re.compile(r"\):[0-9.eE+-]+;$")


def split_newick(text: str) -> list[str]:
    # trees separated by ';', a tree may span lines
    return [t.strip() + ";" for t in text.split(";") if t.strip()]


def map_leaves(newick: str, fn) -> str:
    return _LEAF.sub(lambda m: fn(m.group(0)), newick)


def relabel_leaves(newick: str, mapping: dict[str, str]) -> str:
    return map_leaves(newick, lambda name: mapping.get(name, name))


def make_unique_labels(taxa: list[str]) -> dict[str, str]:
    # a label that is a substring of another gets an X suffix (GRAMPA crashes otherwise)
    taxa = sorted(set(taxa), key=len)
    clash = [t for i, t in enumerate(taxa) if any(t in other for other in taxa[i + 1:])]
    mapping: dict[str, str] = {}
    for taxon in clash:
        others = [t for t in taxa if t != taxon] + list(mapping.values())
        for candidate in (f"{taxon}X", f"{taxon}Z", f"X{taxon}", f"Z{taxon}"):
            if not any(candidate in o or o in candidate for o in others):
                mapping[taxon] = candidate
                break
        else:
            raise ValueError(f"no safe unique label for {taxon}")
    final = [mapping.get(t, t) for t in taxa]
    if any(a in b for a in final for b in final if a != b):
        raise ValueError("could not make labels free of substrings")
    return mapping


def restore_labels(newick: str, mapping: dict[str, str]) -> str:
    return relabel_leaves(newick, {new: old for old, new in mapping.items()})


def _resolve(tree: Tree, rng: random.Random) -> None:
    for node in tree.traverse():
        while len(node.children) > 2:
            a, b = rng.sample(node.children, 2)
            joined = node.add_child(dist=0.0)
            for child in (a, b):
                joined.add_child(child.detach())


def clean_gene_trees(lines: list[str], seed: int = 0) -> list[str]:
    rng = random.Random(seed)
    out = []
    for line in lines:
        text = _INTERNAL.sub(")", line.strip())
        text = _ROOT_LENGTH.sub(");", text)
        fmt = 5 if ":" in text else 9
        try:
            tree = Tree(text, format=fmt)
        except Exception as err:  # ete3 raises several error types for bad Newick
            raise ValueError(f"cannot parse gene tree: {line[:60]}") from err
        _resolve(tree, rng)
        out.append(tree.write(format=fmt))
    return out
