# Stand-in for PhyNetPy's MP_SUGAR: returns a star network over the taxon map keys and logs each call
import json
from pathlib import Path


class _Net:
    def __init__(self, newick):
        self.newick = newick

    def to_newick(self):
        return self.newick


def MP_SUGAR(nexus_path, taxon_map, iter_ct, seed):
    text = Path(nexus_path).read_text()
    missing = [l for labels in taxon_map.values() for l in labels if l not in text]
    if missing:
        raise ValueError(f"taxon map labels not in the trees: {missing}")
    log = Path(nexus_path).with_name("mp_sugar_calls.json")
    calls = json.loads(log.read_text()) if log.is_file() else []
    log.write_text(json.dumps(calls + [{"taxon_map": taxon_map, "iter_ct": iter_ct, "seed": seed}]))
    return {_Net("(" + ",".join(sorted(taxon_map)) + ");"): 7}
