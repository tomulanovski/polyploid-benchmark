# Runs PhyNetPy MP_SUGAR and writes every chain's networks as "Newick: ..." lines.
# Usage: python run_mpsugar.py <trees.nex> <taxon_map.json> <params.json> <out.txt>
import json
import random
import sys
from pathlib import Path

from PhyNetPy.MPSugar import MP_SUGAR

trees, taxon_map, params_json, out = sys.argv[1:]
params = json.loads(Path(params_json).read_text())
taxon_map = json.loads(Path(taxon_map).read_text())
first_seed = params["seed"] if params["seed"] is not None else random.randint(0, 1000000)

lines = [f"MP-SUGAR iterations {params['iterations']}, chains {params['chains']}"]
for chain in range(1, int(params["chains"]) + 1):
    seed = first_seed + chain - 1
    networks = MP_SUGAR(trees, taxon_map, iter_ct=int(params["iterations"]), seed=seed)
    if not networks:
        sys.exit(f"MP_SUGAR returned no network for chain {chain}")
    lines.append(f"Chain {chain} (seed: {seed})")
    for i, (net, score) in enumerate(networks.items(), 1):
        lines += [f"Network {i}:", f"Newick: {net.to_newick()}", f"Score: {score}"]
Path(out).write_text("\n".join(lines) + "\n")
