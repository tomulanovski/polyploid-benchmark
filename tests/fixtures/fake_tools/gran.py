# Stand-in for GRANDMA's gran.py: echoes the species tree back as GRAMPA-style output
import argparse
import json
import re
import sys
from pathlib import Path

p = argparse.ArgumentParser()
for flag in ("-g", "-s", "-o", "-p", "-v"):
    p.add_argument(flag)
p.add_argument("-m", "--mode")
p.add_argument("-c", "--cap")
p.add_argument("-i", "--iter")
p.add_argument("--overwrite", action="store_true")
a = p.parse_args()
out = Path(a.o)
out.mkdir(parents=True, exist_ok=True)
species = Path(a.s).read_text().strip()
if not all(re.fullmatch(r"\d+_\S+", leaf) for leaf in re.findall(r"(?<=[(,])[^(),:;]+", Path(a.g).read_text())):
    sys.exit("gene tree labels are not <copy>_<species>")
marked = re.sub(r"(?<=[(,])([^(),:;]+)", r"\1*", species)
marked = re.sub(r"\)", ")<1>", marked)
(out / "history.json").write_text(json.dumps({"(0, 0)": {"best_mt": marked}}))
(out / "final_multree.tre").write_text(species)
(out / "argv.json").write_text(json.dumps(sys.argv[1:]))
