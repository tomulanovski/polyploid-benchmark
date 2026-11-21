# Stand-in for Polyphest.py: writes the consensus multiset as a star MUL-tree
import argparse
import json
import sys
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--gene_tree_file", required=True)
p.add_argument("--consensus_multiset_file", required=True)
p.add_argument("--output_dir", required=True)
p.add_argument("--filter_strategy", required=True)
a, _ = p.parse_known_args()
taxa = Path(a.consensus_multiset_file).read_text().split()
name = Path(a.gene_tree_file).name.rsplit(".", 1)[0] + "-polyphest.txt"
Path(a.output_dir, name).write_text(f"multree: ({','.join(taxa)});\nnetwork: x;\n")
Path(a.output_dir, "argv.json").write_text(json.dumps(sys.argv[1:]))
