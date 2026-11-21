import json
import sys
from pathlib import Path

import yaml

# method.yaml bool_flags: how a boolean param reaches the tool
#   value:  true -> --key True, false -> left out (Polyphest reads any non-empty string as True)
#   switch: true -> --key, false -> left out (argparse store_true, e.g. GRANDMA)
BOOL_FLAGS = ("value", "switch")


def cli_flags(params: dict, bool_flags: str = "value") -> list[str]:
    # {"percentile": 50, "m": "split"} -> --percentile 50 -m split, one-letter keys become -k
    if bool_flags not in BOOL_FLAGS:
        raise ValueError(f"bool_flags must be one of {list(BOOL_FLAGS)}, got {bool_flags!r}")
    out = []
    for key, value in params.items():
        flag = f"-{key}" if len(key) == 1 else f"--{key}"
        if value is False:
            continue
        if value is True:
            out += [flag] if bool_flags == "switch" else [flag, "True"]
        elif value is None or isinstance(value, (dict, list)):
            raise ValueError(f"param {key}: cannot pass {value!r} on the command line")
        else:
            out += [flag, str(value)]
    return out


# python -m polybench.prep.params flags <params.json> <method.yaml> [key to leave out ...]   one token per line
# python -m polybench.prep.params get <params.json> <key>
if __name__ == "__main__":
    if len(sys.argv) < 4 or sys.argv[1] not in ("flags", "get"):
        sys.exit("usage: python -m polybench.prep.params {flags <params.json> <method.yaml> [key ...] | "
                 "get <params.json> <key>}")
    sys.stdout.reconfigure(newline="\n")  # run.sh reads this, keep \r out on Windows
    cmd, path, *rest = sys.argv[1:]
    params = json.loads(Path(path).read_text())
    if cmd == "flags":
        method_yaml, *skip = rest
        style = (yaml.safe_load(Path(method_yaml).read_text()) or {}).get("bool_flags", "value")
        print("\n".join(cli_flags({k: v for k, v in params.items() if k not in skip}, style)))
    else:
        if len(rest) != 1 or rest[0] not in params:
            sys.exit(f"param {rest} not found in {path}")
        value = params[rest[0]]
        # 1e8 in a params file is the float 100000000.0, tools such as BEAST need the integer
        print(int(value) if isinstance(value, float) and value.is_integer() else value)
