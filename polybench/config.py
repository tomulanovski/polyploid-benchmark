import os
import re
from pathlib import Path

import yaml

TOOL_KEY = re.compile(r"[A-Za-z0-9_]+")  # tools.<key> reaches run.sh as POLYBENCH_<KEY>


def load_config(path) -> dict:
    cfg = yaml.safe_load(Path(path).read_text())
    if not isinstance(cfg, dict):
        raise ValueError(f"{path}: config must be a mapping")
    for key in ("data_root", "results_root", "tools", "executor"):
        if key not in cfg:
            raise ValueError(f"{path}: missing key '{key}'")
    if not isinstance(cfg["tools"], dict):
        raise ValueError(f"{path}: tools must be a mapping")
    bad = sorted(str(k) for k in cfg["tools"] if not TOOL_KEY.fullmatch(str(k)))
    if bad:
        raise ValueError(f"{path}: tool keys {bad} must contain only letters, digits and '_'")
    ex = cfg["executor"]
    if ex.get("type") not in ("local", "slurm"):
        raise ValueError(f"{path}: executor.type must be 'local' or 'slurm', got {ex.get('type')!r}")
    ex.setdefault("slurm_header", [])
    ex.setdefault("conda_activate", "")
    # relative paths are relative to the config file; a tool given as a bare name is looked up on PATH
    base = Path(path).resolve().parent
    for key in ("data_root", "results_root"):
        cfg[key] = _absolute(base, cfg[key])
    for tool, value in cfg["tools"].items():
        if value is not None and any(sep in str(value) for sep in "/\\"):
            cfg["tools"][tool] = _absolute(base, value)
    return cfg


def _absolute(base: Path, value) -> str:
    return str(value) if os.path.isabs(value) else str((base / value).resolve())
