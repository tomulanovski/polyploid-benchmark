import json
import subprocess
import sys

import pytest

from polybench.prep.params import cli_flags


def test_cli_flags_value_style():
    assert cli_flags({"m": "split", "cap": 15, "use_near_isomorphic": True, "use_clusters": False}) == \
        ["-m", "split", "--cap", "15", "--use_near_isomorphic", "True"]


def test_cli_flags_switch_style():
    assert cli_flags({"mode": "split", "plot": True, "debug": False}, "switch") == ["--mode", "split", "--plot"]


def test_cli_flags_rejects_none():
    with pytest.raises(ValueError, match="percentile"):
        cli_flags({"percentile": None})


def test_cli_flags_rejects_unknown_style():
    with pytest.raises(ValueError, match="bool_flags"):
        cli_flags({"a": True}, "maybe")


def test_cli_reads_bool_style_from_method_yaml(tmp_path):
    p = tmp_path / "params.json"
    p.write_text(json.dumps({"mode": "split", "plot": True}))
    m = tmp_path / "method.yaml"
    m.write_text("bool_flags: switch\n")
    res = subprocess.run([sys.executable, "-m", "polybench.prep.params", "flags", str(p), str(m)],
                         capture_output=True, text=True, check=True)
    assert res.stdout.split() == ["--mode", "split", "--plot"]
    m.write_text("name: x\n")
    res = subprocess.run([sys.executable, "-m", "polybench.prep.params", "flags", str(p), str(m)],
                         capture_output=True, text=True, check=True)
    assert res.stdout.split() == ["--mode", "split", "--plot", "True"]
