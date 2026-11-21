import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path


def parse_slurm_time(s: str) -> int:
    days, _, clock = s.rpartition("-")
    h, m, sec = (int(x) for x in clock.split(":"))
    return int(days or 0) * 86400 + h * 3600 + m * 60 + sec


def time_budget(limit_s: int, elapsed_s: float) -> float:
    # time left for the method so the runner can write status.json before SLURM kills the job:
    # 10 minutes of margin, or 5% of the limit when the limit is under 3 hours
    margin = 0.05 * limit_s if limit_s < 3 * 3600 else 600
    return max(1.0, limit_s - elapsed_s - margin)


def bash() -> str:
    exe = shutil.which("bash")
    if exe is None:
        raise FileNotFoundError("bash not found on PATH")
    return exe


def run_command(spec, cfg: dict, args: list[str]) -> list[str]:
    # run.sh, wrapped in the conda activation line when the method has an env
    script = (spec.path / "run.sh").as_posix()
    if spec.env is None:
        return [bash(), script, *args]
    activate = cfg["executor"]["conda_activate"]
    if not activate:
        raise ValueError(f"method {spec.name} needs conda env '{spec.env}' but executor.conda_activate is not set")
    line = f"{activate} {spec.env} && bash {shlex.quote(script)} " + " ".join(shlex.quote(a) for a in args)
    return [bash(), "-c", line]


def tool_env(cfg: dict) -> dict:
    env = dict(os.environ, POLYBENCH_PYTHON=sys.executable)
    for tool, value in cfg["tools"].items():
        if value is not None:
            env[f"POLYBENCH_{tool.upper()}"] = str(value)
    return env


def _kill_tree(proc: subprocess.Popen) -> None:
    if os.name == "nt":
        subprocess.run(["taskkill", "/T", "/F", "/PID", str(proc.pid)], capture_output=True)
    else:
        try:
            os.killpg(proc.pid, 9)
        except ProcessLookupError:
            pass
    proc.wait()


def run_with_timeout(cmd: list[str], env: dict, log_path: Path, timeout_s: float) -> tuple[int | None, float]:
    # returns (exit code or None on timeout, seconds), stdout and stderr are appended to log_path
    extra = {} if os.name == "nt" else {"start_new_session": True}
    with open(log_path, "a") as log:
        start = time.monotonic()
        proc = subprocess.Popen(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, **extra)
        try:
            code = proc.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            _kill_tree(proc)
            code = None
        return code, time.monotonic() - start


def render_slurm(spec, cfg: dict, n: int, label_dir: Path, method_path: Path, params_path: Path,
                 list_path: Path, cfg_path: Path, label: str) -> str:
    if any(c.isspace() for c in str(label_dir)):
        raise ValueError(f"results_root/label contains whitespace, which #SBATCH --output cannot hold: {label_dir}")
    r = spec.resources
    q = lambda p: shlex.quote(p.as_posix() if isinstance(p, Path) else p)
    lines = ["#!/usr/bin/env bash", *cfg["executor"]["slurm_header"],
             f"#SBATCH --time={r['time']}",
             f"#SBATCH --mem={r['memory_gb']}G",
             f"#SBATCH --cpus-per-task={r['threads']}",
             f"#SBATCH --array=0-{n - 1}",
             f"#SBATCH --job-name={label}",
             f"#SBATCH --output={(label_dir / 'logs' / '%A_%a.out').as_posix()}",
             "",
             '"${POLYBENCH_PYTHON:-python}" -m polybench.runner --one'
             f" --method {q(method_path)} --label {q(label)} --params {q(params_path)}"
             f" --dataset-list {q(list_path)} --index $SLURM_ARRAY_TASK_ID"
             f" --config {q(cfg_path)}", ""]
    return "\n".join(lines)
