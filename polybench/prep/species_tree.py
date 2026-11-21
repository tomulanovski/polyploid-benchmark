import subprocess


def astral_command(gene_trees_path, out_path, astral_bin: str, threads: int) -> list[str]:
    return [astral_bin, "-i", str(gene_trees_path), "-o", str(out_path), "-t", str(threads)]


def run_astral(gene_trees_path, out_path, astral_bin: str = "astral4", threads: int = 1) -> None:
    subprocess.run(astral_command(gene_trees_path, out_path, astral_bin, threads), check=True)
