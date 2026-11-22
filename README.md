# polyploid-benchmark

`polybench` runs phylogenetic network inference methods for polyploids on simulated and empirical gene tree
datasets and scores their networks against the truth. It covers the whole benchmark of the accompanying
paper: SimPhy simulation of gene trees from model networks, one adapter per method with a common input and
output format, evaluation with the reticulation descendants and reticulation sisters measures and the
mu-distance, and the figures and tables of the paper. New methods are added as a folder with a YAML file and a
shell script, and are scored the same way as the published ones.

Paper: [citation to be added on publication]

## Install

```bash
conda env create -f environment/core.yml   # env "polybench", Python 3.11
conda activate polybench
pip install -e .
bash environment/install_tools.sh          # downloads the external tools into tools/ and prints the tools: block
```

If the conda solve runs out of memory, create a bare env and let pip install the rest:

```bash
conda create -n polybench -c conda-forge python=3.11 openjdk=11
conda activate polybench
pip install -e . "numpy<2" "pandas<3" scipy==1.12.0 psutil tqdm pyvis
```

Copy `config.example.yaml` to `config.yaml`, set `data_root` and `results_root` and paste the `tools:` block.
Relative paths in a config are relative to the config file. Only the tools of the methods you run are needed.
Methods with their own conda env (`env` in `method.yaml`) need it created from `environment/`:

```bash
conda env create -f environment/polyphest.yml
conda env create -f environment/mpsugar.yml
conda env create -f environment/alloppnet.yml
```

and `executor.conda_activate` set in the config (e.g. `source /path/to/conda.sh && conda activate`, the env
name is appended). `run.sh` scripts need bash, and `install_tools.sh` fetches x86-64 Linux or macOS builds.
Tools that need a manual download (PADRE, the AlloppDT scripts, PhyloNet) are listed by `install_tools.sh` with
instructions.

## Quickstart on the toy example

The toy example (`examples/toy/`) has two small simulated datasets and one empirical dataset, each with a
species tree, so it needs no external tool. The method template (`methods/_template`) returns the species tree
as the network. From the repository root:

```bash
polybench validate --config examples/toy/config.yaml --datasets 'simulated/*/*/rep*' 'empirical/*'
polybench run --config examples/toy/config.yaml --method methods/_template --datasets 'simulated/*/*/rep*' 'empirical/*'
polybench evaluate --config examples/toy/config.yaml --methods methods/_template
polybench paper --summary examples/toy/results/summary --networks examples/toy/data/model_networks --out toy_paper --only fig5 table_s2 --configs toy
```

Results go to `examples/toy/results/my_method/`, the summary CSVs to `examples/toy/results/summary/` and the
completion figure and table to `toy_paper/`. The commands work from any directory when the paths to the
config, the method and the summary are adjusted.

## Reproduce the paper

The Dryad package (doi:10.5061/dryad.XXXXXXX) holds the input data (`data/`, in the layout below), the result
of every method on every dataset (`results/<label>/`) and the summary (`results/summary/`). The simulated data
(`data/simulated/`) and the results on it (`results/<label>/simulated/`) are packed per configuration in
`archives/*.tar.gz`; extract them in the package folder before evaluating:

```bash
(cd dryad && for f in archives/*.tar.gz; do tar -xzf "$f"; done)
```

The figures, tables and statistics of the paper are made from the summary alone, which needs no extraction:

```bash
polybench paper --summary dryad/results/summary --networks dryad/data/model_networks --out paper
```

As in the paper, a GRAMPAIter run that failed or timed out is scored with the GRAMPA result of the same
dataset. `evaluate` scores these GRAMPA networks with GRAMPAIter's settings and writes them to
`summary/fallback.csv`. `--no-grampa-iter-fallback` leaves those runs unscored instead. `--only` makes a subset
of the outputs (`fig3_table3 table2 fig5 fig6_fig7 table4 stats section33 table_s2`).

To score the results again, extract the archives, write a config with `data_root: dryad/data` and
`results_root: dryad/results` and run

```bash
polybench evaluate --config config.yaml --workers 16 --nrf conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M
```

A full rerun simulates the data with `polybench simulate` (SimPhy 1.0.2, `--alignments` adds AliSim
alignments with IQ-TREE 3.0.1; the substitution rate and GTR parameter files are in the Dryad package) and runs
every method on a SLURM cluster. SimPhy can fail to simulate all gene trees of a large network in one run.
`--batches N` then runs SimPhy N times per replicate with `--n-trees`/N trees each, with one seed per batch
recorded in `simphy_params.txt`, and joins the gene trees in batch order. The simulated data of the Dryad
package were made this way for the networks where a single 1000-tree run failed (`--batches 1000`, one tree per
run); their `gene_trees.nwk` keep the tree order the methods received (see the Dryad README). A failing batch stops that network (no automatic retries), rerun it with another `--seed`.

```bash
polybench run --config config.yaml --method grampa_iter --datasets 'simulated/*/*/rep*' 'empirical/*' --executor slurm
polybench run --config config.yaml --method polyphest --params methods/polyphest/examples/p50.yaml --label polyphest_p50 \
    --datasets 'simulated/*/*/rep*' 'empirical/*' --executor slurm
```

Each call submits one job array, with time, memory and CPUs taken from the method's `method.yaml` and the other
`#SBATCH` lines from `executor.slurm_header`. The jobs run `$POLYBENCH_PYTHON -m polybench.runner` (or `python`
if unset) and sbatch passes on the submitting shell's environment, so submit from a shell where the
`polybench` env is active. `--dry-run` writes the sbatch script without submitting it. Polyphest is run once per
percentile (`p50`, `p70`, `p90`), the paper reports the lowest percentile that completed. The GRAMPA results of
the paper are the first iteration of the GRAMPAIter runs (`grampa` method, `mode: iter_first`).

## Adding a method

See [docs/adding_a_method.md](docs/adding_a_method.md): copy `methods/_template`, fill in `method.yaml`, write
`run.sh`, run it on the toy example and compare it with the published methods.

## Data layout

```
<data_root>/model_networks/<network>.tre            true networks of the simulations, as MUL-trees
<data_root>/simulated/<config>/<network>/rep<k>/    one simulated dataset
<data_root>/empirical/<dataset>/                    one empirical dataset
<results_root>/<label>/<dataset_id>/                network.nwk, status.json, log.txt
<results_root>/summary/                             inventory.csv, comparisons.csv, fallback.csv,
                                                    empirical_pairwise.csv, aggregated.csv, nrf.csv
```

A dataset folder holds

- `gene_trees.nwk`: rooted gene trees, one Newick per line. A leaf label is the taxon name, and the copies of a
  polyploid's genes repeat the name (`(A,(B,(A,C)));`). Labels cannot contain whitespace or `: , ( ) ; [ ]`.
- `ploidy.tsv` (optional): `taxon<TAB>copies`, 1 for a diploid, 2 for a tetraploid. Inferred from gene copy
  numbers when a method needs it and the file is absent.
- `species_tree.nwk` (optional): a singly labelled species tree. Inferred with ASTRAL IV (`tools.astral`) when a
  method needs it and the file is absent.
- `alignments/<gene>.fasta` (only for AlloppNET) and optionally `alignments/taxon_map.tsv`
  (`sequence_id<TAB>taxon`). A sequence id that is not in the map is taken as the taxon name.
  An optional third column `individual:A` or `individual:B` sets the homeolog pairing of a tetraploid's
  sequences (all of a taxon's sequences or none, unpaired genomes get a missing partner), otherwise
  sequences are paired gene by gene in sorted id order.
- `published_network.nwk` (empirical, optional): a published network the inferred ones are compared with.

`polybench validate --config config.yaml --datasets ... [--method <name>]` checks these rules and lists every
problem. `polybench prep` writes the inputs a method would get to a separate folder, for inspection.

The empirical gene trees of the paper were inferred with IQ-TREE and rooted, with polytomies resolved, with
NOTUNG 2.9.1.5. Where a gene tree held more copies of a taxon than its ploidy, the most divergent copies were
kept, which is done by

```bash
python -m polybench.prep.copy_selection <gene_trees.nwk> <ploidy.tsv> <out.nwk>
```

(in each gene tree, as many copies of a taxon as `ploidy.tsv` gives, chosen to be far apart in the tree). These
steps were done outside the pipeline, and the Dryad package holds their output.

## Methods

| Method | Name | Version | Inputs | Threads | Reference |
|---|---|---|---|---|---|
| PADRE | `padre` | command-line 2.0 | gene trees | 1 | Lott et al. 2009, Bioinformatics 25:1199 |
| GRAMPA | `grampa` | GRANDMA 86abd14, first iteration (GRAMPA v1.4.1 with `mode: original`) | gene trees, species tree | 8 | Thomas et al. 2017, Syst. Biol. 66:1007 |
| GRAMPAIter | `grampa_iter` | GRANDMA 86abd14 (github.com/Roshex/Grandma) | gene trees, species tree | 8 | github.com/Roshex/Grandma |
| Polyphest | `polyphest` | unpinned, upstream main (github.com/NakhlehLab/Polyphest) | gene trees, ploidy | 1 | Yan and Nakhleh, github.com/NakhlehLab/Polyphest |
| MPAllopp | `mpsugar` | MP-SUGAR in PhyNetPy 0.0.8 | gene trees | 1 | Yan et al. 2022, Syst. Biol. 71:706 |
| AlloppNET | `alloppnet` | BEAST 1.8.2, TreeAnnotator 1.8.2, AlloppDT | alignments, ploidy (diploids and tetraploids only) | 8 | Jones et al. 2013, Syst. Biol. 62:467 |
| MPAllopp (PhyloNet) | `mpallopp_phylonet` | PhyloNet 3.8.2 | gene trees | 12 | Yan et al. 2022, Syst. Biol. 71:706 (algorithm), Wen et al. 2018, Syst. Biol. 67:735 (PhyloNet) |

Species trees are inferred with ASTRAL IV from ASTER v1.23. In the paper MPAllopp was run through MP-SUGAR
(`mpsugar`). The `mpallopp_phylonet` adapter runs PhyloNet's `InferNetwork_MP_Allopp` instead; it was not used
in the paper and has not been tested against real PhyloNet output. Polyphest is run with a near-identical
subtree threshold of 0.2 and percentiles 50, 70 and 90 (`methods/polyphest/examples/`). AlloppNET runs BEAST
for a chain of 10^8 iterations or until the time limit and summarizes the sampled MUL-trees with TreeAnnotator
after 10% burn-in. Every method has a runtime limit of 5 days.

## License

MIT, see [LICENSE](LICENSE).

## How to cite

[citation to be added on publication]

Software: doi:10.5281/zenodo.XXXXXXX. See also [CITATION.cff](CITATION.cff).
