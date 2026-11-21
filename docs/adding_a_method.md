# Adding a method

A method is a folder under `methods/` with a `method.yaml` and a `run.sh`. Nothing else in the package needs to
change. The steps below add a method called `my_method`, run it on the toy example and put it next to the
published methods. Commands are run from the repository root.

## 1. Copy the template

```bash
cp -r methods/_template methods/my_method
```

## 2. Fill in method.yaml

Every key is described in the template. Set `name` to the folder name (results are labelled with it and
`evaluate` finds the method by it), then:

- `requires`: the inputs the method reads, any of `gene_trees`, `ploidy`, `species_tree`, `alignments`.
  A dataset without `ploidy.tsv` gets one inferred from gene copy numbers, a dataset without
  `species_tree.nwk` gets an ASTRAL IV tree (`tools.astral` in the config).
- `output`: `multree` if `network.nwk` holds a MUL-tree, `network` if it holds an extended Newick network.
- `partial_match`: `true` if the method infers only a limited number of reticulations (e.g. one, as GRAMPA),
  so the reticulation measures do not penalize the true reticulations it leaves unmatched. Leave it `false` for
  methods that infer full networks.
- `env`: a conda env name to run `run.sh` in, or `null` for the current environment.
- `resources`: threads, memory and the runtime limit. The limit applies to local runs as well as SLURM.
- `params`: default parameters. `polybench run --params file.yaml` overrides them per run.

## 3. Write run.sh

`run.sh` is called as

```
run.sh <dataset_dir> <out_dir> <params_json> <threads>
```

`<dataset_dir>` holds `gene_trees.nwk` and the inputs listed in `requires` (`ploidy.tsv`, `species_tree.nwk`,
`alignments/`), copied from the dataset or derived; no other file of the dataset is passed on. The script must
write `<out_dir>/network.nwk` with leaf labels equal to the taxon names, copies of a polyploid repeated. Anything
else it writes should go under `<out_dir>/work/`. A non-zero exit or a missing `network.nwk` is recorded as a failed run.

Tool paths from the config's `tools:` block arrive as `POLYBENCH_<TOOL>` environment variables
(`POLYBENCH_ASTRAL`, `POLYBENCH_PHYLONET_JAR`, ...) and `POLYBENCH_PYTHON` is the python that runs polybench.
A method with its own binary gets a new key in the config, e.g. `mytool: /path/to/mytool` under `tools:`,
which arrives as `POLYBENCH_MYTOOL` (keys use letters, digits and `_`).
A parameter is read with

```bash
"$POLYBENCH_PYTHON" -m polybench.prep.params get "$params_json" <key>
```

and all parameters as command line flags with `polybench.prep.params flags` (see `methods/grampa_iter/run.sh`).

Two optional hooks live next to `run.sh`:

- `prepare.py` with `prepare(dataset_dir, work_dir, params)`, called before `run.sh` to write method specific
  input files, usually under `<work_dir>/inputs` (see `methods/padre/prepare.py`).
- `validate.py` with `extra_checks(dataset_dir, params)` returning a list of problems, run by
  `polybench validate --method` (see `methods/polyphest/validate.py`).

## 4. Run it on the toy example

```bash
polybench validate --config examples/toy/config.yaml --datasets 'simulated/*/*/rep*' 'empirical/*' --method my_method
polybench run --config examples/toy/config.yaml --method my_method --datasets 'simulated/*/*/rep*' 'empirical/*'
```

Each dataset gets `examples/toy/results/my_method/<dataset_id>/` with `network.nwk`, `status.json` and `log.txt`.
To see the inputs the method gets without running it:

```bash
polybench prep --config examples/toy/config.yaml --datasets 'empirical/*' --method my_method --out prep_check
```

## 5. Evaluate

```bash
polybench evaluate --config examples/toy/config.yaml --labels my_method
```

This scores `my_method` and writes its rows to `examples/toy/results/summary/` (`inventory.csv`,
`comparisons.csv`, `empirical_pairwise.csv`, `aggregated.csv`). With `--labels` and an existing summary, only
the rows of the given labels are replaced and the rows of other labels are kept as they are. Without
`--labels` the summary is rebuilt from every folder in `results_root`. A method kept outside `methods/` is run
with `--method <folder>` and evaluated with `polybench evaluate --config ... --methods <folder>`.

```bash
polybench paper --summary examples/toy/results/summary --networks examples/toy/data/model_networks \
    --out toy_paper --only fig5 --configs toy
```

`toy_paper/fig5_completion.csv` now has a row for `my_method`.

## 6. Compare with the published methods

Download the Dryad package (doi:10.5061/dryad.XXXXXXX) and point a config at it, with `data_root` set to its
`data/` folder and `results_root` to its `results/` folder, which holds one folder per published method
label and their summary. Run the new method (on a cluster: add `--executor slurm`), add its scores to the
published summary and make the figures and tables:

```bash
polybench run --config config.yaml --method my_method --datasets 'simulated/*/*/rep*' 'empirical/*'
polybench evaluate --config config.yaml --labels my_method --workers 16
polybench paper --summary <results_root>/summary --networks <data_root>/model_networks --out paper
```

The published methods keep their scores from the package, and the empirical pairwise comparisons of the new
method against them are added. The new method appears in Fig. 5 and Table 4 after the published methods. Run
only some outputs with `--only`, e.g. `--only table4 fig5`.
