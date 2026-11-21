# [TITLE OF THE DATASET]

[Authors]. [One-paragraph description: what the data are, the paper they support (citation, DOI), and what
can be done with them.]

Data for: [paper citation and DOI]. Code: polyploid-benchmark, doi:10.5281/zenodo.XXXXXXX (snapshot in `code/`).

## Contents

```
README.md
data/
  model_networks/<network>.tre          the 21 model networks of the simulations, as MUL-trees scaled to a
                                        root-to-tip height of 10,000,000 generations
  empirical/<dataset>/                  one folder per empirical dataset
archives/
  simulated_<config>.tar.gz             data/simulated/<config>/<network>/rep<k>/, one archive per configuration
  results_<config>.tar.gz               results/<label>/simulated/<config>/..., one archive per configuration
results/
  <label>/empirical/<dataset>/          network.nwk, status.json of each method on each empirical dataset
  summary/                              inventory.csv, comparisons.csv, fallback.csv, aggregated.csv,
                                        empirical_pairwise.csv, nrf.csv: the input of every figure and table
distributions/
  substitution_rates.pkl                empirical substitution rates SimPhy samples from
  gtr_parameters_all.pkl                empirical GTR parameters AliSim samples from
supplement/
  Supplementary_Methods.pdf, Table_S1.csv, Table_S2.csv, Fig_S1.pdf ... Fig_S5.pdf
code/
  polyploid-benchmark-<version>.zip     the tagged code release
```

Extract the archives in this folder (`for f in archives/*.tar.gz; do tar -xzf "$f"; done`) to get the full
`data/` and `results/` trees. Configurations (`<config>`): `conf_ils_{low,medium,high}_10M` (ILS only, Ne 200,000,
1,000,000 and 2,000,000), `conf_dup_loss_{low,medium,high}_10M[_ne1M|_ne2M]` (duplication and loss rates 1e-9,
1e-8 and 1e-7 per gene per generation, at Ne 200,000, 1,000,000 or 2,000,000) and
`conf_dup_loss_medium_10M_ne1M_fix{075,050,025}` (homeolog retention 0.75, 0.50 and 0.25 applied to the medium
ILS, medium duplication and loss gene trees, Morales-Briones_2021 not included). 1000 gene trees and 5 replicates
per network and configuration.

Results labels (`<label>`): `padre`, `polyphest_p50`, `polyphest_p70`, `polyphest_p90` (Polyphest at percentile
50, 70, 90), `grampa` (the first GRAMPA iteration of GRAMPAIter), `grampa_iter`, `mpsugar` (MPAllopp, run with
MP-SUGAR), `alloppnet` [+ `mpallopp_phylonet` if packed].

### Which file feeds which result

| Paper output | Files |
|---|---|
| Fig 5 (completion rates) | `results/summary/inventory.csv` |
| Fig 6 (accuracy under ILS and duplication/loss) | `results/summary/comparisons.csv`, `results/summary/fallback.csv` |
| Fig 7 (accuracy under homeolog loss) | `results/summary/comparisons.csv`, `results/summary/fallback.csv` (the `_fix` configurations) |
| Table 2 (simulation settings, gene tree discordance) | `results/summary/nrf.csv`, from `simphy_*.nwk` of the ILS-only configurations |
| Table 4 (tetraploid networks) | `results/summary/inventory.csv`, `results/summary/comparisons.csv`, `results/summary/fallback.csv` |
| Section 3.3 (empirical data) | `results/summary/inventory.csv`, `results/summary/empirical_pairwise.csv`, `results/<label>/empirical/` |
| Table S2 (model network statistics) | `data/model_networks/`, `supplement/Table_S2.csv` |
| Statistical tests of the Results | `results/summary/comparisons.csv`, `results/summary/fallback.csv`, `results/summary/inventory.csv` |

## Formats

- Newick files hold one tree per line, rooted. A leaf label is the taxon name, and the copies of a polyploid
  repeat the name, in gene trees and in MUL-trees (`(A,(B,(A,C)));`). Labels contain no whitespace or
  `: , ( ) ; [ ]`. `network.nwk` is the method output: a MUL-tree, or an extended Newick network (`#H1` nodes)
  for MPAllopp.
- `data/simulated/<config>/<network>/rep<k>/`
  - `gene_trees.nwk`: the SimPhy gene trees with taxon names as leaf labels, in the order the methods of the
    paper received them. That is the shell glob order of the SimPhy files: where SimPhy failed on 1000 trees at
    once it was run 1000 times for one tree each (`batch_1` ... `batch_1000`), and `batch_10` comes before
    `batch_2`. Branch lengths are written with 6 significant figures; `simphy_gene_trees.nwk` holds the raw
    SimPhy trees.
  - `simphy_gene_trees.nwk`, `simphy_species_tree.nwk`: the SimPhy output with SimPhy's labels
    (`<taxon>_<copy>` in the species tree, `<species tree leaf>_<locus>_<individual>` in the gene trees).
  - `simphy_params.txt`: `key=value` lines: `ne`, `dup`, `loss` (per gene per generation), `height`,
    `subst_rate`, `config`, `network`, `replicate`, `n_trees`, `retention` (`None` without homeolog loss),
    `base_config` and `retention_seed` for the retention configurations, `legacy_base_seed` (the SimPhy base seed).
  - `alignments/gene_<n>.fasta` and `alignments/taxon_map.tsv` (`sequence_id<TAB>taxon<TAB>individual:genome`,
    the third column only for tetraploids: the AlloppNET homeolog pairing of the paper runs), only for the datasets
    AlloppNET was run on (the 8 networks with only diploids and tetraploids).
- `data/empirical/<dataset>/`: `gene_trees.nwk`, `ploidy.tsv` (`taxon<TAB>copies`, the number of subgenome
  copies: 1 for a diploid, 2 for a tetraploid, 3 or more allowed for every method but AlloppNET) when copy numbers
  were given or inferred, `species_tree.nwk` (ASTRAL species tree) when the methods
  used one, `published_network.nwk` when the study published a network, `alignments/` for AlloppNET.
- `status.json`: `{"status": "completed" | "failed", "runtime_s": null, "method", "version", "params", "source"}`.
  `failed` means the method did not produce a network within the 5-day limit or crashed. A dataset without a
  result folder for a label was not run with that method. Runtimes are not recorded (`null`). A GRAMPAIter run
  that did not finish has status `failed` and no `network.nwk`; as in the paper, the analysis scores it with the
  GRAMPA result of the same dataset (`results/summary/fallback.csv`).
- `results/summary/`: `inventory.csv` (`kind,config,network,replicate,label,status,runtime_s`), `comparisons.csv`
  (`kind,config,network,replicate,label,metric,value,status`, one row per measure of each inferred network
  against the model network), `fallback.csv` (same columns: the GRAMPA network of each dataset where GRAMPAIter
  failed, scored as a GRAMPAIter result, label `grampa_iter`, status `FALLBACK`), `aggregated.csv`,
  `empirical_pairwise.csv` (`dataset,label_a,label_b,metric,value`),
  `nrf.csv` (`config,network,replicate,mean_nrf`).

## Reproducing the figures and tables

Install the code (doi:10.5281/zenodo.XXXXXXX, see its README) and run

```bash
polybench paper --summary results/summary --networks data/model_networks --out paper
```

To recompute the summary from the networks, extract the archives, write `config.yaml` with `data_root: data`, `results_root: results`,
`tools: {}` and `executor: {type: local}`, and run
`polybench evaluate --config config.yaml --workers 16 --nrf conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M`.
`polybench simulate` regenerates simulated data from the model networks and `distributions/`, and
`polybench run` runs a method on any dataset folder.

## Empirical datasets

The gene trees of the empirical datasets were inferred with IQ-TREE from the published alignments, or taken
from the study, then rooted, with polytomies resolved, with NOTUNG 2.9.1.5. Copy numbers were inferred by kernel
smoothing of the per-gene copy counts, unless reported by the original study (see the paper).

If the package was built with `pack_dryad.py --no-empirical`, each empirical dataset is instead in
`empirical_as_run/<dataset>/`, holding the `gene_trees/` and `networks/` folders exactly as the analyses used and
produced them.

| Dataset | Citation | Original data source | License/redistribution | Changes made |
|---|---|---|---|---|
| Ren_2024 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Zhao_2021 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Ding_2023 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Koenen_2020 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Marcussen_2012 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Wu_2015 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Lawrence_2016 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1]; one tip kept per monophyletic group of copies (`dryad/empirical_provenance/lawrence_2016_reduce.py` in the code) |
| Popp_2005 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Diaz-Perez_2018 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1]; single-species clades collapsed (`dryad/empirical_provenance/reduce_monophyletic_groups.py`) |
| Marcussen_2015 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1]; labels cut at the first underscore (`dryad/empirical_provenance/removed_after_underscore.py`) |
| Morales-Briones_2021 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Wisecaver_2023 (comparable) | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Bendiksby_2011 | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Brysting_2007 | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Liang_2019 | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Liu_2023 | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Marcussen_2011 | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Sessa_2012b | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Hori_2014 | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Shahrestani_2015 | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |
| Soza_2014 | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] | [from Supplementary Table S1] |

"Comparable" datasets have more than one sequenced gene. The methods were run on them and they enter Section
3.3. The other datasets contribute their published networks to the model networks.

## Supplementary material

- `supplement/Supplementary_Methods.pdf`: [description]
- `supplement/Table_S1.csv`: the empirical studies, their data sources and the published networks [description]
- `supplement/Table_S2.csv`: statistics of the 21 model networks (made with `polybench paper --only table_s2`)
- `supplement/Fig_S1.pdf` ... `Fig_S5.pdf`: [one line each]

## Code

`code/polyploid-benchmark-<version>.zip` is a snapshot of the code release used for this package
(doi:10.5281/zenodo.XXXXXXX). The current version is at [repository URL].

## License

The data in this package are released under CC0 1.0 (public domain dedication), as required by Dryad. Please
cite the paper and the original studies of the empirical datasets when you use them.
