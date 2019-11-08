# Vespasian

*Vespasian* performs genome scale detection of site and branch-site signatures of positive selection by orchestrating the execution and interpretation of evolutionary hypothesis tests. Given a collection of fasta alignments of orthologous gene families and labelled trees, Vespasian infers gene trees from a species tree and evaluates site and lineage-specific models of evolution using PAML. Model testing is CPU intensive but embarrassingly parallel, and may be run locally or using a cluster via Snakemake.

Vespasian is a faster and more user friendly rewrite of [VESPA](https://peerj.com/articles/cs-118/) by Webb et al. (2017).



## Status

### `infer-gene-trees` ✅

- Create per-family trees by pruning supplied input tree

### `codeml-setup` ✅ 🧵 🐍

- Configure model tests

### `report` ✅ *beta*

- Produce tables ofw parameter estimates and likelihood ratio tests between models



### To do

- [ ] Positively selected site visualisation
- [ ] IQ-TREE / ParGenes integration
- [ ] Python API
- [ ] Specify site and/or branch-site models
- [ ] Renaming:
- [ ] `infer-gene-trees` -> `infer-trees`
- [ ] Choice of range of omega values



## Installation

### With `conda`

```bash
# Recommended way until I make a conda package
conda create -n vespasian python=3 paml parallel
conda activate vespasian
pip install vespasian
```

Using bundled environment:
```bash
conda env create -f conda.yml
```
### From PyPI

```bash
# Manually install PAML
pip install vespasian
```

### From tarball

```bash
# Manually install PAML
tar xzf vespasian-0.2.0.tar.gz
pip install /path/to/vespasian-0.2.0/
```



## Usage

### Step 1: gene tree inference from a species tree

- Required input (please read carefully):
  - **`input`** Path to directory containing orthologous gene families as individual nucleotide alignments in fasta format with a `.fasta` or `.fa` extension. These should be in frame and free from stop codons. Fasta headers should contain a taxonomic identifier (mirroring tip labels in the tree file), optionally followed by separator character ('`|`' by default). A minimum of seven taxa must be present.
  - **`tree`** Path to species tree in Newick format. Tip labels must correspond to fasta headers before the separator character.
- Output:
  - Directory (default name `gene-trees`) containing minimal gene trees for each family.

```
$ vespasian infer-gene-trees -h
usage: vespasian infer-gene-trees [-h] [-o OUTPUT] [-s SEPARATOR] [-w] [-p]
                              input tree

Create gene trees by pruning a given species tree

positional arguments:
  input                 path to directory containing gene families
  tree                  path to newick formatted species tree

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        path to output directory (default: 'gene-trees')
  -s SEPARATOR, --separator SEPARATOR
                        character separating taxon name and identifier(s)
                        (default: '|')
  -w, --warnings        show warnings (default: False)
  -p, --progress        show progress bar (default: False)
```



### Step 2: Configure model testing environments

- Required input:

  - **`input`** Path to directory containing *aligned* orthologous gene families as individual fasta files.
  - **`gene-trees`** Path to directory containing minimal gene trees.

- Optional input:

  - **`—-branches BRANCHES`**  Path to yaml file containing a *mapping* of lineages to be labelled for evaluation of lineage-specific evolutionary signal. To label an individual leaf node taxon, specify its name followed by a colon. To label an internal node, choose a suitable name (e.g. `carnivora`) followed by a colon and its corresponding leaf nodes inside square brackets (a *sequence* in yaml) and separated by commas. For internal nodes, all child nodes present in the species tree must be specified, even if they are not present in all of the gene families.

    - ```yaml
      cat:
      carnivora: [cat, dog]
      ```

- Output:

  - Directory (default name `codeml`) containing nested directory structure of models and starting parameters for each gene family.
  - File `codeml-commands.sh` containing list of commands to execute the model tests
  - File `Snakefile` for running the contents of `codeml-commands.sh` locally or using a cluster

  

  N.B. By default, at least two taxa must be present within a given family for a named internal node to be labelled. Use `--strict` to skip named internal nodes unless all child leaf nodes are present. 

```
$ vespasian codeml-setup -h
usage: vespasian codeml-setup [-h] [-b BRANCHES] [-o OUTPUT]
                          [--separator SEPARATOR] [--strict] [-w] [-p]
                          input gene-trees

Create suite of branch and branch-site codeml environments

positional arguments:
  input                 path to directory containing aligned gene families
  gene-trees            path to directory containing gene trees

optional arguments:
  -h, --help            show this help message and exit
  -b BRANCHES, --branches BRANCHES
                        path to yaml file containing branches to be labelled
                        (default: -)
  -o OUTPUT, --output OUTPUT
                        path to output directory (default: 'codeml')
  --separator SEPARATOR
                        character separating taxon name and identifier(s)
                        (default: '|')
  --strict              label only branches with all taxa present in tree
                        (default: False)
  -w, --warnings        show warnings (default: False)
  -p, --progress        show progress bar (default: False)
```



### Step 3: Run models

- Ensure `codeml` binary is present inside `$PATH`
- `cd codeml` (the directory created by `codeml-setup` above)

- *Local execution (for small jobs)* 
- Using Snakemake
  - `cd` into directory generated by ` vespasian codeml-setup`
  - `snakemake --cores 8`
  - Or, using GNU parallel
  - `parallel --bar :::: codeml-commands.sh`
  
- *Cluster execution*
  - `snakemake --cores MAXJOBS --cluster OPTIONS`
- SGE example:
      - `snakemake --jobs 100 --cluster "qsub -cwd -V"`
    - Rescomp: `qsub -cwd -V -P bag.prjc -q short.qc`
  
    - Profiles [are available for other cluster platforms](https://snakemake.readthedocs.io/en/stable/executable.html#profiles)



### Step 4: Analyse model output, find positively selected sites

- ~~For now it is necessary to use the original [VESPA](https://github.com/aewebb80/vespa) for this purpose~~
- NEW: use `vespasian report` to extract positively selected sites from your codeml directory
- Documentation coming soon