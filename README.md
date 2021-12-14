[![DOI](https://zenodo.org/badge/141450954.svg)](https://zenodo.org/badge/latestdoi/141450954)
[![Tests](https://img.shields.io/github/workflow/status/bede/vespasian/tests)](https://github.com/bede/vespasian/actions)
[![PyPI](https://img.shields.io/pypi/v/vespasian.svg?color=brightgreen)](https://pypi.org/project/vespasian/)


# Vespasian

Vespasian performs genome scale detection of site and branch-site signatures of positive selection by orchestrating evolutionary hypothesis tests with PAML. Given a collection of alignments of protein-coding orthologous gene families and labelled trees, Vespasian infers gene trees from a species tree and evaluates site and lineage-specific models of evolution. Model testing is CPU-intensive but embarrassingly parallel, and can be executed on one or many machines with [snakemake](https://github.com/snakemake/snakemake). Vespasian is the pure Python successor to [VESPA](https://peerj.com/articles/cs-118/) by Webb et al. (2017).



## Installation

### With `conda` & `pip`

I currently recommend creating a conda environment with `paml` and `datrie`, and letting pip do the rest. PAML output changes subtly between versions, and I have tested against MacOS builds `h01d97ff_5`, `hb4d813b_6` and Linux builds  `h516909a_5`, `h779adbc_6` from [bioconda](https://anaconda.org/bioconda/paml/files).

```bash
conda create -n vespasian python=3 paml datrie
conda activate vespasian
pip install vespasian
```

### Without `conda`

```bash
# Install PAML manually
pip install vespasian
```

### Development

```bash
conda create -n vespasian python=3 paml && conda activate vespasian
git clone https://github.com/bede/vespasian
pip install --editable vespasian
```



## Usage

### Step 1: gene tree inference from a species tree

e.g. `vespasian infer-gene-trees --warnings --progress input tree `

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



### Step 2: Configure model test environments

e.g. `vespasian codeml-setup --progress --warnings --branches branches.yml input gene-trees`

- Required input:

  - **`input`** Path to directory containing *aligned* orthologous gene families as individual fasta files.
  - **`gene-trees`** Path to directory containing minimal gene trees.

- Optional input:

  - **`—-branches BRANCHES`**  Path to yaml file containing a YAML *mapping* of lineages to be labelled for evaluation of lineage-specific evolutionary signal using branch-site tests. To label an individual leaf node taxon, specify its name followed by a colon. To label an internal node, choose a suitable name (e.g. `carnivora`) followed by a colon and its corresponding leaf nodes inside square brackets (a *sequence* in yaml) and separated by commas. For internal nodes, all child nodes present in the species tree must be specified, even if they are not present in all of the gene families.

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
                              [--separator SEPARATOR] [--strict] [-t THREADS]
                              [-w] [-p]
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
                        (default is >= 2) (default: False)
  -t THREADS, --threads THREADS
                        number of parallel workers (default: 6)
  -w, --warnings        show warnings (default: False)
  -p, --progress        show progress bar (default: False)
```



### Step 3: Run models

e.g. `cd codeml && snakemake --cores 8`

- Ensure `codeml` binary is present inside `$PATH`
- Using PAML version `4.9=h01d97ff_5` from Conda is recommended
- `cd codeml` (the directory created by `codeml-setup` in step 2)
- *Local execution (for small jobs)* 
  - `snakemake -k --cores 8` (recommended)
  - Or, using GNU parallel (*not* recommended – doesn't catch errors!)
    - `parallel --bar :::: codeml-commands.sh`
- *Cluster execution*
  
  - `snakemake -k --cores MAXJOBS --cluster OPTIONS`
  - SGE example:
      - `snakemake -k --jobs 100 --cluster "qsub -cwd -V" --max-status-checks-per-second 0.1`
      - Oxford Rescomp: `qsub -cwd -V -P bag.prjc -q short.qc`
      - Profiles [are available for other cluster platforms](https://snakemake.readthedocs.io/en/stable/executable.html#profiles)



### Step 4: Report model tests and positively selected sites

e.g. `vespasian report --progress input`

- Required input (please read carefully):
  - **`input`** path to directory (default `codeml`) containing models configured in step 2 and executed in step 3
- Output:
  - Directory containing per-gene tables of likelihood ratio test results, model parameters, and positively selected sites from the highest scoring models.

```
$ vespasian report -h
usage: vespasian report [-h] [-o OUTPUT] [--hide] [-p] input

Perform likelihood ratio tests and and report positively selected sites

positional arguments:
  input                 path to codeml-setup output directory

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        path to output directory (default: 'report-codeml')
  --hide                hide gratuitous emperor portrait (default: False)
  -p, --progress        show progress bar (default: False)
```



### Todo

- [ ] Positively selected site visualisation
- [ ] Python API
- [ ] Specify site and/or branch-site models only
- [ ] Renaming:
  - [ ] `infer-gene-trees` -> `infer-trees`
- [ ] Consider B-H correction