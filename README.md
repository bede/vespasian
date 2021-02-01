[![Python](https://img.shields.io/pypi/v/vespasian.svg)](https://pypi.org/project/vespasian/)


# Vespasian

Vespasian performs genome scale detection of site and branch-site signatures of positive selection by orchestrating the execution and interpretation of evolutionary hypothesis tests. Given a collection of alignments of protein-coding orthologous gene families and labelled trees, Vespasian infers gene trees from a species tree and evaluates site and lineage-specific models of evolution using PAML. Model testing is CPU intensive but embarrassingly parallel, and is  orchestrated on either a single machine or a remote compute cluster using [snakemake](https://github.com/snakemake/snakemake).

Vespasian is the pure Python successor to [VESPA](https://peerj.com/articles/cs-118/) by Webb et al. (2017).



## Installation

### With `conda` & `pip`

I currently recommend creating a conda environment with `paml` and `datrie`, and letting pip do the rest. PAML output changes in subtle and horrible ways between versions, and I have tested against builds [`h01d97ff_5`](https://anaconda.org/bioconda/paml/4.9/download/osx-64/paml-4.9-h01d97ff_5.tar.bz2) (Darwin) and [`h516909a_5`](https://anaconda.org/bioconda/paml/4.9/download/linux-64/paml-4.9-h516909a_5.tar.bz2) (Linux)

```bash
conda create -n vespasian python=3 paml
conda activate vespasian
pip install vespasian
```


### Without Conda (doesn't install PAML)

```bash
pip install vespasian
```

Alternatively, from a [PyPI](https://pypi.org/project/vespasian/#history) or Github release tarball 

```bash
tar xzf vespasian-0.3.0.tar.gz
pip install /path/to/vespasian-0.3.0/
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
- [ ] IQ-TREE / ParGenes integration
- [ ] Python API
- [ ] Specify site and/or branch-site models
- [ ] Renaming:
  - [ ] `infer-gene-trees` -> `infer-trees`
- [ ] Choice of range of omega values
- [ ] Benjamini-Hochberg correction of multiple tests