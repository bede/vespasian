# vespa-slim

Genome scale detection of signatures of positive selection from orthologous gene families through automatic model testing. Given a collection of fasta files containing orthologous gene families, infer gene trees from a species tree and evaluate site and lineage specific models of evolution using PAML's codeml. Model testing is CPU intensive but embarrassingly parallel, and may be run locally or using a cluster via Snakemake. 

Vespa-slim is intended to be a faster and more user friendly alternative to [VESPA](https://peerj.com/articles/cs-118/) by Webb et al. (2017). 



## Roadmap

### Commands

#### `infer-gene-trees` ✅

- `infer-trees`

#### `codeml-setup` ✅ 

- `configure-models` 🧵
- `test-models` 🐍

#### `codeml-reader` 🔜 

- `report`
- Still reliant on existing VESPA for reporting on model test output




### Features

- Snakemake orchestrated execution of model testing

  - Cluster execution
- IQ-TREE or ParGenes

  

## Usage

### Step 1: gene tree inference from a species tree

- Required input;
  - **`input`** Path to directory containing *aligned* orthologous gene families as individual fasta files. Fasta headers should contain a taxonomic identifier (mirroring tip labels in the tree file), followed by separator character ('`|`' by default).
  - **`tree`** Path to species tree in Newick format. Tip labels must correspond to fasta headers before the separator character.
- Output:
  - Directory (default name `gene-trees`) containing minimal gene trees for each family.

```
$ vespa infer-gene-trees -h
usage: vespa infer-gene-trees [-h] [-o OUTPUT] [-s SEPARATOR] [-w] [-p]
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

  - **`—-branches BRANCHES`**  Path to yaml file containing a *mapping* of lineages to be labelled for evaluation of lineage-specific evolutionary signal. To label an individual leaf node taxon, specify its name followed by a colon. To label an internal node, choose a suitable name (e.g. `carnivora`) followed by a colon and its corresponding leaf nodes inside square brackets (a *sequence* in yaml) and separated by commas. For internal nodes, all child nodes present in the species tree must be specified ,even if they are not present in all of the gene families.

    - ```yaml
      cat:
      carnivora: [cat, dog]
      ```

- Output:

  - Directory (default name `codeml`) containing nested directory structure of models and starting parameters for each gene family, and file containing a list of commands to execute these model tests.

  

  N.B. By default, at least two taxa must be present within a given family for a named internal node to be labelled. Use `--strict` to skip named internal nodes unless all child leaf nodes are present. 

```
$ vespa codeml-setup -h
usage: vespa codeml-setup [-h] [-b BRANCHES] [-o OUTPUT]
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
