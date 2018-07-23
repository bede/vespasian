import sys
import argh

from Bio import SeqIO

from vespa import vespa


def infer_gene_trees(input_path: 'path to directory containing gene families',
                     tree_path: 'species tree in newick format',
                     output_path):
    vespa.infer_gene_trees(input_path, tree_path, output_path)


def codeml_setup(input: 'path to directory containing gene families',
                 branches: 'path to file containing branches to be labelled as foreground lineages'):
    pass


def main():
    parser = argh.ArghParser()
    parser.add_commands([infer_gene_trees,
                         codeml_setup])
    parser.dispatch()


if __name__ == '__main__':
    main()
