import sys
import argh

from Bio import SeqIO

from vespa import vespa


def infer_genetree(alignment: 'path to directory containing gene families',
                   tree: 'species tree in newick format',
                   output):
    vespa.infer_genetree(alignment, tree, output)


def codeml_setup(input: 'path to directory containing gene families',
                 branches: 'path to file containing branches to be labelled as foreground lineages'):
    pass


def main():
    parser = argh.ArghParser()
    parser.add_commands([infer_genetree,
                         codeml_setup])
    parser.dispatch()


if __name__ == '__main__':
    main()
