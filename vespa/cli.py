import sys
import warnings

import argh

from vespa import vespa



def configure_warnings(show_warnings):
    if show_warnings:
        warnings.filterwarnings('always')
    else:
        warnings.filterwarnings('ignore')



def infer_gene_trees(input_dir: 'path of directory containing gene families',
                     tree: 'path of newick formatted species tree',
                     output_dir: 'path of output directory'):
    '''Create gene trees by pruning a given species tree'''
    vespa.infer_gene_trees(input_dir, tree, output_dir)


def codeml_setup(input: 'path to directory containing aligned gene families',
                 gene_trees: 'path to directory containing gene trees',
                 branches: 'path to yaml file containing branches to be labelled as foreground lineages',
                 output: 'path to output directory',
                 progress: 'show progress bar' = False,
                 warnings: 'show warnings' = False):
    '''Create suite of branch and branch-site codeml environments'''
    configure_warnings(warnings)
    vespa.codeml_setup(input, gene_trees, branches, output, progress)



def main():
    parser = argh.ArghParser()
    parser.add_commands([infer_gene_trees,
                         codeml_setup])
    parser.dispatch()


if __name__ == '__main__':
    main()
