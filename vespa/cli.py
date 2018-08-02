import sys
import warnings

import argh

from vespa import vespa



def configure_warnings(show_warnings):
    '''Show or suppress warnings, mainly for TreeSwift Tree.mrca() operations'''
    if show_warnings:
        warnings.filterwarnings('always')
    else:
        warnings.filterwarnings('ignore')



def infer_gene_trees(input: 'path to directory containing gene families',
                     tree: 'path to newick formatted species tree',
                     output: 'path to output directory' = 'gene-trees'):
    '''CLI: create gene trees by pruning a given species tree'''
    vespa.infer_gene_trees(input, tree, output)


def codeml_setup(input: 'path to directory containing aligned gene families',
                 gene_trees: 'path to directory containing gene trees',
                 branches: 'path to yaml file containing branches to be labelled' = None,
                 output: 'path to output directory' = 'codeml',
                 progress: 'show progress bar' = False,
                 warnings: 'show warnings' = False):
    '''CLI: create suite of branch and branch-site codeml environments'''
    configure_warnings(warnings)
    vespa.codeml_setup(input, gene_trees, branches, output, progress)



def main():
    parser = argh.ArghParser()
    parser.add_commands([infer_gene_trees,
                         codeml_setup])
    parser.dispatch()


if __name__ == '__main__':
    main()
