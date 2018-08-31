import warnings

import argh

from vespa import vespa, util



def configure_warnings(show_warnings):
    '''Show or suppress warnings, mainly for TreeSwift Tree.mrca() operations'''
    if show_warnings:
        warnings.filterwarnings('always')
    else:
        warnings.filterwarnings('ignore')


def infer_gene_trees(input: 'path to directory containing gene families',
                     tree: 'path to newick formatted species tree',
                     output: 'path to output directory' = 'gene-trees',
                     separator: 'character separating taxon name and identifier(s)' = '|',
                     warnings: 'show warnings' = False,
                     progress: 'show progress bar' = False):
    '''Create gene trees by pruning a given species tree'''
    configure_warnings(warnings)
    vespa.infer_gene_trees(input, tree, output, separator, progress)


def codeml_setup(input: 'path to directory containing aligned gene families',
                 gene_trees: 'path to directory containing gene trees',
                 branches: 'path to yaml file containing branches to be labelled' = None,
                 output: 'path to output directory' = 'codeml',
                 separator: 'character separating taxon name and identifier(s)' = '|',
                 warnings: 'show warnings' = False,
                 progress: 'show progress bar' = False):
    '''Create suite of branch and branch-site codeml environments'''
    configure_warnings(warnings)
    vespa.codeml_setup(input, gene_trees, branches, output, separator, progress)



def reformat_environments(input: 'path to directory containing codeml environments'):
    '''Reformat vespa-slim codeml environments for use with legacy codeml_reader'''
    util.reformat_environments(input)


def main():
    argh.dispatch_commands([infer_gene_trees,
                            codeml_setup,
                            reformat_environments])


if __name__ == '__main__':
    main()
