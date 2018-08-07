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
                     output: 'path to output directory' = 'gene-trees',
                     separator: 'character separating taxon name and identifier(s)' = '|',
                     progress: 'show progress bar' = False):
    '''CLI: create gene trees by pruning a given species tree'''
    vespa.infer_gene_trees(input, tree, output, separator, progress)


def codeml_setup(input: 'path to directory containing aligned gene families',
                 gene_trees: 'path to directory containing gene trees',
                 branches: 'path to yaml file containing branches to be labelled' = None,
                 output: 'path to output directory' = 'codeml',
                 separator: 'character separating taxon name and identifier(s)' = '|',
                 warnings: 'show warnings' = False,
                 progress: 'show progress bar' = False):
    '''CLI: create suite of branch and branch-site codeml environments'''
    configure_warnings(warnings)
    vespa.codeml_setup(input, gene_trees, branches, output, separator, progress)



def main():
    parser = argh.ArghParser()
    parser.add_commands([infer_gene_trees,
                         codeml_setup])
    parser.dispatch()


if __name__ == '__main__':
    main()
