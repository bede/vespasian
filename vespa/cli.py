import sys
import argh

from Bio import SeqIO


def infer_genetree():
    pass

def codeml_setup(input: 'path to directory containing gene families',
                 branch_file: 'file containing branches to label'):
    print('Done!')


def main():
    parser = argh.ArghParser()
    parser.add_commands([infer_genetree,
                         codeml_setup])
    parser.dispatch()


if __name__ == '__main__':
    main()
