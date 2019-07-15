import os
import shutil

import tqdm

from Bio import AlignIO, SeqIO

from vespasian import vespasian


def convert_phylip(fasta_path, phylip_path):
    '''Convert fasta to sequential phylip format'''
    chunkify_string = lambda x, n: [x[i:i+n] for i in range(0, len(x), n)]
    records = list(SeqIO.parse(fasta_path, 'fasta'))
    indent_len = max(len(r.id) for r in records) + 2
    indent_fmt = '\n' + ' ' * indent_len
    phylip = f' {len(records)} {len(records[0].seq)}'
    
    for r in records:
        chunks = chunkify_string(str(r.seq), 60)
        phylip += f'\n{r.id:<{indent_len}}{chunks[0]}'
        phylip += indent_fmt + indent_fmt.join(chunks)
    
    with open(phylip_path, 'w+') as phylip_fh:
        phylip_fh.write(phylip)



def reformat_environments(codeml_root):
    '''Reformat vespasian codeml environments for use with legacy codeml_reader'''
    codeml_envs = vespa.list_codeml_dirs(codeml_root)
    for env in codeml_envs:
        convert_phylip(f'{env}/align.fa', f'{env}/align.phy')
        shutil.copy(f'{env}/tree.nwk', f'{env}/tree')
