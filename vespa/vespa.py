import os
import sys

from pathlib import Path

import treeswift

from Bio import SeqIO



for family in families:
    gene_tree = infer_gene_tree(family)
    for 

# Project dir structure?
#   alignments
#   gene_trees
#   codeml_workspace


def infer_gene_tree(alignment_path, tree_path, output_path):
    '''Build gene tree by pruning species tree to contain only nodes present in alignment'''
    records = SeqIO.parse(alignment_path, 'fasta')
    records_headers = set((r.id.partition('|')[0] for r in records))
    species_tree = treeswift.read_tree_newick(tree_path)
    gene_tree = species_tree.extract_tree_with(records_headers)
    gene_tree.write_tree_newick(output_path)


def infer_gene_trees(input_path, tree_path, output_path):
    '''Build gene trees for a directory containing gene family alignments'''
    family_paths = [f'{input_path}/{fn}' for fn in os.listdir(input_path)
                    if fn.endswith(('.fa', '.fasta'))]
    families_paths = {Path(a).stem: a for a in family_paths}
    for family, path in families_paths.items():
        # family_out_path = str(Path(output_path)/family)
        os.makedirs(output_path, exist_ok=True)
        infer_gene_tree(path, tree_path, f'{output_path}/{family}.nwk')



def codeml_setup_family(alignment_path, output_path):
    '''Create codeml workspaces for suite of site and branch-site models for single gene family'''
    
    class ControlFile():
        '''Create the required config to generate a CodeML control (.ctl) file'''
        
        def __init__(self, model_name, model=0, NSsites=8, fix_omega=1, ncatG=10, omega=1):
            self.seqfile = 'align.phy'
            self.treefile = 'tree'
            self.outfile = 'out'
            self.noisy = 3
            self.verbose = 0
            self.runmode = 0
            self.seqtype = 1
            self.CodonFreq = 2
            self.aaDist = 0
            self.aaRatefile = 'wag.dat'
            self.model = model
            self.NSsites = NSsites
            self.icode = 0
            self.fix_kappa = 0
            self.kappa = 3
            self.fix_omega = fix_omega
            self.omega = omega
            self.fix_alpha = 1
            self.alpha = 0
            self.Malpha = 0
            self.ncatG = ncatG
            self.clock = 0
            self.getSE = 0
            self.RateAncestor = 0
            self.Small_Diff = .5e-6
        
        def __str__(self):
            return str(self.__dict__)

        def write(self, path):
            '''Write codeml .ctl file to specified path'''
            ctl = ''
            for k, v in self.__dict__.items():
                ctl += f'{k} = {v}\n'
            with open(path, 'w+') as ctl_fh:
                ctl_fh.write(ctl)

    models_presets = {
        'm0': {'model': 0, 'NSsites': 0, 'fix_omega': 0, 'ncatG': 1},
        'm1Neutral': {'model': 0, 'NSsites': 1, 'fix_omega': 0, 'ncatG': 2},
        'm2Selection': {'model': 0, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 3},
        'm3Discrtk2': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 2},
        'm3Discrtk3': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 3},
        'm7': {'model': 0, 'NSsites': 7, 'fix_omega': 0, 'ncatG': 10},
        'm8': {'model': 0, 'NSsites': 8, 'fix_omega': 0, 'ncatG': 10},
        'm8a': {'model': 0, 'NSsites': 8, 'fix_omega': 1, 'ncatG': 10, 'omega': 1},
        'modelA': {'model': 2, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 4},
        'modelAnull': {'model': 2, 'NSsites': 2, 'fix_omega': 1, 'ncatG': 4, 'omega': 1},
    }



    records = SeqIO.parse(alignment_path, 'fasta')
    record_ids = {r.id for r in records}  # Set comprehension

    for model, presets in models_presets.items():
        if model in models_with_fixed_omega:  # m8a and modelAnull
            os.makedirs(f'{codeml_path}/{model}/Omega1')
            ControlFile(model, **presets).write(f'{codeml_path}/{model}/Omega1/codeml.ctl')
        else:
            for omega in omega_values:
                presets['omega'] = omega
                os.makedirs(f'{codeml_path}/{model}/Omega{omega}')
                ControlFile(model, **presets).write(f'{codeml_path}/{model}/Omega{omega}/codeml.ctl')


def label_branch():
    pass


def label_branches(input_path, tree_path, output_path):
    species_tree = treeswift.read_tree_newick(tree_path)
    for name, children in branches.items():
        # tree = treeswift.read_tree_newick('(Cten,Tcom,(((((lung,(lati,(shar,(lamp,(taki,dani))))),((plat,(elep,(huma,cows))),((turt,(zebr,chic)),(cobr,anol)))),((hyno,((newt,cyno),(axol,atig))),(xeno,(nano,(scin,(dend,(pris,(rani,(amee,allo))))))))),Rbiv),(Muni,Mder)));')
        # print('name:', name)
        if children:  # Has children
            mrca = tree.mrca(set(children))  # Fetch MRCA
            mrca.label = "'#1'"
            # print(tree.newick(), '\n')
        
        else:  # Childless, assumed to be leaf node?? Check this! Could be internal node?
            pass


def codeml_setup_families():
    pass

















def workspace_setup(family_path):
    '''Create a CodeML workspace for performing standard tests on an orthologous family'''
    
    ctl_defaults = {
        seqfile: 'align.phy', 
        treefile: 'tree', 
        outfile: 'out', 
        noisy: 3,
        verbose: 0,
        runmode: 0,
        seqtype: 1,
        CodonFreq: 2,
        aaDist: 0,
        aaRatefile: 'wag.dat',
        model: 0,
        NSsites: 0,
        icode: 0,
        fix_kappa: 0,
        kappa: 3,
        fix_omega: 0,
        omega: 1,
        fix_alpha: 1,
        alpha: 0,
        Malpha: 0,
        ncatG: 1,
        clock: 0,
        getSE: 0,
        RateAncestor: 0,
        Small_Diff: .5e-6
    }

    models_presets = {
        'm0': {'model': 0, 'NSsites': 0, 'fix_omega': 0, 'ncatG': 1},
        'm1Neutral': {'model': 0, 'NSsites': 1, 'fix_omega': 0, 'ncatG': 2},
        'm2Selection': {'model': 0, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 3},
        'm3Discrtk2': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 2},
        'm3Discrtk3': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 3},
        'm7': {'model': 0, 'NSsites': 7, 'fix_omega': 0, 'ncatG': 10},
        'm8': {'model': 0, 'NSsites': 8, 'fix_omega': 0, 'ncatG': 10},
        'm8a': {'model': 0, 'NSsites': 8, 'fix_omega': 1, 'ncatG': 10, 'omega': 1},
        'modelA': {'model': 2, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 4},
        'modelAnull': {'model': 2, 'NSsites': 2, 'fix_omega': 1, 'ncatG': 4, 'omega': 1},
    }

    models_with_fixed_omega = set('m8a', 'modelAnull')  # Models with fixed omega value of 1

    omega_values = [0, 1, 2, 10]
    
    model_dirs = [f'{family_path}/{m}' for m in models_presets.keys()]

#     for model, presets in models_presets.items():
#         os.mak
#         ctl_config = ctl_defaults**
#         ctl_config.update(presets)
#         if model not in models_with_fixed_omega:
#             for omega in omega_values:
#                     ctl_config['omega'] = omega


# for model, presets in models_presets.items():
#     if model not in models_with_fixed_omega:
#         for omega in omega_values:
#             ctl = ControFile(model, presets*)


