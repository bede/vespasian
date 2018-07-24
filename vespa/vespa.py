import os
import sys
import yaml

from pathlib import Path

import treeswift

from Bio import AlignIO, SeqIO



def infer_gene_tree(alignment_path, tree_path, output_path):
    '''Build gene tree by pruning species tree to contain only nodes present in alignment'''
    records = SeqIO.parse(alignment_path, 'fasta')
    records_headers = set((r.id.partition('|')[0] for r in records))
    species_tree = treeswift.read_tree_newick(tree_path)
    gene_tree = species_tree.extract_tree_with(records_headers)
    gene_tree.write_tree_newick(output_path)


def infer_gene_trees(input_path, tree_path, output_path):
    '''Build gene trees for a directory containing aligned gene families'''
    family_paths = [f'{input_path}/{fn}' for fn in os.listdir(input_path)
                    if fn.endswith(('.fa', '.fasta'))]
    families_paths = {Path(a).stem: a for a in family_paths}
    for family, path in families_paths.items():
        os.makedirs(output_path, exist_ok=True)
        infer_gene_tree(path, tree_path, f'{output_path}/{family}.nwk')


def parse_branch_file(branch_file):
    '''Parse YAML file containing foreground lineage definitions for codeml analysis'''
    with open(branch_file, 'r') as stream:
        try:
            branches = yaml.load(stream)
        except yaml.YAMLError:
            print('Problem parsing branch file')
            raise
    return branches



def label_branch(gene_tree_path, ):
    '''Returns modified treeswift.Tree with labelled MRCA of specified child nodes'''
    gene_tree = treeswift.read_tree_newick(gene_tree_path)

    return True


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




def codeml_setup(families_dir, gene_trees_dir, branch_file, output_dir):
    '''Bringing everything together'''
    family_paths = [f'{families_dir}/{fn}' for fn in os.listdir(families_dir)
                    if fn.endswith(('.fa', '.fasta'))]
    families_paths = {Path(a).stem: a for a in family_paths}

    for family, path in families_paths.items():
        suite_path = f'{output_dir}/{family}'
        alignment  = AlignIO.parse(path, 'fasta')
        os.makedirs(f'{suite_path}', exist_ok=True)
        setup_site_models(family, path, suite_path)
    
    # branches = parse_branch_file(branch_file)
    # if branches:
    #     for branch, leaves in branches.items():
    #         for label




def setup_site_models(family_name, alignment_path, output_path):
    '''Create codeml workspaces for site models'''
    
    site_models = {
        'm0': {'model': 0, 'NSsites': 0, 'fix_omega': 0, 'ncatG': 1},
        'm1Neutral': {'model': 0, 'NSsites': 1, 'fix_omega': 0, 'ncatG': 2},
        'm2Selection': {'model': 0, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 3},
        'm3Discrtk2': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 2},
        'm3Discrtk3': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 3},
        'm7': {'model': 0, 'NSsites': 7, 'fix_omega': 0, 'ncatG': 10},
        'm8': {'model': 0, 'NSsites': 8, 'fix_omega': 0, 'ncatG': 10},
        'm8a': {'model': 0, 'NSsites': 8, 'fix_omega': 1, 'ncatG': 10}
    }

    models_omega = {
        'm0': [0, 1, 2, 10],
        'm1Neutral': [0, 1, 2, 10],
        'm2Selection': [0, 1, 2, 10],
        'm3Discrtk2': [0, 1, 2, 10],
        'm3Discrtk3': [0, 1, 2, 10],
        'm7': [0, 1, 2, 10],
        'm8': [0, 1, 2, 10],
        'm8a': [1],
        'modelA': [0, 1, 2, 10],
        'modelAnull': [1]
    }

    records = SeqIO.parse(alignment_path, 'fasta')
    record_ids = {r.id for r in records}  # Set comprehension

    for model, params in site_models.items():
        for omega in models_omega[model]:
            params['omega'] = omega
            os.makedirs(f'{output_path}/{family_name}/{model}/Omega{omega}', exist_ok=True)
            ControlFile(model, **params).write(f'{output_path}/{family_name}/{model}/Omega{omega}/codeml.ctl')



# Project dir structure?
#   alignments
#   gene_trees
#   codeml_workspace


def write_phylip(alignment):
    pass

def command_list():
    pass


def fasta_to_phylip(input_path, output_path):
    alignment = AlignIO.parse(input_path, 'fasta')
    AlignIO.write(alignment, output_path, 'fasta')







def codeml_setup_family(alignment_path, output_path):
    '''Create codeml workspaces for suite of site and branch-site models for single gene family'''
    
    class ControlFile():
        '''Create the required config to generate a CodeML control (.ctl) file'''
        
        def __init__(self, model_name, model=0, NSsites=8, fix_omega=1, ncatG=10, omega=1):
            self.seqfile = 'align.phy'
            self.treefile = 'tree.nwk'
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

    site_models = {
        'm0': {'model': 0, 'NSsites': 0, 'fix_omega': 0, 'ncatG': 1},
        'm1Neutral': {'model': 0, 'NSsites': 1, 'fix_omega': 0, 'ncatG': 2},
        'm2Selection': {'model': 0, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 3},
        'm3Discrtk2': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 2},
        'm3Discrtk3': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 3},
        'm7': {'model': 0, 'NSsites': 7, 'fix_omega': 0, 'ncatG': 10},
        'm8': {'model': 0, 'NSsites': 8, 'fix_omega': 0, 'ncatG': 10},
        'm8a': {'model': 0, 'NSsites': 8, 'fix_omega': 1, 'ncatG': 10, 'omega': 1},
    }

    branch_site_models = {
        'modelA': {'model': 2, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 4},
        'modelAnull': {'model': 2, 'NSsites': 2, 'fix_omega': 1, 'ncatG': 4, 'omega': 1},
    }

    records = SeqIO.parse(alignment_path, 'fasta')
    record_ids = {r.id for r in records}  # Set comprehension

    models_with_fixed_omega = set('m8a', 'modelAnull')  # Models with fixed omega value of 1
    omega_values = [0, 1, 2, 10]

    def setup_branch_site_models():
        pass


    for model, presets in models_presets.items():
        if model in models_with_fixed_omega:  # m8a and modelAnull
            os.makedirs(f'{codeml_path}/{model}/Omega1')
            ControlFile(model, **presets).write(f'{codeml_path}/{model}/Omega1/codeml.ctl')
        else:
            for omega in omega_values:
                presets['omega'] = omega
                os.makedirs(f'{codeml_path}/{model}/Omega{omega}')
                ControlFile(model, **presets).write(f'{codeml_path}/{model}/Omega{omega}/codeml.ctl')



    models_with_fixed_omega = set('m8a', 'modelAnull')  # Models with fixed omega value of 1

    omega_values = [0, 1, 2, 10]




class ControlFile():
    '''Create the required config to generate a CodeML control (.ctl) file'''
    
    def __init__(self, model_name, model=0, NSsites=8, fix_omega=1, ncatG=10, omega=1):
        self.seqfile = 'align.phy'
        self.treefile = 'tree.nwk'
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


