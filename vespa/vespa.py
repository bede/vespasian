import os
import shutil
import warnings

import yaml
import tqdm
import treeswift


from pathlib import Path

from Bio import AlignIO, SeqIO

from vespa import util



def parse_branch_file(branch_file):
    '''Parse YAML file containing foreground lineage definitions for codeml analysis'''
    with open(branch_file, 'r') as stream:
        try:
            branches = yaml.load(stream)
        except yaml.YAMLError:
            print('Problem parsing branch file')
            raise
    return branches


def infer_gene_tree(alignment_path, tree_path, output_path, separator='|'):
    '''Build gene tree by pruning species tree to contain only nodes present in alignment'''
    records = list(SeqIO.parse(alignment_path, 'fasta'))
    stems_names = {r.id.partition(separator)[0]: r.id for r in records}  # Map prefix to full name
    stems = set(stems_names.keys())
    species_tree = treeswift.read_tree_newick(tree_path)
    gene_tree = species_tree.extract_tree_with(stems)
    gene_tree.rename_nodes(stems_names)
    gene_tree.write_tree_newick(output_path)
    return stems


def infer_gene_trees(input_path, tree_path, output_path, separator='|', progress=False):
    '''Build gene trees for a directory containing aligned gene families'''
    family_paths = [f'{input_path}/{fn}' for fn in os.listdir(input_path)
                    if fn.endswith(('.fa', '.fasta'))]
    families_paths = {Path(a).stem: a for a in family_paths}
    
    observed_stems = set()
    for family, path in tqdm.tqdm(families_paths.items(), disable=not progress):
        os.makedirs(output_path, exist_ok=True)
        stems = infer_gene_tree(path, tree_path, f'{output_path}/{family}.nwk')
        observed_stems |= stems

    # Ensure species tree contains all taxa observed within alignments
    species_tree_labels = treeswift.read_tree_newick(tree_path).labels(internal=False)
    species_tree_stems = set(l.partition(separator)[0] for l in species_tree_labels)
    if observed_stems.difference(species_tree_stems):
        raise NameError(f'The following taxa were not found in the species tree: \n'
                        f'{observed_stems.difference(species_tree_stems)}')




def label_branch(tree_path, branch_label, leaf_labels, separator='|'):
    '''Return Tree with codeml labelled ancestral branch or leaf node if children absent'''
    tree = treeswift.read_tree_newick(tree_path)
    stems_names = {n.label.partition(separator)[0]: n.label for n in tree.traverse_leaves()}
    if leaf_labels:  # internal node
        expanded_leaf_labels = (stems_names.get(l) for l in leaf_labels)  # generate long names
        mrca = tree.mrca(set(filter(None, expanded_leaf_labels)))  # fetch MRCA else throws RunTimeError
        mrca.label = "'#1'"
    else:  # leaf node
        labels_nodes = tree.label_to_node()
        if branch_label in labels_nodes:  # leaf present
            labels_nodes[branch_label].label += '#1'
    return tree


class ControlFile():
    '''Create the required config to generate a CodeML control (.ctl) file'''
    def __init__(self, model_name, model=0, NSsites=8, fix_omega=1, ncatG=10, omega=1):
        self.seqfile = 'align.fa'
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


def setup_site_models(family_name, family_path, alignment_path, gene_tree_path):
    '''Configure codeml site model environments'''
    models = {
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
        'm8a': [1]
    }

    alignment = list(AlignIO.parse(alignment_path, 'fasta'))
    AlignIO.write(alignment[0], f'{family_path}/align.fa', 'fasta')
    util.convert_phylip(f'{family_path}/align.fa', f'{family_path}/align.phy')
    shutil.copy(gene_tree_path, f'{family_path}/tree.nwk')

    for model, params in models.items():
        for omega in models_omega[model]:
            params['omega'] = omega
            run_dir = f'{family_path}/{family_name}/{model}/Omega{omega}'
            os.makedirs(run_dir, exist_ok=True)
            shutil.copy(f'{family_path}/tree.nwk', run_dir)
            shutil.copy(f'{family_path}/align.fa', f'{run_dir}/align.fa')
            ControlFile(model, **params).write(f'{run_dir}/codeml.ctl')


def setup_branch_site_models(family_name, family_path, alignment_path, gene_tree_path, branches):
    '''Configure codeml branch-site model environments'''
    models = {
        'modelA': {'model': 2, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 4},
        'modelAnull': {'model': 2, 'NSsites': 2, 'fix_omega': 1, 'ncatG': 4, 'omega': 1},
    }

    models_omega = {
        'modelA': [0, 1, 2, 10],
        'modelAnull': [1]
    }

    alignment = list(AlignIO.parse(alignment_path, 'fasta'))
    AlignIO.write(alignment[0], f'{family_path}/align.fa', 'fasta')
    shutil.copy(gene_tree_path, f'{family_path}/tree.nwk')

    if branches:
        for branch, leaves in branches.items():
            branch_path = f'{family_path}/{family_name}_{branch}'
            try:
                labelled_tree = label_branch(f'{family_path}/tree.nwk', branch, leaves)
            except RuntimeError:  # RuntimeError: no MRCA for these nodes
                continue
            for model, params in models.items():
                for omega in models_omega[model]:
                    params['omega'] = omega
                    run_dir = f'{branch_path}/{model}/Omega{omega}'
                    os.makedirs(run_dir, exist_ok=True)
                    labelled_tree.write_tree_newick(f'{run_dir}/tree.nwk')
                    shutil.copy(f'{family_path}/align.fa', f'{run_dir}/align.fa')
                    ControlFile(model, **params).write(f'{run_dir}/codeml.ctl')


def list_codeml_dirs(path):
    '''Recursively finds codeml directory paths within a target directory tree'''
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            if 'Omega' in entry.name:
                yield entry.path
            else:
                yield from list_codeml_dirs(entry.path)


def list_codeml_commands(path, codeml_binary_path):
    '''Generates list of codeml commands to be executed'''
    codeml_dirs = list_codeml_dirs(path)
    codeml_dirs_no_prefix = [Path(p).relative_to(path) for p in codeml_dirs]  # Relative to codeml/
    return [f'cd {c} && {codeml_binary_path}' for c in codeml_dirs_no_prefix]


def codeml_setup(families_dir, gene_trees_dir, branch_file, output_dir, separator, progress=False):
    '''Configure site and branch-site test environment given alignments, gene trees and branches'''
    alignment_paths = [f'{families_dir}/{fn}' for fn in os.listdir(families_dir)
                    if fn.endswith(('.fa', '.fasta'))]
    alignments_paths = {Path(a).stem: a for a in alignment_paths}
    branches = parse_branch_file(branch_file) if branch_file else None

    for family, alignment_path in tqdm.tqdm(alignments_paths.items(), disable=not progress):
        family_path = f'{output_dir}/{family}'
        os.makedirs(f'{family_path}', exist_ok=True)
        gene_tree_path = f'{gene_trees_dir}/{family}.nwk'
        setup_site_models(family, family_path, alignment_path, gene_tree_path)
        setup_branch_site_models(family, family_path, alignment_path, gene_tree_path, branches)
    
    cmds = list_codeml_commands(output_dir, 'codeml')  # hack
    with open(f'{output_dir}/codeml-commands.sh', 'w+') as cmds_fh:
        cmds_fh.write('\n'.join(cmds))
