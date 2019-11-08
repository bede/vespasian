import os
import re
import shutil
import textwrap
import warnings

import yaml
import tqdm
import parmap
import treeswift

import numpy as np
import pandas as pd

from itertools import permutations
from collections import defaultdict
from pathlib import Path

from Bio import AlignIO, SeqIO
from scipy.stats import chi2

from vespasian import util



def parse_branch_file(branch_file):
    '''Parse YAML file containing foreground lineage definitions for codeml analysis'''
    with open(branch_file, 'r') as stream:
        try:
            branches = yaml.load(stream)
        except yaml.YAMLError:
            print('Problem parsing branch file')
            raise RuntimeError
    return branches


def unroot(tree):
    '''
    Trifurcate bifurcating Treeswift Tree (as required by codeml)
    Could possibly be replaced by TreeSwift's new deroot function?
    '''
    for node in tree.traverse_levelorder():
        if node.is_root():
            if node.num_children() > 2:
                return tree  # Root already trifurcating
            children = node.child_nodes()
            if children[0].num_children() == 2:
                donor = children[0]
            elif children[1].num_children() == 2:
                donor = children[1]
            else:
                raise RuntimeError('Tree assumptions broken')
            donor_children = donor.child_nodes()
            node.add_child(donor_children[0])
            node.add_child(donor_children[1])
            node.remove_child(donor)
            return tree


def infer_gene_tree(alignment_path, tree_path, output_path, separator='|'):
    '''Build gene tree by pruning species tree to contain only nodes present in alignment'''
    records = list(SeqIO.parse(alignment_path, 'fasta'))
    stems_names = {r.id.partition(separator)[0]: r.id for r in records}  # Map prefix to full name
    stems = set(stems_names.keys())
    if len(records) == 0:
        raise RuntimeError(f'No records found in {alignment_path}')
    elif len(records) != len(stems):
        raise RuntimeError(f'Duplicate taxa present in {alignment_path}')
    elif len(records) < 7:
        warnings.warn(f'Fewer than 7 sequences found in {alignment_path}. Consider exclusion.')
    species_tree = treeswift.read_tree_newick(tree_path)
    gene_tree = species_tree.extract_tree_with(stems)
    gene_tree.rename_nodes(stems_names)
    gene_tree.is_rooted = False
    gene_tree.ladderize()
    unroot(gene_tree).write_tree_newick(output_path)
    return stems


def infer_gene_trees(input_path, tree_path, output_path, separator='|', progress=False):
    '''Build gene trees for a directory containing aligned gene families'''
    family_paths = [f'{input_path}/{fn}' for fn in os.listdir(input_path)
                    if fn.endswith(('.fa', '.fasta'))]
    families_paths = {Path(a).stem: a for a in family_paths}
    if len(family_paths) == 0:
        raise RuntimeError(f'No fasta files found in {input_path}')

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


def prune_redundant_branches(branches, tree_path, separator):
    '''
    Infer branch hierarchy and remove redundant supertree branches
    i.e. those with identical subtrees
    Return dictionary of branches and leaves without redundancy 
    '''
    tree = treeswift.read_tree_newick(tree_path)
    taxa_present = {t.label.partition(separator)[0] for t in tree.traverse_leaves()}
    
    branches_leaves = {b: set(l) for b, l in branches.items()} 
    branches_leaves_present = {b: l.intersection(taxa_present) for b, l in branches_leaves.items()}
    
    supertrees_subtrees = {i: j for i, j in permutations(branches_leaves, 2)
                                if branches_leaves[i].issuperset(branches_leaves[j])
                                and branches_leaves_present[i] == branches_leaves_present[j]}

    for supertree, subtree in supertrees_subtrees.items():
        del branches_leaves[supertree]
        warnings.warn(f'Skipped labelling redundant supertree {supertree} for {tree_path}')
        print(f'Skipped labelling redundant supertree {supertree} for {tree_path}')

    return  {b: list(l) for b, l in branches_leaves.items()}  # values() from sets to lists


def label_branch(tree_path, branch_label, leaf_labels, separator='|', strict=False):
    '''
    Return Tree with codeml labelled ancestral branch or leaf node if children absent
    Default behaviour is to require >= 2 specified taxa to be present
    Strict behaviour requires all taxa to be present in tree
    Throws RuntimeError() when a branch should be skipped
    Requires leaf nodes be specified without values
    '''
    tree = treeswift.read_tree_newick(tree_path)
    taxa_longnames = {n.label.partition(separator)[0]: n.label for n in tree.traverse_leaves()}
    taxa = set(taxa_longnames.keys())

    if leaf_labels:  # Branch is internal node
        tree_leaf_intersection = set(taxa_longnames.keys()).intersection(leaf_labels)  # Coerces to set
        intersection_size = len(tree_leaf_intersection)
        required_taxa = len(leaf_labels) if strict else 2  # Ignore singletons by default
        
        if intersection_size < required_taxa: 
            warnings.warn(f'Insufficient taxa present to label branch {branch_label} in {tree_path}')
            print(f'Insufficient taxa present to label branch {branch_label} in {tree_path}')
            raise RuntimeError()

        expanded_leaf_labels = (taxa_longnames.get(l) for l in leaf_labels)  # Generate long names
        mrca = tree.mrca(set(filter(None, expanded_leaf_labels)))
        mrca.label = "'#1'"

    else:  # Branch is leaf node
        if branch_label not in taxa:
            warnings.warn(f'Insufficient taxa present to label leaf node {branch_label} in {tree_path}')
            print(f'Insufficient taxa present to label leaf node {branch_label} in {tree_path}')
            raise RuntimeError()
        labels_nodes = tree.label_to_node()
        labels_nodes[taxa_longnames[branch_label]].label += '#1'

    if '#1' not in tree.newick():  # Catch-all for unlabelled branches
        warnings.warn(f'Skipped labelling {branch_label} in {tree_path}')
        print(f'Skipped labelling {branch_label} in {tree_path}')
        raise RuntimeError()
    return tree


####################################################################################################


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

    shutil.copy(alignment_path, f'{family_path}/align.fa')
    util.convert_phylip(f'{family_path}/align.fa', f'{family_path}/align.phy')
    shutil.copy(gene_tree_path, f'{family_path}/tree.nwk')

    for model, codeml_params in models.items():
        for omega in models_omega[model]:
            codeml_params['omega'] = omega
            run_dir = f'{family_path}/{family_name}/{model}/Omega{omega}'
            os.makedirs(run_dir, exist_ok=True)
            if not os.path.exists(f'{run_dir}/tree.nwk'):
                os.symlink('../../../tree.nwk', f'{run_dir}/tree.nwk')  # Relative symlink
                # shutil.copy(f'{family_path}/tree.nwk', f'{run_dir}/tree.nwk')
            if not os.path.exists(f'{run_dir}/align.fa'):
                os.symlink('../../../align.fa', f'{run_dir}/align.fa')  # Relative symlink
                # shutil.copy(f'{family_path}/align.fa', f'{run_dir}/align.fa')
            ControlFile(model, **codeml_params).write(f'{run_dir}/codeml.ctl')


def setup_branch_site_models(family_name, family_path, alignment_path, gene_tree_path, branches,
                             separator, strict):
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
        childless_branches = {b: l for b, l in branches.items() if not l}
        childbearing_branches = {b: l for b, l in branches.items() if l}
        childbearing_branches_pruned = prune_redundant_branches(childbearing_branches,
                                                                f'{family_path}/tree.nwk',
                                                                separator)
        
        branches_pruned = {**childless_branches, **childbearing_branches_pruned}
        
        for branch, leaves in branches_pruned.items():
            branch_path = f'{family_path}/{family_name}_{branch}'
            try:
                labelled_tree = label_branch(f'{family_path}/tree.nwk',
                                             branch, leaves, separator, strict)
            except RuntimeError:  # label_branch() throws RuntimeError to skip a families
                continue
            for model, codeml_params in models.items():
                for omega in models_omega[model]:
                    codeml_params['omega'] = omega
                    run_dir = f'{branch_path}/{model}/Omega{omega}'
                    os.makedirs(run_dir, exist_ok=True)
                    labelled_tree.write_tree_newick(f'{run_dir}/tree.nwk')
                    if not os.path.exists(f'{run_dir}/align.fa'):
                        os.symlink('../../../align.fa', f'{run_dir}/align.fa')  # Relative symlink
                        # shutil.copy(f'{family_path}/align.fa', f'{run_dir}/align.fa')
                    ControlFile(model, **codeml_params).write(f'{run_dir}/codeml.ctl')


def gather_codeml_dirs(path):
    '''Recursively finds codeml directory paths within a target directory tree'''
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=True):
            if 'Omega' in entry.name:
                yield entry.path
            else:
                yield from gather_codeml_dirs(entry.path)


def gather_codeml_commands(path):
    '''Generates list of codeml commands to be executed'''
    codeml_dirs = gather_codeml_dirs(path)
    codeml_dirs_no_prefix = [Path(p).relative_to(path) for p in codeml_dirs]  # Relative to codeml/
    return [f'cd {c} && codeml' for c in codeml_dirs_no_prefix]


def generate_codeml_commands(output_dir):
    '''Gather and write codeml commands to file'''
    cmds = gather_codeml_commands(output_dir)
    with open(f'{output_dir}/codeml-commands.sh', 'w+') as cmds_fh:
        cmds_fh.write('\n'.join(cmds))


def generate_snakefile(output_dir):
    '''Write Snakemake Snakefile inside codeml root directory'''
    snakefile =  '''
        with open('codeml-commands.sh', 'r') as fh:
            commands = [l.strip() for l in fh if l]
            run_dirs = [c[3:].partition(' &&')[0] for c in commands]

        rule all:
            input:
                expand('{od}/out', od=run_dirs)

        rule run:
            output:
                '{od}/out'
            shell:
                'cd {wildcards.od} && codeml'
    '''
    with open(f'{output_dir}/Snakefile', 'w+') as snakefile_fh:
        snakefile_fh.write(textwrap.dedent(snakefile).strip())


def setup_family(family, alignment_path, output_dir, gene_trees_dir, branches, separator, strict):
    '''Parallelisable family level setup function'''
    family_path = f'{output_dir}/{family}'
    os.makedirs(f'{family_path}', exist_ok=True)
    gene_tree_path = f'{gene_trees_dir}/{family}.nwk'
    setup_site_models(family, family_path, alignment_path, gene_tree_path)
    setup_branch_site_models(family, family_path, alignment_path, gene_tree_path, branches,
                             separator, strict)


def codeml_setup(families_dir, gene_trees_dir, branch_file, output_dir, separator, strict,
                 threads=os.cpu_count(), progress=False):
    '''Configure site and branch-site test environment given alignments, gene trees and branches'''
    alignment_paths = [f'{families_dir}/{fn}' for fn in os.listdir(families_dir)
                    if fn.endswith(('.fa', '.fasta'))]
    alignments_paths = {Path(a).stem: a for a in alignment_paths}
    branches = parse_branch_file(branch_file) if branch_file else None

    # Was a nightmare to debug until pm_parallel logic was added to NOT use multiprocessing.Pool
    # When `--threads 1`, allowing switch to run in single process
    parmap.starmap(setup_family,
                   alignments_paths.items(),
                   output_dir,
                   gene_trees_dir,
                   branches,
                   separator,
                   strict,
                   pm_pbar=progress,
                   pm_parallel=False if threads == 1 else True,
                   pm_processes=threads,
                   pm_chunksize=10)

    generate_codeml_commands(output_dir)
    generate_snakefile(output_dir)


####################################################################################################


def gather_codeml_output(path):
    '''Recursively finds codeml output paths within a target directory tree'''
    for entry in os.scandir(path):
        if entry.is_file(follow_symlinks=True):
            if entry.name == 'out':
                yield entry.path
        else:
            yield from gather_codeml_output(entry.path)


def parse_result(path):
    '''Parse codeml output'''
    record = path.split('/')[-5:-1]
    
    result = dict(family=record[0],
                  tree=record[1],
                  model=record[2],
                  w_t0=int(''.join(i for i in record[3] if i.isdigit())),
                  path=path,
                  neb_sites=[],
                  beb_sites=[])
   
    model = result['model']
    
    params = {}
    
    floats_re = re.compile(r"[-+]?\d*\.\d+|[-+]?\d+")
    
    try:
        with open(path) as result_fh:
            result_contents = result_fh.read()
        result_lines = result_contents.splitlines()
        for line in result_lines:
            if line.startswith('lnL'):
                lnl_record = re.split(r'\s+', line.strip())  # Split by contiguous whitespace
                result['lnl'] = float(lnl_record[-2])
                assert result['lnl'] <= 0  # lnL should be negative
            elif model == 'm7':
                if line.startswith(' p ='):
                    m7_pq_record = floats_re.findall(line)
                    params['p'], params['q'] = m7_pq_record
            elif model == 'm8' or model == 'm8a':
                if line.startswith('  p0 ='):
                    m8_p0pq_record = floats_re.findall(line)
                    params['p0'], params['p'], params['q'] = m8_p0pq_record[1:]
                elif line.startswith(' (p1 ='):
                    m8_p1w_record = floats_re.findall(line)
                    params['p1'], params['w'] = m8_p1w_record[1:]
            elif model == 'modelA' or model == 'modelAnull':
                if line.startswith('proportion'):
                    p_record = floats_re.findall(line)
                    p_labels = tuple(range(len(p_record)))
                    p_params = tuple(map(float, p_record))
                    p_labels_params = {f'p{p_label}': p_param for p_label, p_param
                                       in zip(p_labels, p_params)}
                    params = {**params, **p_labels_params}
                elif line.startswith('foreground w'):
                    w_record = floats_re.findall(line)
                    w_labels = tuple(range(len(w_record)))
                    w_params = tuple(map(float, w_record[:3]))
                    w_labels_params = {f'w{w_label}': w_param for w_label, w_param
                                       in zip(w_labels, w_params)}
                    params = {**params, **w_labels_params}
            else:
                if line.startswith('omega'):
                    omega_record = floats_re.findall(line)
                    params['w'] = float(omega_record[0])
                elif line.startswith('p:'):
                    p_record = floats_re.findall(line)
                    p_labels = tuple(range(len(p_record)))
                    p_params = {f'p{p_label}': p_param for p_label, p_param
                                in zip(p_labels, tuple(map(float, p_record)))}
                    params = {**params, **p_params}
                elif line.startswith('w:'):
                    w_record = floats_re.findall(line)
                    w_labels = tuple(range(len(w_record)))
                    w_params = {f'w{w_label}': w_param for w_label, w_param
                                in zip(w_labels, tuple(map(float, w_record)))}
                    params = {**params, **w_params}

        # Get NEB and BEB positive sites for the site models where they are reported
        if model in ('m0','m2Selection', 'm3Discrtk2', 'm3Discrtk3', 'm8'):
            neb_pos_site_lines = re.findall(r'Naive Empirical Bayes \(NEB\).*?mean \+\- SE for w(.*?)\n\n\n', result_contents, re.DOTALL)
            beb_pos_site_lines = re.findall(r'Bayes Empirical Bayes \(BEB\).*?mean \+\- SE for w(.*?)\n\n\n', result_contents, re.DOTALL)
            if neb_pos_site_lines:
                neb_lines = neb_pos_site_lines[0].strip().replace('*','').split('\n')
                neb_lines = list(filter(None, neb_lines))  # Cull empty strings
                # print(path, 'site_neb', neb_lines)
                neb_records = [{'position': int(r[0]), 'residue': r[1], 'p': float(r[2])}
                               for r in [s.split() for s in neb_lines]]
                result['neb_sites'] = neb_records
            if beb_pos_site_lines:
                beb_lines = beb_pos_site_lines[0].strip().replace('*','').split('\n')
                beb_lines = list(filter(None, beb_lines))  # Cull empty strings
                # print(path, 'site_beb', beb_lines)
                beb_records = [{'position': int(r[0]), 'residue': r[1], 'p': float(r[2])}
                               for r in [s.split() for s in beb_lines]]
                result['beb_sites'] = beb_records
            # Hack to prevent overread with confusing empty NEB immediately followed by non-empty BEB
            # A non-greedy shortest-of-either three newlines or seeing 'Bayes' to terminate the regex search would be better
            if '(NEB) analysis\nBayes Empirical Bayes (BEB)' in result_contents:
                result['neb_sites'] = []

        # Get NEB and BEB positive sites for branch-site modelA
        if model == 'modelA':
            # Match records after non-greedy ('?') expansion of Naive|Bayes Empirical Bayes * Prob(w>1):
            neb_pos_site_lines = re.findall(r'Naive Empirical Bayes \(NEB\).*?Prob\(w\>1\)\:\n(.*?)\n\n', result_contents, re.DOTALL)
            beb_pos_site_lines = re.findall(r'Bayes Empirical Bayes \(BEB\).*?Prob\(w\>1\)\:\n(.*?)\n\n', result_contents, re.DOTALL)
            if neb_pos_site_lines:
                neb_lines = neb_pos_site_lines[0].strip().replace('*','').split('\n')
                neb_lines = list(filter(None, neb_lines))  # Cull empty strings
                # print(path, 'branch_site_neb', neb_lines)
                neb_records = [{'position': int(r[0]), 'residue': r[1], 'p': float(r[2])}
                               for r in [s.strip('\n').split() for s in neb_lines]]
                result['neb_sites'] = neb_records
            if beb_pos_site_lines:
                beb_lines = beb_pos_site_lines[0].strip().replace('*','').split('\n')
                beb_lines = list(filter(None, beb_lines))  # Cull empty strings
                # print(path, 'branch_site_beb', beb_lines)
                beb_records = [{'position': int(r[0]), 'residue': r[1], 'p': float(r[2])}
                               for r in [s.split() for s in beb_lines]]
                result['beb_sites'] = beb_records


    except Exception as e:
        print(f'Problem parsing codeml output {path}')
        raise(e)

    result['params'] = params
    return result


def filter_results(results):
    '''Returns the family-tree-model results with the highest log likelihood'''
    names_lnls = defaultdict(lambda: float('-inf'))
    names_records = {}
    for r in results:
        # We only want the best of each family-tree-model permutation
        # name = f"{r['family']}_{r['tree']}_{r['model']}"
        key = (r['family'], r['tree'], r['model'])  # Tuple-keyed dict of family-tree-model
        if r['lnl'] > names_lnls[key]:
            names_lnls[key] = r['lnl']
            names_records[key] = r
    return names_records


def parse_results(input_dir):
    output_paths = list(gather_codeml_output(input_dir))
    if len(output_paths) == 0:
        print('No codeml output files detected')
    results = [parse_result(path) for path in output_paths]
    filtered_results = filter_results(results)
    return filtered_results


def summarise(family_results):
    def fmt_sites(sites):
        sites_fmt = []
        if sites:
            for s in sites:
                sites_fmt.append(f"{s['position']}{s['residue']}:{s['p']}")
        return ' '.join(sites_fmt)

    df = pd.DataFrame.from_dict(family_results.values()).fillna(0)
    df.sort_values('model', axis=0, inplace=True)
    df['params'] = df['params'].apply(lambda p: ' '.join([f'{k}={v}' for k, v in p.items()]))
    df['neb_sites'] = df['neb_sites'].apply(fmt_sites)
    df['beb_sites'] = df['beb_sites'].apply(fmt_sites)
    df['positive_selection'] = np.where(df['beb_sites'], 'yes', 'no')
    df = df.reindex(['family', 'tree', 'model', 'w_t0', 'lnl', 'params', 'positive_selection',
                     'neb_sites', 'beb_sites', 'path'], axis=1)
    return df


def test_likelihood_ratios(family_results):
    '''Docstring'''
    def lr_prob(likelihood_ratio, df):
        return chi2.sf(likelihood_ratio, df)

    def modelanull_vs_modela_prob(likelihood_ratio):
        return lr_prob(likelihood_ratio)/2

    site_lrts = {('m0', 'm3Discrtk2'): {'df': 2, 'cv': 5.991},
                 ('m1Neutral', 'm2Selection'): {'df': 2, 'cv': 5.991},
                 ('m3Discrtk2', 'm3Discrtk3'): {'df': 1, 'cv': 1.000},
                 ('m7', 'm8'): {'df': 2, 'cv': 5.991},
                 ('m8', 'm8a'): {'df': 2, 'cv': 2.706}}

    branch_site_lrts = {('m1Neutral', 'modelA'): {'df': 2, 'cv': 5.991},
                        ('modelAnull', 'modelA'): {'df': 2, 'cv': 3.841}}
    
    families_branches = {(r[:2]) for r in family_results.keys()}

    lrts = []
    for family, tree in families_branches:
        if family == tree:  # Site models
            for test, meta in site_lrts.items():
                null_lnl = family_results[(family, tree, test[0])]['lnl']
                alt_lnl = family_results[(family, tree, test[1])]['lnl']
                lrt = abs(meta['df']*(null_lnl-alt_lnl))
                lrt_record = dict(tree=tree,
                      lrt=f'{test[0]} vs. {test[1]}',
                      null_model_lnl=null_lnl,
                      alt_model_lnl=alt_lnl,
                      lrt_result=lrt,
                      critical_value=meta['cv'],
                      p=lr_prob(lrt, 1),
                      null_rejected=True if lrt > meta['cv'] else False)
                lrts.append(lrt_record)

        else:  # Branch-site models
            for test, meta in branch_site_lrts.items():
                if test[0] == 'm1Neutral':  # m1Neutral is a site model
                    null_lnl = family_results[(family, family, test[0])]['lnl']
                else:
                    null_lnl = family_results[(family, tree, test[0])]['lnl']
                alt_lnl = family_results[(family, tree, test[1])]['lnl']
                lrt = abs(meta['df']*(null_lnl-alt_lnl))
                lrt_record = dict(tree=tree,
                                  lrt=f'{test[0]} vs. {test[1]}',
                                  null_model_lnl=null_lnl,
                                  alt_model_lnl=alt_lnl,
                                  lrt_result=lrt,
                                  critical_value=meta['cv'],
                                  p=lr_prob(lrt, 1),
                                  null_rejected=True if lrt > meta['cv'] else False)
                lrts.append(lrt_record)

    df = pd.DataFrame(lrts)
    df = df.reindex(('tree', 'lrt', 'null_model_lnl', 'alt_model_lnl', 'lrt_result', 'p',
                     'critical_value', 'null_rejected'), axis=1)
    df = df.sort_values(['tree', 'lrt'])
    return df


def report(input_dir, output_dir):
    '''Perform likelihood ratio tests and and report positively selected sites'''
    filtered_results = parse_results(input_dir)    
    families_results = {f: {} for f in (k[0] for k in filtered_results.keys())}
    for name, result in filtered_results.items():
        families_results[name[0]][name] = result
    
    for family, family_results in families_results.items():
        summary_df = summarise(family_results)
        lrts_df = test_likelihood_ratios(family_results)
        os.makedirs(f'{output_dir}/{family}', exist_ok=True)
        summary_df.to_csv(f'{output_dir}/{family}/summary.tsv', sep='\t', index=False)
        lrts_df.to_csv(f'{output_dir}/{family}/lrts.tsv', sep='\t', index=False)

    return families_results
