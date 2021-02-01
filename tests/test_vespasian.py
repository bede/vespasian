import os
import subprocess

import pytest

from vespasian import vespasian


cwd = os.getcwd()
data_dir = 'tests/data'  # Test from project root

def run(cmd, cwd=cwd):
    return subprocess.run(cmd,
                          cwd=cwd,
                          shell=True,
                          check=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def test_setup_pipeline():
    run('rm -rf gene-trees codeml codeml-unlabelled', cwd=data_dir)  # Start afresh


def test_infer_genetree():  # default --output is gene-trees
    run_cmd = run('vespasian infer-gene-trees frog-families frogs.nwk', cwd=data_dir)
    print(run_cmd.stdout, run_cmd.stderr)


def test_infer_genetree_with_errant_taxa():
    with pytest.raises(NameError) as exec_info:
        vespasian.infer_gene_trees(f'{data_dir}/frog-families-errant-taxa',
                               f'{data_dir}/frogs.nwk',
                               f'{data_dir}/gene-trees-errant-taxa')
    assert 'alloallo' in str(exec_info.value)


def test_infer_genetree_with_duplicate_taxa():
    with pytest.raises(RuntimeError) as exec_info:
        vespasian.infer_gene_trees(f'{data_dir}/duplicate-taxa/aln',
                               f'{data_dir}/duplicate-taxa/saiga.nwk',
                               f'{data_dir}/gene-trees-duplicate-taxa')
    assert 'Duplicate taxa present' in str(exec_info.value)


def test_codeml_setup():  # default --output is codeml
    run_cmd = run('vespasian codeml-setup frog-families gene-trees --branches frogs.yaml', cwd=data_dir)
    print(run_cmd.stdout, run_cmd.stderr)


def test_codeml_setup_unlabelled():
    run_cmd = run('vespasian codeml-setup frog-families gene-trees --output codeml-unlabelled', cwd=data_dir)
    print(run_cmd.stdout, run_cmd.stderr)


def test_report():
    # vespasian.report(f'{data_dir}/report_testing/codeml/1082_7', f'{data_dir}/report_testing/codeml/1082_7')
    run('rm -rf report_testing/codeml_report')
    run_cmd = run('vespasian report report_testing/codeml/1082_7 --output report_testing/codeml_report', cwd=data_dir)


# def test_snakefile_dryrun():
#     run(f'snakemake --dryrun', cwd=f'{data_dir}/codeml')

def test_peter_fusions():
    run_cmds = []
    run_cmds.append(run('rm -rf gene-trees codeml report-codeml',
                    cwd=f'{data_dir}/peter'))
    run_cmds.append(run('vespasian infer-gene-trees aln F6935_Domain1_nuc_regions.nwk',
                        cwd=f'{data_dir}/peter'))
    run_cmds.append(run('vespasian codeml-setup -b branches.yaml aln gene-trees',
                        cwd=f'{data_dir}/peter'))
    run_cmds.append(run('vespasian report snakemaked',
                        cwd=f'{data_dir}/peter'))

def test_vlad_report():
    run_cmds = []
    run_cmds.append(run('rm -rf run_00/report-codeml run_01/report-codeml',
                    cwd=f'{data_dir}/vlad'))
    run_cmds.append(run('vespasian report codeml',
                        cwd=f'{data_dir}/vlad/run_00'))
    run_cmds.append(run('vespasian report codeml',
                        cwd=f'{data_dir}/vlad/run_01'))


def test_report_frog_modela_neb():
    run_cmd = run('vespasian report frog-modela-neb/family --output frog-modela-neb/report-codeml', cwd=data_dir)
    run('rm -rf frog-modela-neb/report-codeml', cwd=data_dir)


# Slow tests, run with --slow


@pytest.mark.slow
# def test_codeml_single_family():
#     sample_workspace = 'codeml/s_4787/s_4787/m3Discrtk3/Omega0'
#     run(f'codeml', cwd=f'{data_dir}/{sample_workspace}')
#     run(f'rm 2NG* rst* rub out', cwd=f'{data_dir}/{sample_workspace}')

@pytest.mark.slow
def test_vlad_3():
    '''
    Catch labelled root. Gymnophiona should not be labelled as it is the root and meaning no background
    This test has been repurposed as the scenario in for which it was conceived *should* no longer be possible
    '''
    sample_workspace = 'vlad/3'
    run('vespasian infer-gene-trees -w PNSAF3/ the_tree_short.nwk', cwd=f'{data_dir}/{sample_workspace}')
    run('vespasian codeml-setup -w -b nodes.yaml PNSAF3/ gene-trees/', cwd=f'{data_dir}/{sample_workspace}')
    assert not os.path.exists(f'{data_dir}/{sample_workspace}/codeml/PNSA_1/PNSA_1_Gymnophiona/')
    run(f'rm -rf gene-trees codeml', cwd=f'{data_dir}/{sample_workspace}')


