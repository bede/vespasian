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


# def setup_pipeline():
#     run_cmd('rm -rf gene-trees codeml codeml-unlabelled', cwd=data_dir)  # Start afresh


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
    with pytest.raises(NameError) as exec_info:
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

# Slow tests, run with --slow


@pytest.mark.slow
def test_codeml_single_family():
    sample_workspace = 'codeml/s_4787/s_4787/m3Discrtk3/Omega0'
    run(f'codeml', cwd=f'{data_dir}/{sample_workspace}')
    run(f'rm 2NG* rst* rub out', cwd=f'{data_dir}/{sample_workspace}')
