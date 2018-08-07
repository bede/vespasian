import os
import subprocess


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


def setup_pipeline():
    run_cmd('rm -rf gene-trees codeml codeml-unlabelled', cwd=data_dir)  # Start afresh


def test_infer_genetree():  # default --output is gene-trees
    run_cmd = run('vespa infer-gene-trees frog-families tree2.tre', cwd=data_dir)
    print(run_cmd.stdout, run_cmd.stderr)


def test_codeml_setup():  # default --output is codeml
    run_cmd = run('vespa codeml-setup frog-families gene-trees --branches branches.yaml', cwd=data_dir)
    print(run_cmd.stdout, run_cmd.stderr)


def test_codeml_setup_unlabelled():
    run_cmd = run('vespa codeml-setup frog-families gene-trees --output codeml-unlabelled', cwd=data_dir)
    print(run_cmd.stdout, run_cmd.stderr)
