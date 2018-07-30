import os
import sys
import subprocess


cwd = '/Users/bede/Research/Tools/vespa-slim/tests/data'

def run(cmd, cwd=cwd):
    return subprocess.run(cmd,
                          cwd=cwd,
                          shell=True,
                          check=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def setup_pipeline():
    run_cmd('rm -rf gene-trees codeml', cwd=cwd)  # Start afresh


def test_infer_genetree():
    run_cmd = run('vespa infer-gene-trees frog-families tree2.tre gene-trees', cwd=cwd)
    print(run_cmd.stdout, run_cmd.stderr)


def test_codeml_setup():
    run_cmd = run('vespa codeml-setup frog-families gene-trees branches.yaml codeml', cwd=cwd)
    print(run_cmd.stdout, run_cmd.stderr)
