import os
import sys

import subprocess

from Bio import SeqIO

subprocess.run('vespa infer-gene-trees'
			   ' /Users/bede/Research/Tools/vespa-slim/tests/data/frog_families'
			   ' /Users/bede/Research/Tools/vespa-slim/tests/data/tree2.tre'
			   ' ~/Desktop/out')

# vespa infer-gene-trees /Users/bede/Research/Tools/vespa-slim/tests/data/frog_families /Users/bede/Research/Tools/vespa-slim/tests/data/tree2.tre ~/Desktop/out

# vespa codeml-setup /Users/bede/Research/Tools/vespa-slim/tests/data/frog_families /Users/bede/Desktop/out /Users/bede/Research/Notebooks/res/2018-07-18/branches.yaml /Users/bede/Desktop/outty