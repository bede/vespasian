import re
import sys

from setuptools import setup

import vespasian


with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(name = 'vespasian',
      version = vespasian.__version__,
      description = 'Genome scale evolutionary hypothesis testing',
      # long_description=long_description,
      # long_description_content_type='text/markdown',
      # url = 'http://github.com/bede/konstel',
      author = 'Bede Constantinides',
      author_email = 'bedeabc@gmail.com',
      license = 'LICENSE',
      packages=['vespasian'],
      zip_safe=True,
      python_requires='>=3.7',
      install_requires=['biopython>=1.78',
                        'treeswift==1.1.14',
                        'tqdm',
                        'argh',
                        'pyyaml',
                        'parmap',
                        'snakemake',
                        'pandas>=1.2.4',
                        'numpy>=1.20.2',
                        'scipy>=1.6.2'],
      entry_points = {'console_scripts':['vespasian=vespasian.cli:main']},
      classifiers=['Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Natural Language :: English'])
