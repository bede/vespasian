import re
import sys

from setuptools import setup


__version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        open('vespa/__init__.py').read()).group(1)


if sys.version_info[0] < 3:
      sys.exit('Requires Python >= 3(.5)')


CLASSIFIERS = ['Environment :: Console',
               'Environment :: MacOS X',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
               'Natural Language :: English',
               'Operating System :: POSIX :: Linux',
               'Operating System :: MacOS :: MacOS X',
               'Programming Language :: Python :: 3.5',
               'Programming Language :: Python :: 3.6',
               'Programming Language :: Python :: 3.7',
               'Topic :: Scientific/Engineering :: Bio-Informatics']


setup(name = 'vespa-slim',
      version = __version__,
      description = 'Selective pressure analysis toolkit',
      url = '',
      author = "Bede Constantinides, Mary O'Connell",
      author_email = 'bedeabc@gmail.com',
      license = 'LICENSE',
      packages=['vespa'],
      zip_safe=True,
      install_requires=['biopython>=1.72',
                        'treeswift>=1.0.57',
                        'tqdm>=4.24.0',
                        'argh>=0.26.2',
                        'pyyaml>=3.13',
                        'six==1.11.0',
                        'ete3==3.1.1'],
      entry_points = {'console_scripts':['vespa=vespa.cli:main']})

# MAFFT, Muscle, IQ-TREE