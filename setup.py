import re
import sys

from setuptools import setup


__version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        open('vespasian/__init__.py').read()).group(1)


if sys.version_info < (3,6):
      sys.exit('Requires Python >= 3.6')


CLASSIFIERS = ['Environment :: Console',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
               'Natural Language :: English',
               'Operating System :: POSIX :: Linux',
               'Operating System :: MacOS :: MacOS X',
               'Programming Language :: Python :: 3.6',
               'Programming Language :: Python :: 3.7',
               'Topic :: Scientific/Engineering :: Bio-Informatics']


setup(name = 'vespasian',
      version = __version__,
      description = 'Genome scale evolutionary hypothesis testing',
      url = '',
      author = "Bede Constantinides, Mary O'Connell",
      author_email = 'bedeabc@gmail.com',
      license = 'LICENSE',
      packages=['vespasian'],
      zip_safe=True,
      install_requires=['biopython',
                        'treeswift',
                        'tqdm',
                        'argh',
                        'pyyaml',
                        'parmap',
                        'pytest',
                        'snakemake',
                        'pandas',
                        'numpy',
                        'scipy'],
      entry_points = {'console_scripts':['vespasian=vespasian.cli:main']})
