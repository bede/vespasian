import re
import sys

from setuptools import setup


__version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        open('vespa/__init__.py').read()).group(1)


if sys.version_info < (3,6):
      sys.exit('Requires Python >= 3.6')


CLASSIFIERS = ['Environment :: Console',
               'Environment :: MacOS X',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
               'Natural Language :: English',
               'Operating System :: POSIX :: Linux',
               'Operating System :: MacOS :: MacOS X',
               'Programming Language :: Python :: 3.6',
               'Programming Language :: Python :: 3.7',
               'Topic :: Scientific/Engineering :: Bio-Informatics']


setup(name = 'vespa-slim',
      version = __version__,
      description = 'Toolkit for genome scale evolutionary hypothesis testing',
      url = '',
      author = "Bede Constantinides, Mary O'Connell",
      author_email = 'bedeabc@gmail.com',
      license = 'LICENSE',
      packages=['vespa'],
      zip_safe=True,
      install_requires=['biopython==1.73',
                        'treeswift==1.1.0',
                        'tqdm==4.31.1',
                        'argh==0.26.2',
                        'pyyaml==3.13',
                        'six==1.12.0',
                        'ete3==3.1.1',
                        'parmap==1.5.1',
                        'pytest'],
      entry_points = {'console_scripts':['vespa=vespa.cli:main']})
