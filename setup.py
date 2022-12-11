#!/usr/bin/env python

from setuptools import setup, find_packages

with open("khipu/__init__.py") as f:
    exec([x for x in f.readlines() if '__version__' in x][0])

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
  name='khipu-metabolomics',
  version=__version__,

  author='Shuzhao Li',
  author_email='shuzhao.li@gmail.com',
  description='Common utilities for interpreting mass spectrometry data',
  long_description_content_type="text/markdown",
  long_description=long_description,
  url='https://github.com/shuzhao-li/khipu',
  license='BSD',

  keywords='chemistry, bioinformatics, mass spectrometry',

  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3.7',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  # changed from earlier setuptools
  packages=find_packages(
    include=['*', '']
  ),
  include_package_data=True,
  install_requires=[
    'mass2chem',
    'treelib',
    'networkx'
  ],

  python_requires='>=3.7',



)
