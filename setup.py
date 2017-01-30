"""
A setuptools-based setup module.
Mostly rewritten from https://github.com/pypa/sampleproject/blob/master/setup.py
and https://packaging.python.org/distributing
"""

from setuptools import setup

setup(name='kmers',
      version='0.1.0',
      description='Bioinformatic k-mer analysis',
      long_description='Calculating k-mer distributions of (sets of) sequences, classifying sequences and calculating distances',
      url='https://github.com/synedraacus/kmers',
      author='A.A. Morozov', author_email='morozov@lin.irk.ru',
      license=''
      )