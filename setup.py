"""
A setuptools-based setup module.
Mostly rewritten from https://github.com/pypa/sampleproject/blob/master/setup.py
and https://packaging.python.org/distributing
"""

from setuptools import setup, find_packages

setup(name='kmers',
      version='0.2.0',
      description='Bioinformatic k-mer analysis',
      long_description='Calculating k-mer distributions of (sets of) sequences, classifying sequences and calculating distances',
      url='https://github.com/synedraacus/kmers',
      author='A.A. Morozov', author_email='morozov@lin.irk.ru',
      license='MIT',
      classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'Intended Audience :: Education',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            #  Why do I even have to write it manually?
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
      keywords='bioinformatics kmers classification',
      packages=find_packages(),
      install_requires=['biopython']
      )
