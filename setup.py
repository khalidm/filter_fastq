#!/usr/bin/python

from distutils.core import setup

LONG_DESCRIPTION = \
'''Filter FASTQ file'''

setup(
    name='filter_fastq',
    version='0.1.0.0',
    author='Khalid Mahmood',
    author_email='kmahmood@unimelb.edu.au',
    packages=['filter_fastq'],
    package_dir={'filter_fastq': 'filter_fastq'},
    entry_points={
        'console_scripts': ['filter_fastq = filter_fastq.filter_fastq:main']
    },
    url='https://github.com/khalidm/filter_fastq',
    license='LICENSE',
    description=('Filter out reads from FASTQ file'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["biopython"],
)
