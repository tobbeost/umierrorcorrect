#!/usr/bin/env python3
from setuptools import setup
from setuptools import find_packages

pack=find_packages()

setup(name='umierrorcorrect',
      version='0.1',
      description='UMI error correct',
      long_description = open('README.md').read(),
      url='http://github.com/tobbeost/umierrorcorrect',
      author='Tobias Osterlund',
      author_email='tobias.osterlund@gu.se',
      packages=pack,
      package_data={'umierrorcorrect': ['README.md']
                   },
      include_package_data=True,
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics'],
      scripts=['run_umierrorcorrect.py', 
                'preprocess.py',
                'run_mapping.py',
                'umi_error_correct.py'
                'filter_bam.py'
                'filter_cons.py'],
      zip_safe=False)
