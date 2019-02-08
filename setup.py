#!/usr/bin/env python3
from setuptools import setup
from setuptools import find_packages

pack = find_packages()

install_requires = ["pysam>=0.8.4"]


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
      install_requires=install_requires,
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics'],
      scripts=['umierrorcorrect/run_umierrorcorrect.py',
                'umierrorcorrect/preprocess.py',
                'umierrorcorrect/run_mapping.py',
                'umierrorcorrect/umi_error_correct.py',
                'umierrorcorrect/filter_bam.py',
                'umierrorcorrect/filter_cons.py'],
      zip_safe=False)
