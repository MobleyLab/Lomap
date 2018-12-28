"""
Setup script for Lomap

You can install lomap with

python setup.py install
"""

import sys

from setuptools import setup, find_packages

if sys.argv[-1] == 'setup.py':
    print("To install, run 'python setup.py install'")
    print()

if sys.version_info[:2] < (3, 6):
    print("Lomap requires Python 3.6 or later (%d.%d detected)." %
          sys.version_info[:2])
    sys.exit(-1)


description = """
The Lead Optimization Mapper (LOMAP) is an automated algorithm
to plan efficient relative free energy calculations between
potential ligands within a substantial of compounds'
"""

setup(
    name='lomap',
    version='0.0.0',
    description='Lead Optimization Mapper',
    long_description=description,
    url='https://github.com/MobleyLab/Lomap',
    author='Gaetano Calabro and David Mobley',
    author_email='gcalabro -at- eyesopen.com',
    license='MIT',
    platforms=['Linux-64', 'Mac OSX-64'],
    packages=find_packages()+['tests'],
    include_package_data=True,
    zip_safe=False
)

