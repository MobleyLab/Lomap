"""
Setup script for Lomap2

You can install lomap with

python setup.py install
"""

import sys,os
from os.path import relpath, join

from setuptools import setup, find_packages

if sys.argv[-1] == 'setup.py':
    print("To install, run 'python setup.py install'")
    print()

if sys.version_info[:2] < (2, 7):
    print("Lomap requires Python 2.7 or later (%d.%d detected)." %
          sys.version_info[:2])
    sys.exit(-1)


descr = """
The Lead Optimization Mapper (LOMAP) is an automated algorithm
to plan efficient relative free energy calculations between  
potential ligands within a substantial of compounds'
"""

setup(
    name                 = 'lomap', 
    version              = '0.0.0', 
    description          = 'Lead Optimization Mapper 2',
    long_description     = descr,
    url                  = 'https://github.com/nividic/Lomap',
    author               = 'Gaetano Calabro and David Mobley',
    author_email         = 'gcalabro -at- uci.edu',
    license              = 'LGPL',
    platforms            = ['Linux-64', 'Mac OSX-64', 'Unix-64'],
    packages             = find_packages()+['test'],
    include_package_data = True,
      
    entry_points         = {'console_scripts':['lomap=lomap.dbmol:startup']},
    zip_safe             = False
)

