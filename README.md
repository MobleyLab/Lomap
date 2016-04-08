[![Anaconda Badge](https://anaconda.org/nividic/lomap/badges/version.svg)](https://anaconda.org/nividic/lomap)

# Lomap
Alchemical free energy calculations hold increasing promise 
as an aid to drug discovery efforts. However, applications of 
these techniques in discovery projects have been relatively 
few, partly because of the difficulty of planning and setting up 
calculations. The lead optimization mapper (LOMAP) was 
introduced as an automated algorithm to plan efficient relative 
free energy calculations between potential ligands within 
a substantial of compounds. The original LOMAP code was mainly
based on commercial APIs such as OpenEye and Schrodinger. The aim 
of this project is to deveop a new version of LOMAP based on free
avalaible APIs such as RDKit offering the scientific community a 
free tool to plan in advance binding free energy calculations


## Prerequisites
* RDKit Release 2015.09.2
* Graphviz 2.38
* pygraphviz 1.12
* NetworkX 1.11
* Matplotlib 1.5.1

Authors
-------
* Gaetano Calabro' <gcalabro@uci.edu>
* David Mobley <dmobley@uci.edu>

## Installation

Add to the conda channels:

conda config --add channels nividic

and then:

conda install lomap

Usage
-----

import lomap

db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files")

strict, loose = db_mol.build_matrices()

nx_graph = db_mol.build_graph() 



## Issues
* Lomap is in debugging stage and it has been tested on Ubuntu 14.04 and OSX Yosemite


## Disclaimers
* This code is currently in alpha release status. Use at your own risk. We will almost certainly be making changes to the API in the near future.