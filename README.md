[![Anaconda Badge](https://anaconda.org/nividic/lomap/badges/version.svg)](https://anaconda.org/nividic/lomap)

[![Build Status](https://travis-ci.org/MobleyLab/Lomap.svg?branch=master)](https://travis-ci.org/MobleyLab/Lomap)

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
* RDKit Release >2015.09.2
* Graphviz 2.38
* pygraphviz
* NetworkX 
* Matplotlib 
* PyQt 4.11

Authors
-------
* Gaetano Calabro' <gcalabro@uci.edu>
* David Mobley <dmobley@uci.edu>

## Installation

Add to the conda channels:

conda config --add channels nividic

and then:

conda install lomap

or

Add to the conda channels:

conda config --add channels mobleylab

and then:

conda install lomap


Usage
-----
```python
import lomap

# Generate the molecule database starting from a directory containing .mol2 files

db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True)

# Calculate the similarity matrix betweeen the database molecules. Two molecules are generated
# related to the scrict rule and loose rule 

strict, loose = db_mol.build_matrices()

# Generate the NetworkX graph and output the results
nx_graph = db_mol.build_graph() 


# Calculate the Maximum Common Subgraph (MCS) between 
# the first two molecules in the molecule database 
# ignoring hydrogens and depicting the mapping in a file
    
MC = lomap.MCS.getMapping(db_mol[0].getMolecule(), db_mol[1].getMolecule(), hydrogens=False, fname='mcs.png')


# Alchemical transformation are usually performed between molecules with
# the same charges. However, it is possible to allow this transformation
# manually setting the electrostatic score for the whole set of molecules 
# producing a connected graph. The electrostatic scrore must be in the 
# range [0,1]


db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, ecrscore=0.1)
strict, loose = db_mol.build_matrices()
nx_graph = db_mol.build_graph() 
```



## Issues
* Lomap is in debugging stage and it has been tested on Ubuntu 14.04 and OSX Yosemite
* Lomap has been developed in python 2.7 and 3.4

## Disclaimers
* This code is currently in alpha release status. Use at your own risk. We will almost certainly be making changes to the API in the near future.
