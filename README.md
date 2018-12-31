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
* Python 3.6
* RDKit Release >= 2018.09.1.0
* (optional) Openeye-toolkits >= 2018.10.1



Authors
-------
* Gaetano Calabro' <gcalabro@uci.edu>
* David Mobley <dmobley@uci.edu>

## Installation

by using the local env file:

* conda env create -f env.yaml -n lomap

* source activate lomap

inside the lomap directory type"

* python setup.py develop



Usage
-----
```python
import lomap

# If multiple toolkits are present RDK will be the default one
db_rdk = lomap.DBMolecules("./tests/data")

# Checking RDK mol
print(db_rdk[0].molecule)

# Switching Toolkit

try:
    lomap.toolkits.set_default("OE")

    db_oe = lomap.DBMolecules("./tests/data")

    # Checking OE mol
    print(db_oe[0].molecule)
except:
    pass
```



## Issues
* Lomap is in debugging stage and it has been tested on Ubuntu 16.04 and OSX High Sierra
* Lomap has been developed in python 3.6

## Disclaimers
* This code is currently in alpha release status. Use at your own risk. We will almost certainly be making changes to the API in the near future.
