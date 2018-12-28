# *****************************************************************************
# Lomap: A toolkit to plan alchemical relative binding affinity calculations
# Copyright 2015 - 2019  UC Irvine and the Authors
#
# Authors: Dr Gaetano Calabro' and Dr David Mobley
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the MIT License along with this library;
# if not, see https://opensource.org/licenses/MIT
# *****************************************************************************

from lomap.utils import Molecule

import logging

import os

from glob import iglob

from rdkit import Chem


def read_molecules(directory):
    """
    Read in all the molecules files

    Parameters
    ----------
    directory: str
        the directory containing the molecule files to read in

    Returns
    -------
    rdkmol_list : list of Molecule objects
        the container list of all the allocated Molecule objects
    """

    logging.warning("RDK reading is restricted to single Tripos mol2 files")

    fns = [f for f in iglob(directory+'/**/*.mol2', recursive=True) if os.path.isfile(f)]
    fns.sort()

    # This list is used as container to handle all the molecules read in by using OpenEye.
    # All the molecules are instances of  Molecule class
    rdkmol_list = []

    # List of molecule that failed to load in
    mol_error_list_fn = []

    for fn in fns:

        # The RDkit molecule object reads in as mol2 file. The molecule is not sanitized and
        # all the hydrogens are kept in place
        mol = Chem.MolFromMol2File(fn, sanitize=False, removeHs=False)

        # Reading problems
        if mol is None:
            logging.error('Error reading the file: {}'.format(os.path.basename(fn)))
            mol_error_list_fn.append(os.path.basename(fn))
            continue

        mol_toolkit = Molecule(mol, os.path.basename(fn))
        rdkmol_list.append(mol_toolkit)

    logging.info('Finish reading input files. {} structures in total....skipped {}'.format(len(rdkmol_list),
                                                                                           len(mol_error_list_fn)))
    if mol_error_list_fn:
        logging.warning('\nSkipped molecules:')
        logging.warning(30*'-')
        for fn in mol_error_list_fn:
            logging.warning(fn)
        print(30*'-')

    return rdkmol_list
