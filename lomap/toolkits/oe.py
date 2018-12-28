# *****************************************************************************
# Lomap: A toolkit to plan alchemical relative binding affinity calculations
# Copyright 2015 - 2018  UC Irvine and the Authors
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


from openeye import oechem
import os
from glob import iglob
import logging


def read_molecules(directory):
    """
    Read in all the molecules files

    Parameters
    ----------
    directory: str
        the directory containing the molecule files to read in

    Returns
    -------
    oemol_list : list of Molecule objects
        the container list of all the allocated Molecule objects
    """

    from lomap.utils import Molecule

    fns = [f for f in iglob(directory+'/**', recursive=True) if os.path.isfile(f)]
    fns.sort()

    # This list is used as container to handle all the molecules read in by using OpenEye.
    # All the molecules are instances of  Molecule class
    oemol_list = []

    # List of molecule that failed to load in
    mol_error_list_fn = []

    for fn in fns:
        with oechem.oemolistream(fn) as ifs:
            ifs_test = oechem.oemolistream(fn)
            mol_test = oechem.OEMol()
            if oechem.OEReadMolecule(ifs_test, mol_test):
                for mol in ifs.GetOEMols():
                    title = mol.GetTitle()
                    if not title:
                        title = os.path.basename(fn)
                    mol_toolkit = Molecule(mol, title)
                    oemol_list.append(mol_toolkit)
                ifs.close()
            else:
                logging.warning('Error reading the file: {}'.format(fn))
                mol_error_list_fn.append(fn)
                continue

            ifs_test.close()

    logging.info('Finish reading input files. {} structures in total....skipped {}'.format(len(oemol_list),
                                                                                           len(mol_error_list_fn)))
    if mol_error_list_fn:
        logging.warning('\nSkipped molecules:')
        logging.warning(30*'-')
        for fn in mol_error_list_fn:
            logging.warning(fn)
        print(30*'-')

    return oemol_list
