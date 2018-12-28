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


"""
LOMAP
=====

Alchemical free energy calculations hold increasing promise as an aid to drug
discovery efforts. However, applications of these techniques in discovery
projects have been relatively few, partly because of the difficulty of planning
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an
automated algorithm to plan efficient relative free energy calculations between
potential ligands within a substantial of compounds.
"""

import logging

from lomap import toolkits


class DBMolecules(object):
    """

    This class is used as a container for all the Molecules

    """

    # Initialization function
    def __init__(self, directory, verbose='info'):

        """
        Initialization of  the Molecule Database Class

        Parameters
        ----------
        directory : str
           the mol2 directory file name
        verbose : bool
           verbose mode
        """

        # Set the Logging
        if verbose == 'off':
            logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

        if verbose == 'info':
            logging.basicConfig(format='%(message)s', level=logging.INFO)

        if verbose == 'pedantic':
            logging.basicConfig(format='%(message)s', level=logging.DEBUG)

        # Internal list container used to store the loaded molecule objects
        self.__list = toolkits.read_molecules(directory)

        # Index used to perform index selection by using __iter__ function
        self.__ci = 0

        self.__toolkit_type = toolkits.DEFAULT

    @property
    def nums(self):
        """
        This function recovers the total number of molecules currently stored in
        the molecule database
        """
        return len(self.__list)

    def __iter__(self):
        """
        Index generator
        """
        return self

    def __next__(self):
        """
        Select the molecule during an iteration
        """

        if self.__ci > len(self.__list) - 1:
            self.__ci = 0
            raise StopIteration
        else:
            self.__ci = self.__ci + 1
            return self.__list[self.__ci - 1]

    def __getitem__(self, index):
        """
        Slicing and index selection function
        """

        return self.__list[index]
