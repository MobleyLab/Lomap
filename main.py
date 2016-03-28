#******************
# MODULE DOCSTRING
#******************

"""

LOMAP: Main Function
====================

Alchemical free energy calculations hold increasing promise as an aid to drug 
discovery efforts. However, applications of these techniques in discovery 
projects have been relatively few, partly because of the difficulty of planning 
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an 
automated algorithm to plan efficient relative free energy calculations between 
potential ligands within a substantial of compounds.

"""

#*****************************************************************************
# Lomap2: A toolkit to plan alchemical relative binding affinity calculations
# Copyright 2015 - 2016  UC Irvine and the Authors
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
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see http://www.gnu.org/licenses/
#*****************************************************************************

#****************
# MODULE IMPORTS
#****************

import dbmol 
import argparse
import logging
import os        

# Start up function
def startup() :
    """
    
    The startup function, which handles the command line user interface
    
    """
                
    # Command line user interface

    parser = argparse.ArgumentParser(description='Lead Optimization Mapper 2. A program to plan alchemical relative binding affinity calculations', prog='LOMAPv1.0')
    parser.add_argument('directory', action=check_dir, help='The mol2 file directory')
    parser.add_argument('-t', '--time', default=20, action=check_int, type=int, help='Set the maximum time in seconds to perform the mcs search between pair of molecules.')
    parser.add_argument('-p', '--parallel', default=1, action=check_int, type=int, help='Set the parallel mode. If an integer number N is specified, N processes will be executed to build the similarity matrices.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Turn on the verbose mode.' )
    
    parser.add_argument('-o', '--output', default=False, action='store_true', help='Generates output files')
    parser.add_argument('-n', '--name', default='out', help='File name prefix used to generate the output files')

    parser.add_argument('-d', '--display', default=False, action='store_true', help='Display the generated graph by using Matplotlib')
    parser.add_argument('-m', '--max', default=6, action=check_int ,type=int, help='The maximum distance used to cluster the graph nodes')
    parser.add_argument('-c', '--cutoff', default=0.4 ,type=float, help='The Minimum Similariry Score (MSS) used to build the graph')
    
    # Options and arguments passed by the user
    options=parser.parse_args()
    
    # Set the Logging 
    logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.INFO)
    
    # Molecule DataBase initialized with the passed user options
    db_mol = dbmol.DBMolecules(options)
   
    # Similarity score matrix generation
    db_mol.build_matrices()
        
    # Graph generation based on the similarity score matrix
    db_mol.build_graph()   

    
# Classes used to check some of the passed user options

# Class used to check the input directory 
class check_dir(argparse.Action):
    def __call__(self, parser, namespace, directory, option_string=None):
        if not os.path.isdir(directory):
            raise argparse.ArgumentTypeError('The directory name is not a valid path: %s' % directory)
        if os.access(directory, os.R_OK):
            setattr(namespace,self.dest, directory)
        else:
            raise argparse.ArgumentTypeError('The directory name is not readable: %s' % directory)
    
# Class used to check the parallel and time user options
class check_int(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if value < 1:
            raise argparse.ArgumentTypeError('%s is not a positive integer number' % value)
        setattr(namespace, self.dest, value)

    


# Main function         
if ("__main__" == __name__) :
    startup()
