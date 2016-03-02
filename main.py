#############################################################################
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
#############################################################################

"""
Main program
"""

import dbmol 
import glob
from rdkit import Chem
from rdkit import rdBase
import argparse
import sys, os
import logging
        
def build_database(molid_list, options) :
    """
    This function creates the molecules database and generates the matrix score.
    
    molid_list: A list of molecules instance of the defined Molecule class
    opt: the parameters selected to perform the reseacrh of the MCS between 
    a pair of molecules 

    """

    db_mol = dbmol.DBMolecules(molid_list,options)

    db_mol.build_matrices()
    
    #print db_mol.strict_mtx
    #print db_mol.loose_mtx
    db_mol.build_graph()
    

def startup() :
    """
    The startup function, which will handle the command line interface.
    """
                
    class check_dir(argparse.Action):
        def __call__(self, parser, namespace, directory, option_string=None):

            if not os.path.isdir(directory):
                raise argparse.ArgumentTypeError('The directory name is not a valid path: %s' % directory)

            if os.access(directory, os.R_OK):
                setattr(namespace,self.dest, directory)
            else:
                raise argparse.ArgumentTypeError('The directory name is not readable: %s' % directory)
    

    class check_int(argparse.Action):
        def __call__(self, parser, namespace, value, option_string=None):
            
            if value < 1:
                raise argparse.ArgumentTypeError('%s is not a positive integer number' % value)
                
            setattr(namespace, self.dest, value)

    

    parser = argparse.ArgumentParser(description='Lead Optimization Mapper 2. A program to plan alchemical relative binding affinity calculations', prog='LOMAPv1.0')

    parser.add_argument('directory', action=check_dir, help='The mol2 file directory')
    parser.add_argument('-t', '--time', default=20, action=check_int, type=int, help='Set the maximum time in seconds to perform the mcs search between pair of molecules.')
    parser.add_argument('-p', '--parallel', default=1, action=check_int, type=int, help='Set the parallel mode. If an integer number N is specified, N processes will be executed to build the similarity matrices.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Turn on the verbose mode.' )
    parser.add_argument('-o', '--output', default='out', help='Output file base name')
    parser.add_argument('-g', '--graph', default=False, action='store_true', help='Display the generated graph by using Matplotlib')
    parser.add_argument('-m', '--max', default=6, action=check_int ,type=int, help='The maximum distance used to cluster the graph nodes')
    parser.add_argument('-c', '--cutoff', default=0.4 ,type=float, help='The Minimum Similariry Score (MSS) used to build the graph')
    
    #Options and arguments passed by the user
    args=parser.parse_args()
    

    logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.INFO)
    
    
    # This list is used as container to handle all the molecules read in by using RdKit.
    # All the molecules are instances of the allocated class Molecules
    molid_list = []
    mol_error_list_fn = []
    
    print('\nReading structure mol2 files from directory: %s' % args.directory)
    print(30*'-')

    # The .mol2 file format is the only supported so far
    mol_fnames = glob.glob(args.directory + "/*.mol2" )
    
    print_cnt = 0
    
    for fname in mol_fnames :
        # The RDkit molecule object is read in as mol2 file. The molecule is not sanitized and 
        # all the hydrogens are kept in place
        
        rdkit_mol = Chem.MolFromMol2File(fname, sanitize=False, removeHs=False)
        
        if rdkit_mol == None:
            logging.warning('Error reading the file: %s' % os.path.basename(fname))
            mol_error_list_fn.append(os.path.basename(fname))
            continue
            
        mol = dbmol.Molecule(rdkit_mol, os.path.basename(fname))
        
        if args.verbose:
            logging.info('ID %s\t%s' % (mol.getID(), os.path.basename(fname)))
        
        else:
            if print_cnt == 15:
                logging.info('ID %s\t%s' % (mol.getID(), os.path.basename(fname)))
                print(3*'\t.\t.\n')
            
            if print_cnt < 15 or print_cnt == (len(mol_fnames) - 1):
                logging.info('ID %s\t%s' % (mol.getID(), os.path.basename(fname)))
                
                
        print_cnt+= 1
        
        molid_list.append(mol)

    print(30*'-')

    print('Finish reading structure input files. %d structures in total....skipped %d\n' % (len(molid_list), len(mol_error_list_fn)))
    
    if mol_error_list_fn:
        print('Skipped molecules:')
        print(30*'-')
        for fn in  mol_error_list_fn:
            logging.warning('%s'% fn)    
        print(30*'-')

    if (len( mol_fnames ) > 1) :
        build_database(molid_list, args)

    
        
# main function         
if ("__main__" == __name__) :
    startup()
