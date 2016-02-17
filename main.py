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

"""Main program
"""

import dbmol 
import glob
from rdkit import Chem
from rdkit import rdBase
from optparse import OptionParser
import sys, os

        
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
    
    rdk_ver  = rdBase.rdkitVersion

    print rdk_ver

    parser = OptionParser( usage = "Usage: %prog [options] <structure-file-dir>", version = "%prog v0.0" )
    
    parser.add_option("-t", "--time", default = 20 , help = " Set the maximum time to perform the mcs search between pair of molecules")
    
    parser.add_option("-p", "--parallel", default = 1, type='int' , help = " Set the parallel mode on. If an integer number N is specified, N processes will be executed")
    

    #All the following parameters have been disabled and they need to be re-implemented as in the original Lomap code if possible
    # parser.add_option( "--debug", default = False, action = "store_true", help = "turn on debugging mode." )
    # parser.add_option( "-m", "--mcs", metavar = "FILE",
    #                    help = "read MCS searching results directly from FILE and avoid searching again. " \
    #                           "FILE should be a Schrodinger canvasMCS output file in the CSV format." )
    # parser.add_option( "-o", "--output", metavar = "BASENAME", default = "simimap",
    #                    help = "output files' base name. The following files will be written: <basename>.dot, and "
    #                    "<basename>.pkl." )
    # parser.add_option( "-s", "--siminp", metavar = "BASENAME",
    #                    help = "simulation input files' base name. When this option is specified, a number of input files "
    #                    "for FEP simulations will be written out." )
    # parser.add_option( "-g", "--graph", metavar = "FILENAME", help = "use the graph as saved in file FILENAME." )
    # parser.add_option( "-b", "--build",default = False, action = "store_true" , help = "build score matrix before doing graph planning")
    # parser.add_option( "-t", "--siminp_type", metavar = "TYPE", default = "mae",
    #                    help = "simulation input file type [mae | gro]" )
    # parser.add_option( "-r", "--receptor", default = 0, metavar = "N", type = "int",
    #                    help = "specify the initial N structures as the common receptor. This option is needed when "
    #                    "you want to write out structure input files for relative binding free energy calculations." )
    # parser.add_option( "--save",  default = False, action = "store_true", help = "do not delete temporary files." )
    
    
    #A tuple of options and arguments passed by the user
    (opt, args) = parser.parse_args()

    if (len( args ) == 0) :
        parser.print_help()
        sys.exit( 0 )

         
    # This list is used as container to handle all the molecules read in by using RdKit.
    # All the molecules are instances of the allocated class Molecules
    molid_list = []

    for a in args :
        print( "Reading structures from '%s'..." % a )
        n = 0
        # The .mol2 file format is the only supported so far
        mol_fnames = glob.glob( a + "/*.mol2" )

            
        for fname in mol_fnames :
            # The RDkit molecule object is read in as mol2 file. The molecule is not sanitized and 
            # all the hydrogens are kept in place
            try:
                rdkit_mol = Chem.MolFromMol2File(fname, sanitize=False, removeHs=False)
            except:
                print('Error reading the file %s', os.path.basename( fname))
            
            mol = dbmol.Molecule(rdkit_mol, os.path.basename( fname ))
            print  '    %s ID %s' % ( os.path.basename( fname ), mol.getID() ) 
            molid_list.append(mol)

            

    print( "--------------------------------------------" )
    print( "Finish reading structure input files. %d structures in total" % len( molid_list ) )
    
    if (len( mol_fnames ) > 1) :
        build_database(molid_list, opt)

    
        
# main function         
if ("__main__" == __name__) :
    startup()
