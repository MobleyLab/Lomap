"""Main program
"""

from  db_mol import *
import glob
import sys, os
import logging
from rdkit import Chem
from optparse import OptionParser
from rdkit import rdBase

# Logger setup
logger = logging.getLogger()
if (logger.handlers) :
    for handler in logger.handlers :
        logger.removeHandler( handler )
        
logging.basicConfig( format  = '%(asctime)s: %(message)s',
                     datefmt = '%m/%d/%y %I:%M:%S',
                     level   = logging.INFO )



def build_database(molid_list, options) :
    """
    This function creates the molecules database and generates the matrix score.  
    
    molid_list: A list of molecules instance of the defined Molecule class
    opt: the parameters selected to perform the reseacrh of the MCS between a pair of molecules 

    """

    molecules = DB_Molecules(molid_list)

    strict_mtx, loose_mtx =  molecules.build_matrix(options)
    
    print strict_mtx
    print loose_mtx


def startup() :
    """
    The startup function, which will handle the command line interface.
    """
    
    rdk_ver  = rdBase.rdkitVersion

    print rdk_ver

    parser = OptionParser( usage = "Usage: %prog [options] <structure-file-dir>", version = "%prog v0.0" )
    parser.add_option("-t", "--time", default = 20 , help = " Set the maximum time to perform the mcs search between pair of molecules")
    
    parser.add_option( "--debug", default = False, action = "store_true", help = "turn on debugging mode." )

    #All the following parameters have been disabled and they need to be re-implemented as in the original Lomap code if possible

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

    if (opt.debug) :
        logger.setLevel( logging.DEBUG )
        logging.debug( "Debugging mode is on." )
     
    # This list is used as container to handle all the molecules read in by using RdKit.
    # All the molecules are instances of the allocated class Molecules
    molid_list = []

    for a in args :
        logging.info( "Reading structures from '%s'..." % a )
        n = 0
        # The .mol2 file format is the only supported so far
        mol_fnames = glob.glob( a + "/*.mol2" )

            
        for fname in mol_fnames :
            if (n < 8) :
                logging.info( "    %s" % os.path.basename( fname ) )
            elif (n == 8) :
                logging.info( "    (more)..." )
                break
                n += 1
                logging.info( "    %d files found." % len( mol_fnames ) )
        
            # The RDkit molecule object is read in as mol2 file. The molecule is not sanitized and 
            # all the hydrogens are kept in place
            mol = Molecule(Chem.MolFromMol2File(fname, sanitize=False, removeHs=False))

            molid_list.append(mol)


    logging.info( "--------------------------------------------" )
    logging.info( "Finish reading structure input files. %d structures in total" % len( molid_list ) )
    
    if (len( mol_fnames ) > 1) :
        build_database(molid_list, opt)

    
        
# main function         
if ("__main__" == __name__) :
    startup()
