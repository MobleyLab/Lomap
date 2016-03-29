#******************
# MODULE DOCSTRING
#******************

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

from rdkit import Chem
import numpy as np
import mcs
import graphgen
import sys,os
import math
import multiprocessing
import networkx as nx
import logging
import glob
import argparse


__all__ = ['DBMolecules', 'SMatrix', 'Molecule']


#*************************
# Molecule Database Class
#*************************

class DBMolecules(object):
    """
    
    This class is used as a container for all the Molecules  
    
    """

    # Initialization function
    def __init__(self, dir_name, time_mcs=20, parallel_mode=1, verbose_mode='off', 
                 output_mode=False, out_name='out', display_mode=False, 
                 max_graph=6, cutoff_graph=0.4): 

        """
        Initialization of  the Molecule Database Class
    
        Parameters
        ----------
        dir_name : str 
           the mol2 directory file name
        time_mcs : int
           the maximum time in seconds used to perform the MCS search
        parallel_mode : int
           the number of cores used to generate the similarity score matrices
        verbose_mode : bool
           verbose mode
        output_mode : bool
           a flag used to generate or not the output files
        out_name : str
           the file name prefix used to produce the output files
        display_mode : bool
           a flag used to display or not a network made by using matplotlib
        max_graph : int
           the maximum distance used to cluster the graph nodes
        cutoff_graph : float
           the Minimum Similarity Score (MSS) used to build the graph

        """

        #########################TO DO =>CHECK INPUTS###########################


        # Set the Logging 
        if verbose_mode == 'off':
            logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.CRITICAL)
   
        if verbose_mode == 'info':
            logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.INFO)
        
        if verbose_mode == 'pedantic':
            logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.DEBUG)
                    

        # Options to buid the MCS and other parameters
        self.options = argparse.Namespace(directory=dir_name, time=time_mcs, 
                                          parallel=parallel_mode, 
                                          verbose=verbose_mode, output=output_mode,
                                          name=out_name, display=display_mode, 
                                          max=max_graph, cutoff=cutoff_graph)
        

        # Internal list container used to store the loaded molecule objects
        self.__list = self.read_mol2_files()
    

        # Dictionary which holds the mapping between the generated molecule IDs and molecule file names
        self.dic_mapping = {}
        
        for mol in self.__list:
            self.dic_mapping[mol.getID()]=mol.getName()

        # Index used to perform index selection by using __iter__ function
        self.__ci = 0

        # Symmetric matrices used to store the mcs scoring. The matrices are subclasses of numpy 
        self.strict_mtx = SMatrix(shape=(0,))
        self.loose_mtx = SMatrix(shape=(0,))

        
        # Empty pointer to the compound graph 
        self.Graph = nx.Graph() 


    
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


    def __setitem__(self, index, molecule):
        """
        Index setting function
        
        Parameters
        ----------
        index : int 
           the molecule index
        molecule : Molecule obj
           the molecule to assign to the molecule database by selecting the index: 
           DB[index] = molecule
        
        """
        
        if not isinstance(molecule, Molecule):
            raise ValueError('The passed molecule is not a Molecule object')
        
        self.__list[index] = molecule

        
    def __add__(self, molecule):        
        """
        Add a new molecule to the molecule database    
        
        Parameters
        ----------
        molecule : Molecule obj 
           the molecule to append into the molecule database
        """
        
        if not isinstance(molecule, Molecule):
            raise ValueError('The passed molecule is not a Molecule object')
        
            self.__list.append(molecule)

    
    def nums(self):
        """
        This function recovers the total number of molecules currently stored in
        the molecule database
        """
        return len(self.__list)



    def read_mol2_files(self):
        """
        Read in all the mol2 files

        Returns
        -------
        molid_list : list of Molecule objects
           the container list of all the allocated Molecule objects

        """
        
        # This list is used as container to handle all the molecules read in by using RdKit.
        # All the molecules are instances of  Molecule class
        molid_list = []

        # List of molecule that failed to load in
        mol_error_list_fn = []
    
        logging.info(30*'-')

        # The .mol2 file format is the only supported so far
        mol_fnames = glob.glob(self.options.directory + "/*.mol2" )
    
        if (len( mol_fnames ) < 2) :
            raise ValueError('The directory %s must contain at least two mol2 files' % self.options.directory)
        
        print_cnt = 0

        for fname in mol_fnames :
        
            # The RDkit molecule object reads in as mol2 file. The molecule is not sanitized and 
            # all the hydrogens are kept in place
            rdkit_mol = Chem.MolFromMol2File(fname, sanitize=False, removeHs=False)
        
            # Reading problems
            if rdkit_mol == None:
                logging.warning('Error reading the file: %s' % os.path.basename(fname))
                mol_error_list_fn.append(os.path.basename(fname))
                continue
            
            # The Rdkit molecule is stored in a Molecule object
            mol = Molecule(rdkit_mol, os.path.basename(fname))
        

            # Cosmetic printing and status
            logging.info('ID %s\t%s' % (mol.getID(), os.path.basename(fname)))
        
            
            if print_cnt == 15:
                logging.info('ID %s\t%s' % (mol.getID(), os.path.basename(fname)))
                logging.info(3*'\t.\t.\n')
            
            if print_cnt < 15 or print_cnt == (len(mol_fnames) - 1):
                logging.info('ID %s\t%s' % (mol.getID(), os.path.basename(fname)))
                
                
            print_cnt+= 1
        
            molid_list.append(mol)

        logging.info(30*'-')

        logging.info('Finish reading input files. %d structures in total....skipped %d\n' % (len(molid_list), len(mol_error_list_fn)))
    
        if mol_error_list_fn:
            logging.warning('Skipped molecules:')
            loggign.warning(30*'-')
            for fn in  mol_error_list_fn:
                logging.warning('%s'% fn)    
            print(30*'-')
        

        return molid_list

        

        
    def compute_mtx(self, a, b, strict_mtx, loose_mtx):
        """
        Compute a chunk of the similariry score matrices. The chunk is selected 
        by the start index a and the final index b. The matrices are indeed 
        treated as linear array

        Parameters
        ----------
        a : int 
           the start index of the chunk 
        b : int
           the final index of the chunk
        
        strict_mtx: python multiprocessing array
           srict simimarity score matrix. This array is used as shared memory 
           array managed by the different allocated processes. Each process 
           operates on a separate chunk selected by the indexes a and b


        loose_mtx: python multiprocessing array
           loose similarity score matrix. This array is used as shared memory 
           array managed by the different allocated processes. Each process 
           operates on a separate chunk selected by the indexes a and b

        """
        # name = multiprocessing.current_process().name
        # print name
        # print 'a = %d, b = %d' % (a,b)
        # print '\n'  

        
        def ecr(mol_i, mol_j):
            """
            This function computes the similariry score between the passed molecules 
            by using the EleCtrostatic Rule (ECR)

            Parameters
            ----------
            mol_i : Rdkit molecule object 
               the first molecules used to calculate the ECR rule  
            mol_j : Rdkit molecule object 
               the second molecules used to calculate the ECR rule 

            Returns
            -------
            scr_ecr: float
                the calculated similarity score (1 if mol_i and mol_j have the
                same total charges, 0  otherwire)

            """
            
            total_charge_mol_i = 0.0
            
            for atom in mol_i.GetAtoms():
                total_charge_mol_i += float(atom.GetProp('_TriposPartialCharge'))

            total_charge_mol_j = 0.0
            
            for atom in mol_j.GetAtoms():
                total_charge_mol_j += float(atom.GetProp('_TriposPartialCharge'))

            if abs(total_charge_mol_j - total_charge_mol_i) < 1e-3:
                scr_ecr = 1.0
            else:
                scr_ecr = 0.0
            
            return scr_ecr
        
        # Total number of loaded molecules
        n = self.nums()
        
        # Looping over all the elements of the selected matrix chunk
        for k in range(a, b+1):

            # The linear index k is converted into the row and column indexes of
            # an hypothetical bidimensional symmetric matrix
            i = int(n - 2 - math.floor(math.sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5))
            j = int(k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2)
            #print 'k = %d , i = %d , j = %d' % (k,i,j)
    
            # The Rdkit molecules moli and molj are extracted form the molecule database
            moli = self[i].getMolecule()
            molj = self[j].getMolecule()

            #print 'Processing molecules:\n%s\n%s' % (self[i].getName(),self[j].getName())

            # The Electrostatic score rule is calculated
            ecr_score = ecr(moli, molj)

            # The MCS is computed just if the passed molecules have the same charges 
            if ecr_score == 1.0:
                try: 
                    if self.options.verbose == 'pedantic':
                        logging.info(50*'-')
                        logging.info('MCS molecules: %s - %s' % (self[i].getName(), self[j].getName())) 
                    
                    # Maximum Common Subgraph (MCS) calculation    
                    MC = mcs.MCS(moli, molj, options=self.options)

                except Exception as e:
                    if self.options.verbose == 'pedantic':
                        logging.warning('Skipping MCS molecules: %s - %s\t\n\n%s' % (self[i].getName(), self[j].getName(), e))
                        logging.info(50*'-')
                    continue
            else:
                continue


            # The scoring between the two molecules is performed by using different rules.
            # The total score will be the product of all the single rules
               
            tmp_scr = ecr_score * MC.mncar() * MC.mcsr()
            
            strict_scr = tmp_scr *  MC.tmcsr(strict_flag=True) 
            loose_scr = tmp_scr * MC.tmcsr(strict_flag=False) 
        
            strict_mtx[k] = strict_scr
            loose_mtx[k] = loose_scr
    
        return


    def build_matrices(self):
        """
        This function coordinates the calculation of the similarity score matrices
        by distribuiting chunks of the matrices between the allocated processes

        """
        
        logging.info('\nMatrix scoring in progress....\n')   
        
        # The similarity score matrices are defined instances of the class SMatrix
        # which implemets a basic class for symmetric matrices
        self.strict_mtx = SMatrix(shape=(self.nums(),))
        self.loose_mtx = SMatrix(shape=(self.nums(),))

        # The total number of the effective elements present in the symmetric matrix
        l = self.nums()*(self.nums() - 1)/2

        
        if self.options.parallel == 1: # Serial execution
            self.compute_mtx(0, l-1, self.strict_mtx, self.loose_mtx)
        else: # Parallel execution
            
            logging.info('Parallel mode is on')
            
            # Number of selected processes
            np = self.options.parallel

            delta = l/np
            rem = l%np

            if delta < 1:
                kmax = l
            else:
                kmax = np

            proc = []

            # Shared memory array used by the different allocated processes
            strict_mtx = multiprocessing.Array('d', self.strict_mtx)
            loose_mtx =  multiprocessing.Array('d', self.loose_mtx)

            # Chopping the indexes ridistribuiting the remainder
            for k in range(0, kmax):
    
                spc = delta + int(rem/(k+1) > 0)
    
                if k == 0:
                    i = 0
                else:
                    i = j + 1

                if k!= kmax - 1:
                    j = i + spc - 1
                else:
                    j = l - 1

                # Python multiprocessing allocation
                p = multiprocessing.Process(target=self.compute_mtx , args=(i, j, strict_mtx, loose_mtx))
                p.start()
                proc.append(p)
            # End parallel execution        
            for p in proc:
                p.join()
          
            # Copying back the results
            self.strict_mtx[:] = strict_mtx[:]
            self.loose_mtx[:] = loose_mtx[:]
            
        return (self.strict_mtx, self.loose_mtx)

    def build_graph(self):
        """
        This function coordinates the Graph generation

        """
        logging.info('\nGenerating graph in progress....')

        # The Graph is build from an instance of the Class GraphGen by passing
        # the selected user options
        Gr = graphgen.GraphGen(self)

        # Writing the results is files
        if self.options.output:
            try:
                Gr.writeGraph()
            except Exception as e:
                logging.error(str(e))

        # Handle to the the NetworkX generated graph
        self.Graph = Gr.getGraph()

        #print self.Graph.nodes(data=True)
       
        # Display the graph by using Matplotlib
        if self.options.display:
            Gr.draw()

        return self.Graph


    def write_dic(self):
        """
        This function write out a text file with the mapping between the 
        generated molecule indexes and the corresponding molecule file names

        """

        try:
            file_txt = open(self.options.name+'.txt', 'w')
        except Exception:
            raise IOError('It was not possible to write out the mapping file:')

        file_txt.write('#ID\tFileName\n')
        for key in self.dic_mapping:
            file_txt.write('%d\t%s\n' % (key, self.dic_mapping[key]))

        file_txt.close() 

#*************************
# Symmetric  Class
#*************************

class SMatrix(np.ndarray):
    """
    This class implements a "basic" interface for symmetric matrices 
    subclassing ndarray. The class interanlly stores a bidimensional 
    numpy array as a linear array A[k], however the user can still 
    access to the matrix elements by using a two indeces notation A[i,j]
 
    """

    def __new__(subtype, shape, dtype=float, buffer=None, offset=0, strides=None, order=None):
        
        if len(shape) > 2:
            raise ValueError('...0...')

        elif len(shape) == 2:
            if shape[0] != shape[1]:
                raise ValueError('...1...')

        n = shape[0]        
        l = shape[0]*(shape[0] - 1)/2

        shape = (l,)
        
        obj = np.ndarray.__new__(subtype, shape , dtype, buffer, offset, strides, order)

        # Array inizialization
        obj = obj*0.0
        
        return obj
        

    def __getitem__(self, *kargs):
        """
        This function retrieves the selected elements i,j from the symmetric
        matrix A[i,j]
        
        Parameters
        ----------
        *kargs : python tuples
           the passed elements i,j  
        
        Returns
        -------
            : float
            the selected element extracted from the allocated linear array
            
        """
        
        if isinstance( kargs[0], int ):
            k = kargs[0]
            return super(SMatrix, self).__getitem__(k)

        elif len(kargs[0]) > 2:
            raise ValueError('Tuple dimension must be two')
                
        i = kargs[0][0]
        j = kargs[0][1]

        if i == j:
            return 0.0
        
        # Length of the linear array 
        l = self.size
        
        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        n = int((1+math.sqrt(1+8*l))/2)

        if i > n - 1:
            raise ValueError('First index out of bound')
        
        if j > n - 1:
            raise ValueError('Second index out of bound')
      
        if i < j:
            k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
        else:
            k = (n*(n-1)/2) - (n-j)*((n-j)-1)/2 + i - j - 1 
        
        return super(SMatrix, self).__getitem__(k)

   
    def __setitem__(self, *kargs):
        """
        This function set the matrix elements i,j to the passed value
        
        Parameters
        ----------
        *kargs : python tuples
           the passed elements i,j, value to set  
        
            
        """
             
        if isinstance( kargs[0], int ):
            k = kargs[0]
            value = kargs[1]
            return super(SMatrix, self).__setitem__(k,value)

        elif len(kargs[0]) > 2:
            raise ValueError('Tuple dimension must be two')
      
        # Passed indexes and value to set
        i = kargs[0][0]
        j = kargs[0][1]
        value = kargs[1]
        
        # Length of the linear array 
        l = self.size
        
        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        n = int((1+math.sqrt(1+8*l))/2)
        
        if i > n - 1:
            raise ValueError('First index out of bound')
        if j > n - 1:
            raise ValueError('Second index out of bound')
        
        if i < j:
            k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
        else:
            k = (n*(n-1)/2) - (n-j)*((n-j)-1)/2 + i - j - 1 
        
        super(SMatrix, self).__setitem__(k,value)


    def to_numpy_2D_array(self) :
        """
        This function returns the symmetric similarity score numpy matrix 
        generated from the linear array
        
        Returns
        -------
        np_mat : numpy matrix
           the symmetric similarity score numpy matrix built by using the linear
           array 
        
        """

        # Length of the linear array 
        l = self.size
        
        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        n = int((1+math.sqrt(1+8*l))/2)

        np_mat = np.zeros((n,n))

        for i in range (0,n):
            for j in range(0,n):
                np_mat[i,j] = self[i,j]

        return np_mat
        
    def mat_size(self) :
        """
        This function returns the size of the square similarity score matrix 
        
        Returns
        -------
        n : int
           the size of the similarity score matrix
        
        """ 

        # Length of the linear array 
        l =  self.size
        
        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        n = int((1+math.sqrt(1+8*l))/2)

        return n 
        



#*************************
# Molecule Class
#*************************

class Molecule(object):
    """
    This Class stores the Rdkit molecule objects, their identification number 
    and the total number of instantiated molecules 

    """

    # This variable is used to count the current total number of molecules
    # The variable is defined as private
    __total_molecules = 0

    
    def __init__(self, molecule, molname):
        """
        Initialization class function 
        
        Parameters
        ----------
        molecule : Rdkit molecule object
           the molecule
        molname : str
           the molecule file name
        
        """

        #Check Inputs
        if not isinstance(molecule, Chem.rdchem.Mol):
            raise ValueError('The passed molecule object is not a RdKit molecule')


        if not isinstance(molname, str):
             raise ValueError('The passed molecule name must be a string')
        

        # The variable __molecule saves the current RDkit molecule object
        # The variable is defined as private
        self.__molecule = molecule    
        
            
        # The variable __ID saves the molecule identification number 
        # The variable is defined as private
        self.__ID = Molecule.__total_molecules

        
        # The variable __name saves the molecule identification name 
        # The variable is defined as private
        self.__name = molname
    
        Molecule.__total_molecules+=1

    
    def getID(self):
        """
        Get the molecule ID number
 
 
        Returns
        -------
           : int
           the molecule ID number

        """
        return self.__ID

    
    def getMolecule(self):
        """
        Get the Rdkit molecule object

        Returns
        -------
        mol_copy : Rdkit molecule object
           The copy of the RDkit molecule

        """
        mol_copy = Chem.Mol(self.__molecule)
        return mol_copy


    
    def getName(self):
        """
        Get the molecule file name

        Returns
        -------
           : str 
           the molecule string file name

        """
        
        return self.__name


    
    @staticmethod
    def get_mol_num():
        """
        This class function returns the current total number of allocated molecules. 
 

        Returns
        -------
           : int 
           the total number of allocated molecules 

        """
        return Molecule.__total_molecules



# Classes used to check some of the passed user options in the main function

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
    
    # Command line user interface
    parser = argparse.ArgumentParser(description='Lead Optimization Mapper 2. A program to plan alchemical relative binding affinity calculations', prog='LOMAPv1.0')
    parser.add_argument('directory', action=check_dir,\
                        help='The mol2 file directory')
    parser.add_argument('-t', '--time', default=20, action=check_int,type=int,\
                        help='Set the maximum time in seconds to perform the mcs search between pair of molecules')
    parser.add_argument('-p', '--parallel', default=1, action=check_int,type=int,\
                        help='Set the parallel mode. If an integer number N is specified, N processes will be executed to build the similarity matrices')
    parser.add_argument('-v', '--verbose', default='off', type=str,\
                        choices=['off', 'info', 'pedantic'], help='verbose mode selection')
    
    parser.add_argument('-o', '--output', default=False, action='store_true',\
                        help='Generates output files')
    parser.add_argument('-n', '--name', default='out',\
                        help='File name prefix used to generate the output files')

    parser.add_argument('-d', '--display', default=False, action='store_true',\
                        help='Display the generated graph by using Matplotlib')
    parser.add_argument('-m', '--max', default=6, action=check_int ,type=int,\
                        help='The maximum distance used to cluster the graph nodes')
    parser.add_argument('-c', '--cutoff', default=0.4 ,type=float,\
                        help='The Minimum Similariry Score (MSS) used to build the graph')
    
    # Options and arguments passed by the user
    ops= parser.parse_args()
    
    
    # Molecule DataBase initialized with the passed user options
    db_mol = DBMolecules(ops.directory, ops.time, ops.parallel, ops.verbose, 
                        ops.output, ops.name, ops.display, ops.max, ops.cutoff) 

   
    # Similarity score linear array generation
    strict, loose =  db_mol.build_matrices()
    
    # Get the 2D numpy matrices
    # strict.to_numpy_2D_array()
    # loose.to_numpy_2D_array()
   
    # Graph generation based on the similarity score matrix
    nx_graph = db_mol.build_graph()   

    # print nx_graph.nodes(data=True)
    # print nx_graph.edges(data=True)


