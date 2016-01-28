from rdkit import Chem
import numpy as np
import mcs
import sys
import math
import multiprocessing


class DBMolecules(object):
    """
    This class is used a container for all Molecule objects.  
    The class implements indexing and slicing  
    """

    def __init__(self, molecules):
        """
        This function initializes the molecule database
    
        molecules: A list of molecules instance of the defined Molecule class
       
        """
        
        #Check the passed molecules objects
        if not isinstance(molecules, list):
            raise ValueError('The passed parents must be a list')

        for mol in molecules:
             if not isinstance(mol , Molecule) :
                 raise ValueError('The passed molecule list does not contain all Molecule objects')
        
        # list container used to store the loaded molecules
        self.__list = molecules

        # index used to perform index selection by using __iter__ function
        self.__ci = 0

        # symmetric matrices used to store the mcs scoring 
        self.strict_mtx = SMatrix(shape=(0,))
        self.loose_mtx = SMatrix(shape=(0,))

    # index generator
    def __iter__(self):
        return self

    # select the molecule during an iteration
    def next(self):
        if self.__ci > len(self.__list) - 1:
            self.__ci = 0
            raise StopIteration
        else:
            self.__ci = self.__ci + 1
            return self.__list[self.__ci - 1] 
            
    # slicing and index selection function
    def __getitem__(self,index):
         return self.__list[index]

    # index setting function     
    def __setitem__(self,index,molecule):
        
        if not isinstance(molecule, Molecule):
            raise ValueError('The passed molecule is not a Molecule object')
        
        self.__list[index] = molecule

    # add a new molecule to the molecule list    
        
    def add(self,molecule):        
        
        if not isinstance(molecule, Molecule):
            raise ValueError('The passed molecule is not a Molecule object')
        
            self.__list.append(molecule)

    # This function recovers the total number of molecules currently stored in the database
    def nums(self):
        return len(self.__list)



    def compute_mtx(self, a, b, strict_mtx, loose_mtx, options):
        
        # name = multiprocessing.current_process().name
        # print name
        # print 'a = %d, b = %d' % (a,b)
        # print '\n'        
    
        n = self.nums()
    
        for k in range(a,b+1):
            i = int(n - 2 - math.floor(math.sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5))
            j = int(k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2)
            #print 'k = %d , i = %d , j = %d' % (k,i,j)
    
            moli = self[i].getMolecule()
            molj = self[j].getMolecule()

            #print 'Processing molecules:\n%s\n%s' % (self[i].getName(),self[j].getName())

            # The MCS object between moli and molj is created with the passed option parameters
            try:
                MC = mcs.MCS(moli, molj, options)
            except:
                print 'Skipping....'
                continue
                

            # The scoring between the two molecules is performed by using different rules.
            # The total score will be the product of all the single rules
               
            tmp_scr = MC.ecr() * MC.mncar() * MC.mcsr()

            strict_scr = tmp_scr *  MC.tmcsr(strict_flag=True) 
            loose_scr = tmp_scr * MC.tmcsr(strict_flag=False) 
                
            strict_mtx[k] = strict_scr
            loose_mtx[k] = loose_scr
    
        return


    # This function build the matrix score by using the implemented class MCS (Maximum Common Subgraph)
    def build_matrices(self, options):
        
        print 'Matrix scoring in progress....'   
        
        self.strict_mtx = SMatrix(shape=(self.nums(),))
        self.loose_mtx = SMatrix(shape=(self.nums(),))

        #length linear array input
        l = self.nums()*(self.nums() - 1)/2

        
        if options.parallel == 1:#Serial execution
            self.compute_mtx(0, l-1, self.strict_mtx, self.loose_mtx, options)
        else:#Parallel execution
            
            print 'Parallel is on'
            
            #number of processors
            np = options.parallel

            delta = l/np
            rem = l%np

            if delta < 1:
                kmax = l
            else:
                kmax = np

            proc = []

            strict_mtx = multiprocessing.Array('d', self.strict_mtx)
            loose_mtx =  multiprocessing.Array('d', self.loose_mtx)

            #Chopping the indexes ridistribuiting the remainder
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

                p = multiprocessing.Process(target=self.compute_mtx , args=(i, j, strict_mtx, loose_mtx, options))
                p.start()
                proc.append(p)
                    
            for p in proc:
                p.join()
          
                
            self.strict_mtx = strict_mtx[:]
            self.loose_mtx = loose_mtx[:]
            


class SMatrix(np.ndarray):
    """
    This class implements a "basic" interface for symmetric matrices 
    subclassing ndarray. The class interanlly stores a bidimensional 
    numpy array as a linear array, however the user can still access 
    to the matrix elements by using a two indeces notation
    
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

        #Inizialization
        obj = obj*0.0
        
        return obj
        

    def __getitem__(self, *kargs):

        if isinstance( kargs[0], int ):
            k = kargs[0]
            return super(SMatrix, self).__getitem__(k)

        elif len(kargs[0]) > 2:
            raise ValueError('Tuple dimension must be two')
                
        i = kargs[0][0]
        j = kargs[0][1]

        if i == j:
            return 0.0
        
        l = self.size
        
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
             
        if isinstance( kargs[0], int ):
            k = kargs[0]
            value = kargs[1]
            return super(SMatrix, self).__setitem__(k,value)

        elif len(kargs[0]) > 2:
            raise ValueError('Tuple dimension must be two')
                
        i = kargs[0][0]
        j = kargs[0][1]
        value = kargs[1]

        l = self.size
        
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


class Molecule(object):
    """
    This Class stores the Rdkit molecule objects, their identification number 
    and the total number of instantiated molecules 

    """

    # This variable is used to count the current total number of molecules
    # The variable is defined as private
    __total_molecules = 0

    # Initialization method
    def __init__(self, molecule, molname):
        
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

    # This function returns the molecule ID    
    def getID(self):
        return self.__ID

    # This function returns the copy of the RDkit molecule object
    def getMolecule(self):
        mol_copy = Chem.Mol(self.__molecule)
        return mol_copy


    # This function returns the molecule name
    def getName(self):
        return self.__name


    # This class function returns the current total number of instantiated molecules. 
    @staticmethod
    def get_mol_num():
        return Molecule.__total_molecules


