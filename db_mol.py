import copy
import math
import sys
import mcs
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


class DB_Molecules(object):
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

    
    # This function build the matrix score by using the implemented class MCS (Maximum Common Subgraph)
    def build_matrix(self, options):
        
        print 'Matrix scoring in progress....'   

        for i in range(0,self.nums()-1):
            for j in range(i+1,self.nums()):
                
                # moli and molj are effectively RDKit molecule objects
                moli = self[i].getMolecule()
                molj = self[j].getMolecule()

                # The MCS object between moli and molj is created with the passed option parameters
                MC = mcs.MCS(moli, molj, options)

                # This function is for debugging purpose only (may be useful elsewhere?)
                MC.draw()

                # The scoring between the two molecules is performed by using different rules.
                # The total score will be the product of all the single rules
                scr_strict = MC.ecr_score() * MC.tmcsr(MC.strict_mcs_mol) * MC.mcsr(MC.strict_mcs_mol) * MC.mncar(MC.strict_mcs_mol)
                scr_loose =  MC.ecr_score() * MC.tmcsr(MC.loose_mcs_mol) * MC.mcsr(MC.loose_mcs_mol) * MC.mncar(MC.loose_mcs_mol)

                print scr_strict
                print scr_loose


class Molecule(object):
    """
    This Class stores the Rdkit molecule objects, their identification number and the total number of molecules loaded so far.  

    """

    # This variable is used to count the current total number of molecules
    # The variable is defined as private
    __total_molecules = 0

    # Initialization method
    def __init__(self,molecule):
        
        #Check Inputs
        if not isinstance(molecule, Chem.rdchem.Mol):
            raise ValueError('The passed molecule object is not a RdKit molecule')

        # The variable __molecule saves the current RDkit molecule object
        # The variable is defined as private
        self.__molecule = molecule    

            
        # The variable __ID saves the molecule identification number 
        # The variable is defined as private
        self.__ID = Molecule.__total_molecules
    
        Molecule.__total_molecules+=1

    # This function returns the molecule ID    
    def getID(self):
        return self.__ID

    # This function returns the copy of the RDkit molecule object
    def getMolecule(self):
        return copy.copy(self.__molecule)


    # This class function returns the current total number of molecules. 
    @staticmethod
    def get_mol_num():
        return Molecule.__total_molecules


