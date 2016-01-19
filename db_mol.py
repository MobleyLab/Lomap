import copy
import math
import sys
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


class DB_Molecules(object):
    
    def __init__(self,molecules):
        
        if not isinstance(molecules, list):
            raise ValueError('The passed parents must be a list')

        for mol in molecules:
             if not isinstance(mol , Molecule) :
                 raise ValueError('The passed molecule list does not contain all Molecule objects')
        
        # list container used to store the loaded molecules
        self.__list = molecules

        # index used to perform index selection by using __iter__
        self.__ci = 0

    # index generator
    def __iter__(self):
        return self

    # select the molecule duruning an iteration on the molecules
    def next(self):
        if self.__ci > len(self.__list) - 1:
            self.__ci = 0
            raise StopIteration
        else:
            self.__ci = self.__ci + 1
            return self.__list[self.__ci - 1] 
            
    # slicing and index selection funtion
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

    # the total number of molecules
    def nums(self):
        return len(self.__list)

    # score the molecules
    def score(self, options):
       
        #################################### RULES ##########################################

        # ECR rule
        def ecr(moli,molj):
             
            total_charge_moli = 0.0
            
            for atom in moli.GetAtoms():
                total_charge_moli += float(atom.GetProp('_TriposPartialCharge'))

            total_charge_molj = 0.0
            for atom in molj.GetAtoms():
                total_charge_molj += float(atom.GetProp('_TriposPartialCharge'))

            if abs(total_charge_molj - total_charge_moli) < 1e-3:
                scr_ecr = 1.0
            else:
                scr_ecr = 0.0

            
            return scr_ecr


        # MCSR rule
        def mcsr(moli, molj, mcs_patt, beta=0.1):
            
            mcs_mol = Chem.MolFromSmarts(mcs_patt.smartsString) 
        
            # number heavy atoms
            nha_moli = moli.GetNumHeavyAtoms()
            nha_molj = molj.GetNumHeavyAtoms()
            nha_mcs_mol = mcs_mol.GetNumHeavyAtoms()
            
            scr_mcsr = math.exp(-beta*(nha_moli + nha_molj - 2*nha_mcs_mol))

            return scr_mcsr

        # MNACR rule
        def mncar(mcs_patt, ths=3):
            
            mcs_mol = Chem.MolFromSmarts(mcs_patt.smartsString) 
                
            nha_mcs_mol = mcs_mol.GetNumHeavyAtoms()
            
            if nha_mcs_mol > ths:
                scr_mncar = 1.0
            else:
                scr_mncar = 0.0 
                
            return scr_mncar


        def tmcsr(moli, molj, mcs_pat, filename):
            
            mcs_mol = Chem.MolFromSmarts(mcs_pat.smartsString) 

            orig_nha_mcs_mol = mcs_mol.GetNumHeavyAtoms() 

            moli_sub = moli.GetSubstructMatch(mcs_mol)
            molj_sub = molj.GetSubstructMatch(mcs_mol)
            mcs_mol_sub = mcs_mol.GetSubstructMatch(mcs_mol)

            map_mcs_mol_to_moli_sub = dict(zip(mcs_mol_sub,moli_sub))


            moli_ring_set = set()
            
            for atom in moli.GetAtoms():
                 if atom.IsInRing():
                     moli_ring_set.add(atom.GetIdx())

            



                        

            
            # map_list = zip(moli_sub, molj_sub)

            
            # sys.exit(-1)

            # Allchem.Compute2DCoords(moli)
            # Allchem.Compute2DCoords(molj)

            # from rdkit.Chem.Draw.MolDrawing import DrawingOptions
            
            # DrawingOptions.includeAtomNumbers=True

            # img = Draw.MolsToGridImage([moli,molj], molsPerRow=2, subImgSize=(200,200),
            #                            legends=['b','a'], highlightAtomLists=[moli_sub, molj_sub] )
            
            # img.save(filename)


            # print map_list
            
            
        
        
        #######################################################################


        print 'Scoring.....'   

        # parameters used in different rules

        for i in range(0,self.nums()-1):
            for j in range(i+1,self.nums()):
                
                moli = self[i].getMolecule()
                molj = self[j].getMolecule()

                moli_noh = AllChem.RemoveHs(moli)
                molj_noh = AllChem.RemoveHs(molj)
                
                
                strict_mcs = rdFMCS.FindMCS([moli_noh, molj_noh],
                                        timeout=options.time, 
                                        atomCompare=rdFMCS.AtomCompare.CompareAny, 
                                        bondCompare=rdFMCS.BondCompare.CompareAny, 
                                        matchValences=False, 
                                        ringMatchesRingOnly=True, 
                                        completeRingsOnly=True)

                loose_mcs = rdFMCS.FindMCS([moli_noh, molj_noh],
                                       timeout=options.time, 
                                       atomCompare=rdFMCS.AtomCompare.CompareAny, 
                                       bondCompare=rdFMCS.BondCompare.CompareAny, 
                                       matchValences=False, 
                                       ringMatchesRingOnly = False, 
                                       completeRingsOnly=False )
                

                if strict_mcs.canceled:
                    print 'WARNING: timeout reached to find the MCS between molecules: %d and %d' % (moli.getID(),molj.getID())            
                if loose_mcs.canceled:
                    print 'WARNING: timeout reached to find the MCS between molecules: %d and %d' % (moli.getID(),molj.getID())
                if strict_mcs.numAtoms == 0:
                    print 'WARNING: no strict MCS was found between molecules: %d and %d' % (moli.getID(),molj.getID())
                if loose_mcs.numAtoms == 0:
                    print 'WARNING: no loose MCS was found between molecules: %d and %d' % (moli.getID(),molj.getID())
                    

                ecr_score = ecr(moli,molj)
                
                mcsr_score_strict = mcsr(moli_noh ,molj_noh, strict_mcs)
                mcsr_score_loose = mcsr(moli_noh, molj_noh, loose_mcs)
                            
                mnca_score_strict = mncar(strict_mcs)
                mnca_score_loose = mncar(loose_mcs)

                
                #tmcsr_score_strict = tmcsr(moli_noh, molj_noh, strict_mcs, 'strict.png')
                tmcsr_score_loose = tmcsr(moli_noh, molj_noh, loose_mcs, 'loose.png')
                


                # print ecr_score, (mcsr_score_strict, mcsr_score_loose), (mnca_score_strict, mnca_score_loose)
                

                sys.exit(-1)


class Molecule(object):
    # This variable is used to count the current total number of molecules
    # The variable is defined as private
    __total_molecules = 0

    # Initialization method
    def __init__(self,molecule):
        
        #Check Inputs
        if not isinstance(molecule, Chem.rdchem.Mol):
            raise ValueError('The passed molecule object is not a RdKit molecule')

        # The variable _molecule saves the current RDkit molecule object
        # The variable is defined as private
        self.__molecule = molecule    

            
        # The variable _ID saves the molecule identification number 
        # The variable is defined as private
        self.__ID = Molecule.__total_molecules
    
        Molecule.__total_molecules+=1

    # This function returns the molecule ID    
    def getID(self):
        return self.__ID

    # This function returns the copy of the RDkit molecule object
    def getMolecule(self):
        return copy.copy(self.__molecule)

    @staticmethod
    def get_mol_num():
        return Molecule.__total_molecules


