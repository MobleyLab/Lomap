from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
import math
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem import Draw
from itertools import chain
import sys
import copy


class MCS(object):
    """
    This class is used to compute the Maximum Common Subgraph (MCS) between two RDkit molecule objects.  
    
    """

    def __init__(self, moli, molj, options):
        """
        Inizialization function
    
        moli: Rdkit molecule object related to the first molecule used to perform the MCS calculation
        molj: Rdkit molecule object related to the second molecule used to perform the MCS calculation
        options: MCS options

        """

        def map_mcs_mol(mcs_mol):
            """
            This function is used to define a map between the generated mcs, the starting molecules and 
            vice versa
    
            mcs_mol: the mcs molecule generated from the mcs calulation between the two passed molecules
 
            """
   
            # mcs indexes mapped back to the first molecule (moli) used to perform the MCS
            moli_sub = self.__moli_noh.GetSubstructMatch(mcs_mol)
              
            mcsi_sub = mcs_mol.GetSubstructMatch(mcs_mol)
            
            # mcs to moli
            map_mcs_mol_to_moli_sub = zip(mcsi_sub, moli_sub)

            #print  map_mcs_mol_to_moli_sub
           
            # An RDkit atomic property is defined to store the mapping to moli
            for idx in map_mcs_mol_to_moli_sub:
                mcs_mol.GetAtomWithIdx(idx[0]).SetProp('to_moli',str(idx[1]))

            # mcs indexes mapped back to the second molecule (molj) used to perform the MCS
            molj_sub = self.__molj_noh.GetSubstructMatch(mcs_mol)
              
            mcsj_sub = mcs_mol.GetSubstructMatch(mcs_mol)
             
            # mcs to molj
            map_mcs_mol_to_molj_sub = zip(mcsj_sub, molj_sub)
             
            #print map_mcs_mol_to_molj_sub
            
            # An RDkit atomic property is defined to store the mapping to molj
            for idx in map_mcs_mol_to_molj_sub:
                mcs_mol.GetAtomWithIdx(idx[0]).SetProp('to_molj',str(idx[1]))

            # Chirality
            chiral_at_moli_noh = [seq[0] for seq in Chem.FindMolChiralCenters(self.__moli_noh)]
            chiral_at_molj_noh = [seq[0] for seq in Chem.FindMolChiralCenters(self.__molj_noh)]

            chiral_at_mcs_moli_noh = set([seq[0] for seq in map_mcs_mol_to_moli_sub if seq[1] in chiral_at_moli_noh])
            chiral_at_mcs_molj_noh = set([seq[0] for seq in map_mcs_mol_to_molj_sub if seq[1] in chiral_at_molj_noh])

            chiral_at_mcs = chiral_at_mcs_moli_noh | chiral_at_mcs_molj_noh
            
            for i in chiral_at_mcs:
                at = mcs_mol.GetAtomWithIdx(i)
                at.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)

            #Chem.SanitizeMol(mcs_mol)

            #print chiral_at_mcs

            #For each mcs atom is we save its original index in a specified poperty. This is useful if 
            #will be necessary to remove/add atoms from the mcs molecule
            for at in mcs_mol.GetAtoms():
                at.SetProp('org_idx',str(at.GetIdx()))

            return


        # Local pointers to the passed molecules
        self.moli = moli
        self.molj = molj

        # Local pointers to the passed molecules without hydrogens
        # These variables are defined as private
        self.__moli_noh = AllChem.RemoveHs(moli)
        self.__molj_noh = AllChem.RemoveHs(molj)
        
        # Strict MCS pattern calculation 
        self.__strict_mcs = rdFMCS.FindMCS([self.__moli_noh, self.__molj_noh],
                                           timeout=options.time,
                                           atomCompare=rdFMCS.AtomCompare.CompareAny, 
                                           bondCompare=rdFMCS.BondCompare.CompareAny, 
                                           matchValences=False, 
                                           ringMatchesRingOnly=True, 
                                           completeRingsOnly=True,
                                           matchChiralTag=False)
        # Loose MCS pattern calculation
        self.__loose_mcs = rdFMCS.FindMCS([self.__moli_noh, self.__molj_noh],
                                          timeout=options.time, 
                                          atomCompare=rdFMCS.AtomCompare.CompareAny, 
                                          bondCompare=rdFMCS.BondCompare.CompareAny, 
                                          matchValences=False, 
                                          ringMatchesRingOnly = False, 
                                          completeRingsOnly=False, 
                                          matchChiralTag=False)
                                        
        # Checking
        if self.__strict_mcs.canceled:
            print 'WARNING: timeout reached to find the MCS between molecules: %d and %d' % (self.moli.getID(),self.molj.getID())            
        if self.__loose_mcs.canceled:
            print 'WARNING: timeout reached to find the MCS between molecules: %d and %d' % (self.moli.getID(),self.molj.getID())
        if self.__strict_mcs.numAtoms == 0:
            print 'WARNING: no strict MCS was found between molecules: %d and %d' % (self.moli.getID(),self.molj.getID())
        if self.__loose_mcs.numAtoms == 0:
            print 'WARNING: no loose MCS was found between molecules: %d and %d' % (self.moli.getID(), self.molj.getID())
                                
            
        # The found MCS patterns (smart strings) are converted to RDkit molecule objects
        self.strict_mcs_mol = Chem.MolFromSmarts(self.__strict_mcs.smartsString)
        self.loose_mcs_mol = Chem.MolFromSmarts(self.__loose_mcs.smartsString)

        # Mapping between the found MCS molecules and the passed molecules moli and molj for the strict and loose mcs
        map_mcs_mol(self.strict_mcs_mol)
        map_mcs_mol(self.loose_mcs_mol)


        return


    def draw(self, strict_fname = 'strict_mcs.png', loose_fname = 'loose_mcs.png'):
        """
        This function is used to draw the passed molecules and the strict and loose mcs molecules
        At this stage it is used as debugging tools
                    
        """
   
        
        moli_noh = Chem.Mol(self.__moli_noh)
        molj_noh = Chem.Mol(self.__molj_noh)
        
        strict_mcs_mol = Chem.Mol(self.strict_mcs_mol) 
        loose_mcs_mol = Chem.Mol(self.loose_mcs_mol)

        Chem.SanitizeMol(strict_mcs_mol)
        Chem.SanitizeMol(loose_mcs_mol)
        

        strict_moli_sub = moli_noh.GetSubstructMatch(strict_mcs_mol)
        strict_molj_sub = molj_noh.GetSubstructMatch(strict_mcs_mol)

        strict_mcs_sub =  strict_mcs_mol.GetSubstructMatch(strict_mcs_mol)
        
        AllChem.Compute2DCoords(moli_noh)
        AllChem.Compute2DCoords(molj_noh)
        AllChem.Compute2DCoords(strict_mcs_mol)
        AllChem.Compute2DCoords(loose_mcs_mol)
       

        DrawingOptions.includeAtomNumbers=True
        
        moli_fname='Moli'
        molj_fname='Molj'
        mcs_fname = 'Mcs'

        img = Draw.MolsToGridImage([moli_noh, molj_noh, strict_mcs_mol], 
                                   molsPerRow=3, subImgSize=(200,200),
                                   legends=[moli_fname,molj_fname,mcs_fname], 
                                   highlightAtomLists=[strict_moli_sub, strict_molj_sub, strict_mcs_sub] )

        img.save(strict_fname)

        
             
        loose_moli_sub = moli_noh.GetSubstructMatch(loose_mcs_mol)
        loose_molj_sub = molj_noh.GetSubstructMatch(loose_mcs_mol)
        
        loose_mcs_sub =  loose_mcs_mol.GetSubstructMatch(loose_mcs_mol)
        


        img = Draw.MolsToGridImage([moli_noh, molj_noh, loose_mcs_mol], 
                                   molsPerRow=3, subImgSize=(200,200),
                                   legends=[moli_fname,molj_fname,mcs_fname], 
                                   highlightAtomLists=[loose_moli_sub, loose_molj_sub, loose_mcs_sub] )

        img.save(loose_fname)

        return

    ############ RULES ############

    #ECR Rule (Electrostatic rule)
    def ecr_score(self):
         
        total_charge_moli = 0.0
            
        for atom in self.moli.GetAtoms():
            total_charge_moli += float(atom.GetProp('_TriposPartialCharge'))

        total_charge_molj = 0.0
        for atom in self.molj.GetAtoms():
            total_charge_molj += float(atom.GetProp('_TriposPartialCharge'))

        if abs(total_charge_molj - total_charge_moli) < 1e-3:
            scr_ecr = 1.0
        else:
            scr_ecr = 0.0

            
        return scr_ecr


    # MCSR Rule
    def mcsr(self,mcs_mol, beta=0.1):
        
        # number heavy atoms
        nha_moli = self.moli.GetNumHeavyAtoms()
        nha_molj = self.molj.GetNumHeavyAtoms()
        nha_mcs_mol = mcs_mol.GetNumHeavyAtoms()
            
        scr_mcsr = math.exp(-beta*(nha_moli + nha_molj - 2*nha_mcs_mol))

        return scr_mcsr


    # MNACR rule
    def mncar(self,mcs_mol, ths=3):
            
        nha_mcs_mol = mcs_mol.GetNumHeavyAtoms()
    
        if nha_mcs_mol > ths:
            scr_mncar = 1.0
        else:
            scr_mncar = 0.0 
        
        return scr_mncar


    # TMCRS rule (Trim rule)
    def tmcsr(self, mcs_mol, beta=0.1):
        
        def delete_broken_ring(mcs_mol):
            
            def extend_conflict_to_whole_ring(mol, conflict):
                
                ring_info = mol.GetRingInfo()
                
                rings = ring_info.AtomRings()

                whole_rings = []

                for ring in rings:
                    whole_rings.append(set(ring))

                old_len=0

                while(old_len != len(conflict)):
                    old_len=len(conflict)
                    atom_to_add = set()
                    for at in conflict:
                        for i, ring in enumerate(whole_rings):
                            if at in ring:
                                atom_to_add |= set(ring)
                                whole_rings[i]=[]

                    conflict |= atom_to_add
                
                return
            

            def extend_conflict_to_noaromatic_ring(mol, conflict):
                
                ring_info = mol.GetRingInfo()
                
                rings = ring_info.AtomRings()

                no_aromatic_rings  = []

                for ring in rings:
                    no_aromatic_rings.append(set([j for j in ring if not mol.GetAtomWithIdx(j).GetIsAromatic()]))
                    
                old_len=0

                while(old_len != len(conflict)):
                    old_len=len(conflict)
                    atom_to_add = set()
                    for at in conflict:
                        for i, ring in enumerate(no_aromatic_rings):
                            if at in ring:
                                atom_to_add |= set(ring)
                                no_aromatic_rings[i]=[]

                    conflict |= atom_to_add
                
                return

            #Try to rewrite te following section by using rdkit ring info    
                
            moli_ring_at = set()
            
            for atom in self.__moli_noh.GetAtoms():
                if atom.IsInRing():
                    moli_ring_at.add(atom.GetIdx())

            molj_ring_at = set()
            
            for atom in self.__molj_noh.GetAtoms():
                if atom.IsInRing():
                    molj_ring_at.add(atom.GetIdx())


            mcs_mol_no_ring_at = set()
                    

            for atom in mcs_mol.GetAtoms():
                if not atom.IsInRing():
                    mcs_mol_no_ring_at.add(atom.GetIdx())
            

            moli_conflict = set([int(mcs_mol.GetAtomWithIdx(i).GetProp('to_moli')) for i in mcs_mol_no_ring_at]) & moli_ring_at
            molj_conflict = set([int(mcs_mol.GetAtomWithIdx(i).GetProp('to_molj')) for i in mcs_mol_no_ring_at]) & molj_ring_at

            
            if mcs_mol is self.strict_mcs_mol:
                extend_conflict_to_whole_ring(self.__moli_noh, moli_conflict)
                extend_conflict_to_whole_ring(self.__molj_noh, molj_conflict)
            else:
                extend_conflict_to_noaromatic_ring(self.__moli_noh, moli_conflict)
                extend_conflict_to_noaromatic_ring(self.__molj_noh, molj_conflict)

            #MCS atoms that corresponds to atoms in the moli and molj conflict
            moli_to_mcs_conflict = set()
            molj_to_mcs_conflict = set()

            for at in mcs_mol.GetAtoms():
                idx = int(at.GetProp('to_moli'))
                if idx in moli_conflict:
                   moli_to_mcs_conflict.add(at.GetIdx())
                idx = int(at.GetProp('to_molj'))
                if idx in molj_conflict:
                    molj_to_mcs_conflict.add(at.GetIdx())


            conflict_mcs = list(moli_to_mcs_conflict | molj_to_mcs_conflict)
            
            conflict_mcs.sort(reverse=True)

            #Delete conflict atoms from the copy of mcs_mol

            edit_mcs_mol = Chem.EditableMol(mcs_mol)

            #WARNING: atom indexes are changed
            for i in conflict_mcs:
                edit_mcs_mol.RemoveAtom(i)


            # HERE############################# Try avoid copy
            mcs_mol = edit_mcs_mol.GetMol()
        

            return mcs_mol, conflict_mcs


        #copy of the original mcs
        #mcs_mol_copy = copy.copy(mcs_mol)

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        mcs_mol_copy = Chem.Mol(mcs_mol)

        orig_nha_mcs_mol = mcs_mol_copy.GetNumHeavyAtoms() 

        #At this point the mcs_mol_copy is changed 
        mcs_mol_copy, partial_ring = delete_broken_ring(mcs_mol_copy)


        mcs_ring_set = set()
        mcs_chiral_set = set()

        for atom in mcs_mol_copy.GetAtoms():
                if atom.IsInRing():
                    mcs_ring_set.add(atom.GetIdx())
                if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                    mcs_chiral_set.add(atom.GetIdx())
        
        
        #Loop over the mcs chirial atoms to see if they are also ring atoms
        delete_atoms = set()
        
        for atom_idx in mcs_chiral_set:
            
            if atom_idx in mcs_ring_set:
                
                at = mcs_mol_copy.GetAtomWithIdx(atom_idx)
               
                neighs = at.GetNeighbors()
                neighs_set = set()

                for atom in neighs:
                    neighs_set.add(atom.GetIdx())
        
                delete_atoms |= (neighs_set - mcs_ring_set)

            else:
                #If the chiral atom is not a ring atom, we simple delete it
                delete_atom.add(atom_idx)


        delete_atoms = list(delete_atoms)

        delete_atoms.sort(reverse=True)
    
        edit_mcs_mol = Chem.EditableMol(mcs_mol_copy)

        #WARNING atom indexes are changed
        for idx in delete_atoms:
            edit_mcs_mol.RemoveAtom(idx)

        mcs_mol_copy = edit_mcs_mol.GetMol()
        
        fragments = Chem.rdmolops.GetMolFrags(mcs_mol_copy)

        max_idx = 0
        lgt_max = 0
        
        for idx in range(0,len(fragments)):
            lgt = len(fragments[idx])
            if lgt > lgt_max:
                lgt_max = lgt
                max_idx = idx
    
            
        max_frag = fragments[max_idx]
        
        #The number of heavy atoms in the max fragment
        max_frag_num_heavy_atoms = 0
        for idx in max_frag:
            at = mcs_mol_copy.GetAtomWithIdx(idx)
            if at.GetAtomicNum() > 1:
                max_frag_num_heavy_atoms += 1   


        return math.exp(-2*beta*(orig_nha_mcs_mol - max_frag_num_heavy_atoms))     
        
        
if ("__main__" == __name__) :
    
    mola = Chem.MolFromMol2File('data/chiral2.mol2', sanitize=False, removeHs=False)
    molb = Chem.MolFromMol2File('data/chiral1.mol2', sanitize=False, removeHs=False)
    
    MC = MCS(mola,molb,20)

    MC.draw()

    MC.tmcsr(MC.loose_mcs_mol)
