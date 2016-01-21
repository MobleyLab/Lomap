from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem import Draw
from optparse import OptionParser
import sys
import math

class MCS(object):
    """
    This class is used to compute the Maximum Common Subgraph (MCS) between two
    RDkit molecule objects.  
    
    """

    def __init__(self, moli, molj, options):
        """
        Inizialization function
    
        moli: Rdkit molecule object related to the first molecule used to 
        perform the MCS calculation
        molj: Rdkit molecule object related to the second molecule used to 
        perform the MCS calculation
        options: MCS options

        """

        def map_mcs_mol():
            """
            This function is used to define a map between the generated mcs, the
            starting molecules and vice versa
    
            mcs_mol: the mcs molecule generated from the mcs calulation between 
            the two passed molecules
 
            """
   
            # mcs indexes mapped back to the first molecule moli
            moli_sub = self.__moli_noh.GetSubstructMatch(self.mcs_mol)
              
            mcsi_sub = self.mcs_mol.GetSubstructMatch(self.mcs_mol)
            
            # mcs to moli
            map_mcs_mol_to_moli_sub = zip(mcsi_sub, moli_sub)

            #print  map_mcs_mol_to_moli_sub
           
            # An RDkit atomic property is defined to store the mapping to moli
            for idx in map_mcs_mol_to_moli_sub:
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp('to_moli', str(idx[1]))

            # mcs indexes mapped back to the second molecule molj 
            molj_sub = self.__molj_noh.GetSubstructMatch(self.mcs_mol)
              
            mcsj_sub = self.mcs_mol.GetSubstructMatch(self.mcs_mol)
             
            # mcs to molj
            map_mcs_mol_to_molj_sub = zip(mcsj_sub, molj_sub)
             
            #print map_mcs_mol_to_molj_sub
            
            # An RDkit atomic property is defined to store the mapping to molj
            for idx in map_mcs_mol_to_molj_sub:
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp('to_molj', str(idx[1]))

            # Chirality
            chiral_at_moli_noh = [seq[0] for seq in Chem.FindMolChiralCenters(self.__moli_noh)]
            chiral_at_molj_noh = [seq[0] for seq in Chem.FindMolChiralCenters(self.__molj_noh)]

            chiral_at_mcs_moli_noh = set([seq[0] for seq in map_mcs_mol_to_moli_sub if seq[1] in chiral_at_moli_noh])
            chiral_at_mcs_molj_noh = set([seq[0] for seq in map_mcs_mol_to_molj_sub if seq[1] in chiral_at_molj_noh])

            chiral_at_mcs = chiral_at_mcs_moli_noh | chiral_at_mcs_molj_noh
            
            for i in chiral_at_mcs:
                at = self.mcs_mol.GetAtomWithIdx(i)
                at.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)


            if chiral_at_mcs:
                print 'WARNING: Chiral atoms detected'

            #For each mcs atom we save its original index in a specified 
            #property. This could be very usefull in the code development
            for at in self.mcs_mol.GetAtoms():
                at.SetProp('org_idx',str(at.GetIdx()))


            return


        def set_ring_counter(mol):
            
            for at in mol.GetAtoms():
                at.SetProp('rc','0')

            rginfo = mol.GetRingInfo()

            rgs = rginfo.AtomRings()
         
            #print rgs
   
            rgs_set = set([e for l in rgs for e in l])
            
            for idx in rgs_set:
                for r in rgs:
                    if(idx in r):
                        val = int(mol.GetAtomWithIdx(idx).GetProp('rc'))
                        val = val + 1
                        mol.GetAtomWithIdx(idx).SetProp('rc',str(val))
            return
            

        # Local pointers to the passed molecules
        self.moli = moli
        self.molj = molj

        # Local pointers to the passed molecules without hydrogens
        # These variables are defined as private
        self.__moli_noh = AllChem.RemoveHs(moli)
        self.__molj_noh = AllChem.RemoveHs(molj)
        

        # MCS pattern calculation
        self.__mcs = rdFMCS.FindMCS([self.__moli_noh, self.__molj_noh],
                                          timeout=options.time, 
                                          atomCompare=rdFMCS.AtomCompare.CompareAny, 
                                          bondCompare=rdFMCS.BondCompare.CompareAny, 
                                          matchValences=False, 
                                          ringMatchesRingOnly=True, 
                                          completeRingsOnly=False, 
                                          matchChiralTag=False)
                                        
        # Checking
        if self.__mcs.canceled:
            print 'WARNING: timeout reached to find the MCS between molecules: %d and %d' \
            % (self.moli.getID(),self.molj.getID())            
        if self.__mcs.numAtoms == 0:
            print 'WARNING: no strict MCS was found between molecules: %d and %d' \
            % (self.moli.getID(),self.molj.getID())
        

        # The found MCS pattern (smart strings) is converted to a RDkit molecule
        self.mcs_mol = Chem.MolFromSmarts(self.__mcs.smartsString)



        #Sanitize the MCS molecule
        try:
            Chem.SanitizeMol(self.mcs_mol)
        except:
            print 'Sanitization Failed...'
            sanitFail = Chem.SanitizeMol(self.mcs_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY, catchErrors=True)
            if sanitFail: 
                raise ValueError(sanitFail)

        
        # Mapping between the found MCS molecule and moli,  molj
        map_mcs_mol()

        #Set the ring counters for each molecule
        set_ring_counter(self.__moli_noh)
        set_ring_counter(self.__molj_noh)
        set_ring_counter(self.mcs_mol)

        
        # for at in self.mcs_mol.GetAtoms():
        #     print 'at = %d rc = %d' % (at.GetIdx(), int(at.GetProp('rc')))


        return

    def draw_molecule(self,mol,fname='mol.png'):
        
        DrawingOptions.includeAtomNumbers=True
        AllChem.Compute2DCoords(mol)
        Chem.Draw.MolToFile(mol,fname)
        
        for at in mol.GetAtoms():
            print 'atn = %d rc = %d org = %d to_molij = (%d,%d)' \
                % (at.GetIdx(), int(at.GetProp('rc')),  
                   int(at.GetProp('org_idx')),
                   int(at.GetProp('to_moli')), int(at.GetProp('to_molj')))
        return
        


    def draw_mcs(self, fname = 'mcs.png'):
        """
        This function is used to draw the passed molecules and their mcs molecule
        At this stage it is used as debugging tools
                    
        """
   
        #Copy of the molecules
        moli_noh = Chem.Mol(self.__moli_noh)
        molj_noh = Chem.Mol(self.__molj_noh)
        mcs_mol = Chem.Mol(self.mcs_mol) 
        

        try:
            Chem.SanitizeMol(self.mcs_mol)
        except ValueError:
            print "Sanitization failed...."
        

        moli_sub = moli_noh.GetSubstructMatch(self.mcs_mol)
        molj_sub = molj_noh.GetSubstructMatch(self.mcs_mol)

        mcs_sub =  self.mcs_mol.GetSubstructMatch(self.mcs_mol)
        
        AllChem.Compute2DCoords(moli_noh)
        AllChem.Compute2DCoords(molj_noh)
        AllChem.Compute2DCoords(self.mcs_mol)
               
        DrawingOptions.includeAtomNumbers=True
        
        moli_fname='Moli'
        molj_fname='Molj'
        mcs_fname = 'Mcs'

        img = Draw.MolsToGridImage([moli_noh, molj_noh, self.mcs_mol], 
                                   molsPerRow=3, subImgSize=(400,400),
                                   legends=[moli_fname,molj_fname,mcs_fname], 
                                   highlightAtomLists=[moli_sub, molj_sub, mcs_sub] )

        img.save(fname)

        
        return

    ############ RULES ############

    #ECR Rule (Electrostatic rule)
    def ecr(self):
         
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
    def mcsr(self, beta=0.1):
        
        # number heavy atoms
        nha_moli = self.moli.GetNumHeavyAtoms()
        nha_molj = self.molj.GetNumHeavyAtoms()
        nha_mcs_mol = self.mcs_mol.GetNumHeavyAtoms()
            
        scr_mcsr = math.exp(-beta*(nha_moli + nha_molj - 2*nha_mcs_mol))

        return scr_mcsr


    # MNACR rule
    def mncar(self, ths=4):
        
        #This rule has been modified from the rule desribed in the Lomap paper
        #to match the implemented version
 
        nha_mcs_mol = self.mcs_mol.GetNumHeavyAtoms()
        nha_moli = self.moli.GetNumHeavyAtoms()
        nha_molj = self.molj.GetNumHeavyAtoms()
    
        scr_mncar = float((nha_mcs_mol >= ths) or (nha_moli + 3) or (nha_molj + 3))
     
        return scr_mncar


    # TMCRS rule (Trim rule) 
    def tmcsr(self, beta=0.1, strict_flag=True):
        
        def delete_broken_ring():

            #Strict: we cancel all the atoms in conflict in the mcs and 
            #delete all eventually non ring atoms that are left 
            def extend_conflict(mol, conflict):
                
                mcs_conflict = list(conflict)
                mcs_conflict.sort(reverse=True)


                #Editing the mcs molecule deleting all the selected conficting atoms
                edit_mcs_mol = Chem.EditableMol(mol)

                #WARNING: atom indexes are changed
                for i in mcs_conflict:
                    edit_mcs_mol.RemoveAtom(i) 
                
                mcs_mol = edit_mcs_mol.GetMol()

                #Deleting broken ring atoms if the atom rc > 0 and the atom is not 
                #in a ring anymore
                mcs_conflict = [ at.GetIdx()  for at in mcs_mol.GetAtoms() if int(at.GetProp('rc')) > 0 and not at.IsInRing()]
                
                mcs_conflict.sort(reverse=True)

                edit_mcs_mol = Chem.EditableMol(mcs_mol)
                
                #WARNING: atom indexes are changed
                for i in mcs_conflict:
                    edit_mcs_mol.RemoveAtom(i) 
                    
                mcs_mol = edit_mcs_mol.GetMol()


                #Deleting eventually disconnected parts and keep the max fragment left
                fragments = Chem.rdmolops.GetMolFrags(mcs_mol)

                max_idx = 0
                lgt_max = 0
        
                for idx in range(0,len(fragments)):
                    lgt = len(fragments[idx])
                    if lgt > lgt_max:
                        lgt_max = lgt
                        max_idx = idx
    
                        
                max_frag = fragments[max_idx]

                mcs_conflict = [ at.GetIdx() for at in mcs_mol.GetAtoms() if not at.GetIdx() in max_frag ]
        
                mcs_conflict.sort(reverse=True)

                edit_mcs_mol = Chem.EditableMol(mcs_mol)

                #WARNING: atom indexes have changed
                for i in mcs_conflict:
                    edit_mcs_mol.RemoveAtom(i) 
                    
                mcs_mol = edit_mcs_mol.GetMol()

                #self.draw_molecule(mcs_mol)

                return mcs_mol
                
            
            mcs_conflict = set()
            
            for at in self.mcs_mol.GetAtoms():

                moli_idx = int(at.GetProp('to_moli'))
                molj_idx = int(at.GetProp('to_molj'))

                moli_idx_rc =  int(self.__moli_noh.GetAtomWithIdx(moli_idx).GetProp('rc'))
                molj_idx_rc =  int(self.__molj_noh.GetAtomWithIdx(molj_idx).GetProp('rc'))
                
                #Moli atom is a ring atom (rc>0) and its rc is different from 
                #the corresponding mcs rc atom  
                if moli_idx_rc > 0 and (moli_idx_rc != int(at.GetProp('rc'))):
                    if strict_flag:#In strict mode we add the atom
                        mcs_conflict.add(at.GetIdx())
                    else:#In loose mode we add the atom if it is a not aromatic atom only
                        if not at.GetIsAromatic():
                            mcs_conflict.add(at.GetIdx())
                        

                #Molj atom is a ring atom (rc>0) and its rc is different 
                #from the corresponding mcs rc atom 
                if molj_idx_rc > 0 and (molj_idx_rc  != int(at.GetProp('rc'))):
                    if strict_flag:#In strict mode we add the atom
                        mcs_conflict.add(at.GetIdx())
                    else:#In loose mode we add the atom if it is a not aromatic atom only
                        if not at.GetIsAromatic():
                            mcs_conflict.add(at.GetIdx())

            mcs_mol = extend_conflict(self.mcs_mol, mcs_conflict)
            
                        
            return mcs_mol


        mcs_mol_copy = Chem.Mol(self.mcs_mol)

        orig_nha_mcs_mol = mcs_mol_copy.GetNumHeavyAtoms() 


        #At this point the mcs_mol_copy has changed 
        mcs_mol_copy = delete_broken_ring()


        #Deleting Chiral Atoms
        mcs_ring_set = set()
        mcs_chiral_set = set()

        for atom in mcs_mol_copy.GetAtoms():
                if atom.IsInRing():
                    mcs_ring_set.add(atom.GetIdx())
                if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                    mcs_chiral_set.add(atom.GetIdx())
        
        
        #Loop over the mcs chirial atoms to check if they are also ring atoms
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
                #If the chiral atom is not a ring atom, we delete it
                delete_atom.add(atom_idx)


        delete_atoms = list(delete_atoms)

        delete_atoms.sort(reverse=True)
    
        edit_mcs_mol = Chem.EditableMol(mcs_mol_copy)

        #WARNING atom indexes have changed
        for idx in delete_atoms:
            edit_mcs_mol.RemoveAtom(idx)

        mcs_mol_copy = edit_mcs_mol.GetMol()


        #self.draw_molecule(mcs_mol_copy)

        
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
   
    parser = OptionParser( usage = "Usage: %prog [options] <structure-file-dir>", version = "%prog v0.0" )
    parser.add_option("-t", "--time", default = 20 , help = " Set the maximum time to perform the mcs search between pair of molecules")
    
    #A tuple of options and arguments passed by the user
    (opt, args) = parser.parse_args()


    mola = Chem.MolFromMol2File('data/mol1.mol2', sanitize=False, removeHs=False)
    molb = Chem.MolFromMol2File('data/mol2.mol2', sanitize=False, removeHs=False)
    
    MC = MCS(mola,molb,opt)

    MC.draw_mcs()
    
    print 'TMCRS STRICT = %f , TMCRS LOOSE = %f' % (MC.tmcsr(strict_flag=True), MC.tmcsr(strict_flag=False))
    print 'MCSR = ', MC.mcsr()
    print 'MNCAR = ', MC.mncar()
    print 'ECR = ', MC.ecr()
    
    
