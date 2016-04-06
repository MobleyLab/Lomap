#******************
# MODULE DOCSTRING
#******************

"""

LOMAP: Maximum Common Subgraph and scoring calculations
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
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem import Draw
import sys
import math
from rdkit import RDLogger
import logging
import argparse


#*******************************
# Maximum Common Subgraph Class
#*******************************


__all__ = ['MCS']

class MCS(object):
    """

    This class is used to compute the Maximum Common Subgraph (MCS) between two
    RDkit molecule objects and to score their similarity by using defined rules 
    
    """

    def __init__(self, moli, molj, options=argparse.Namespace(time=20, verbose='info')):
        """
        Inizialization function
    
        Parameters
        ----------

        moli : RDKit molecule object 
            the first molecule used to perform the MCS calculation
        molj : RDKit molecule object 
            the second molecule used to perform the MCS calculation
        options : argparse python object 
            the list of user options 
       
        """

        def map_mcs_mol():
            """

            This function is used to define a map between the generated mcs, the
            molecules and vice versa
           
            """
   
            # mcs indexes mapped back to the first molecule moli

            if self.__moli_noh.HasSubstructMatch(self.mcs_mol):
                moli_sub = self.__moli_noh.GetSubstructMatch(self.mcs_mol) 
            else:
                raise ValueError('RDkit MCS Subgraph first molecule search failed')
                
            mcsi_sub = self.mcs_mol.GetSubstructMatch(self.mcs_mol)
            
            # mcs to moli
            map_mcs_mol_to_moli_sub = zip(mcsi_sub, moli_sub)

            #print  map_mcs_mol_to_moli_sub
           
            # An RDkit atomic property is defined to store the mapping to moli
            for idx in map_mcs_mol_to_moli_sub:
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp('to_moli', str(idx[1]))

            # mcs indexes mapped back to the second molecule molj 

            if self.__molj_noh.HasSubstructMatch(self.mcs_mol):
                molj_sub = self.__molj_noh.GetSubstructMatch(self.mcs_mol)
            else:
                raise ValueError('RDkit MCS Subgraph second molecule search failed')
             

            if self.mcs_mol.HasSubstructMatch(self.mcs_mol):
                mcsj_sub = self.mcs_mol.GetSubstructMatch(self.mcs_mol)
            else:
                raise ValueError('RDkit MCS Subgraph search failed')
        
   
            # mcs to molj
            map_mcs_mol_to_molj_sub = zip(mcsj_sub, molj_sub)
             
            #print map_mcs_mol_to_molj_sub
            
            # Map between the two molecules
            self.__map_moli_molj = zip(moli_sub, molj_sub)

            # An RDkit atomic property is defined to store the mapping to molj
            for idx in map_mcs_mol_to_molj_sub:
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp('to_molj', str(idx[1]))

            # Chirality

            # moli chiral atoms
            chiral_at_moli_noh = [seq[0] for seq in Chem.FindMolChiralCenters(self.__moli_noh)]
            # molj chiral atoms
            chiral_at_molj_noh = [seq[0] for seq in Chem.FindMolChiralCenters(self.__molj_noh)]

            chiral_at_mcs_moli_noh = set([seq[0] for seq in map_mcs_mol_to_moli_sub if seq[1] in chiral_at_moli_noh])
            chiral_at_mcs_molj_noh = set([seq[0] for seq in map_mcs_mol_to_molj_sub if seq[1] in chiral_at_molj_noh])
            
            # mcs chiral atoms
            chiral_at_mcs = chiral_at_mcs_moli_noh | chiral_at_mcs_molj_noh
            
            for i in chiral_at_mcs:
                at = self.mcs_mol.GetAtomWithIdx(i)
                at.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)


            if chiral_at_mcs and options.verbose == 'pedantic':
                logging.info('Chiral atom detected')

            # For each mcs atom we save its original index in a specified 
            # property. This could be very usefull in the code development
            # when deletition or atom insertions are performed
            for at in self.mcs_mol.GetAtoms():
                at.SetProp('org_idx',str(at.GetIdx()))

            return


        def set_ring_counter(mol):
            
            """

            This function is used to attach to each molecule atom a ring counter
            rc. This parameter is used to asses if a ring has been broken or not
            during the MCS mapping
         
            Parameters
            ----------
            mol : RDKit Molecule obj
                the molecule used to define the atom ring counters
       

            """
            
            # set to zero the atom ring counters
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

        # Set logging level and format
        logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.INFO)
    
        # Local pointers to the passed molecules
        self.moli = moli
        self.molj = molj
        
        if not options.verbose == 'pedantic':
            lg = RDLogger.logger()
            lg.setLevel(RDLogger.CRITICAL)
        
        # Local pointers to the passed molecules without hydrogens
        # These variables are defined as private
        try:
            self.__moli_noh = AllChem.RemoveHs(moli)
            self.__molj_noh = AllChem.RemoveHs(molj)
        except Exception:
            self.__moli_noh = AllChem.RemoveHs(moli, sanitize=False)
            self.__molj_noh = AllChem.RemoveHs(molj, sanitize=False)
            
            Chem.SanitizeMol(self.__moli_noh, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            Chem.SanitizeMol(self.__molj_noh, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            

        # MCS calculaton. In RDKit the MCS is a smart string. Ring atoms are 
        # always mapped in ring atoms. 
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
            logging.warning('Timeout reached to find the MCS between the molecules')
  
        if self.__mcs.numAtoms == 0:
            raise ValueError('No MCS was found between the molecules')
        

        # The found MCS pattern (smart strings) is converted to a RDKit molecule
        self.mcs_mol = Chem.MolFromSmarts(self.__mcs.smartsString)

                
        try: # Try to sanitize the MCS molecule
            Chem.SanitizeMol(self.mcs_mol)
        except Exception: # if not try to recover the atom aromaticity wich is 
            # important for the ring counter
            sanitFail = Chem.SanitizeMol(self.mcs_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY, catchErrors=True)
            if sanitFail: # if not the MCS is skipped
                raise ValueError('Sanitization Failed...')

        # Mapping between the found MCS molecule and moli,  molj
        try:
            map_mcs_mol()
        except Exception as e:
            raise ValueError(str(e))

        #Set the ring counters for each molecule
        set_ring_counter(self.__moli_noh)
        set_ring_counter(self.__molj_noh)
        set_ring_counter(self.mcs_mol)

        # for at in self.mcs_mol.GetAtoms():
        #     print 'at = %d rc = %d' % (at.GetIdx(), int(at.GetProp('rc')))

        if not options.verbose == 'pedantic':
            lg.setLevel(RDLogger.WARNING)
        
        return

    def getMap(self):
        """

        This function is used to return a list of pairs of atom indexes generated
        by the mapping between the two molecules used to calculate the MCS.

        Returns
        -------
        pair of indexes related to the atom mapping 

        """
        
        return self.__map_moli_molj


    def draw_molecule(self, mol, fname='mol.png'):
      
        """

        This function is used to draw a molecule. The main purpose is for debugging
        
        Parameters
        ----------
        mol : RDkit molecule obj 
            the molecule to draw
        fname : string
            the filename used for the .png file
        
        """
  
        DrawingOptions.includeAtomNumbers=True
        AllChem.Compute2DCoords(mol)
        Chem.Draw.MolToFile(mol,fname)
        
        # Useful info for debugging
        # for at in mol.GetAtoms():
        #     print 'atn = %d rc = %d org = %d to_molij = (%d,%d)' \
        #         % (at.GetIdx(), int(at.GetProp('rc')),  
        #            int(at.GetProp('org_idx')),
        #            int(at.GetProp('to_moli')), int(at.GetProp('to_molj')))
        return
        

    def draw_mcs(self, fname='mcs.png', verbose='off'):
        """
        
        This function is used to draw the passed molecules and their mcs
        The main use is for debugging
       
        Parameters
        ----------
        fname : string
            the filename used for the .png file
             
        """

        #Copy of the molecules
        moli_noh = Chem.Mol(self.__moli_noh)
        molj_noh = Chem.Mol(self.__molj_noh)
        mcs_mol = Chem.Mol(self.mcs_mol) 
        
        if not verbose == 'pedantic':
            lg = RDLogger.logger()
            lg.setLevel(RDLogger.CRITICAL)
        
        try:
            Chem.SanitizeMol(self.mcs_mol)
        except Exception: # if not try to recover the atom aromaticity wich is 
            # important for the ring counter
            sanitFail = Chem.SanitizeMol(self.mcs_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY, catchErrors=True)
            if sanitFail: # if not the MCS is skipped
                raise ValueError('Sanitization Failed...')

        if self.__moli_noh.HasSubstructMatch(self.mcs_mol):
            moli_sub = moli_noh.GetSubstructMatch(mcs_mol)
        else:
            raise ValueError('MCS Subgraph moli failed')
            
        if self.__molj_noh.HasSubstructMatch(self.mcs_mol):
            molj_sub = molj_noh.GetSubstructMatch(mcs_mol)
        else:
            raise ValueError('MCS Subgraph molj failed')

        
        if self.mcs_mol.HasSubstructMatch(self.mcs_mol): 
            mcs_sub =  self.mcs_mol.GetSubstructMatch(mcs_mol)
        else:
            raise ValueError('MCS Subgraph failed')

        AllChem.Compute2DCoords(moli_noh)
        AllChem.Compute2DCoords(molj_noh)
        AllChem.Compute2DCoords(mcs_mol)
               
        DrawingOptions.includeAtomNumbers=True
        
        moli_fname='Moli'
        molj_fname='Molj'
        mcs_fname = 'Mcs'

        img = Draw.MolsToGridImage([moli_noh, molj_noh, mcs_mol], 
                                   molsPerRow=3, subImgSize=(400,400),
                                   legends=[moli_fname,molj_fname,mcs_fname], 
                                   highlightAtomLists=[moli_sub, molj_sub, mcs_sub])

        img.save(fname)

        DrawingOptions.includeAtomNumbers=False

        return

    ############ MCS BASED RULES ############

    # MCSR Rule
    def mcsr(self, beta=0.1):
        
        """
        This rule computes the similarity between the two passed molecules 
        used to compute the MCS
        
        Parameters
        ----------
        beta : float
            a parameter used to refine the exponential function used in the
            scoring

        Returns
        -------
        scr_mcsr : float
            the rule score

             
        """

        # The number of heavy atoms in each molecule
        nha_moli = self.moli.GetNumHeavyAtoms()
        nha_molj = self.molj.GetNumHeavyAtoms()
        nha_mcs_mol = self.mcs_mol.GetNumHeavyAtoms()

        # score
        scr_mcsr = math.exp(-beta*(nha_moli + nha_molj - 2*nha_mcs_mol))

        return scr_mcsr


    # MNACR rule
    def mncar(self, ths=4):
 
        """
        This rule cut the similarity score between two molecules if they do
        not share the selected number of atoms 

        
        Parameters
        ----------
        ths : float
            the minumum number of atoms to share
        
        Returns
        -------
        scr_mncar : float
            the rule score     
        """
       
        # This rule has been modified from the rule desribed in the Lomap paper
        # to match the LOMAP first implementation provided by schrodinger
 
        nha_mcs_mol = self.mcs_mol.GetNumHeavyAtoms()
        nha_moli = self.moli.GetNumHeavyAtoms()
        nha_molj = self.molj.GetNumHeavyAtoms()
    
        scr_mncar = float((nha_mcs_mol >= ths) or (nha_moli + 3) or (nha_molj + 3))
     
        return scr_mncar


    # TMCRS rule (Trim rule) 
    def tmcsr(self, beta=0.1, strict_flag=True):
       
        """
        This rule check if rings have been broken during the MCS mapping 
        and if chiral atoms are presents. If rings are broken all the 
        remaining ring atoms are deleted. Atoms connected to chiral centers
        are deleted as well
        
 
        Parameters
        ----------
        beta : float
            a parameter used to refine the exponential function used 
            in the scoring
            
        stric_flag : bool
            a flag used to select the scrict or loose mode
             
        """
       
        def delete_broken_ring():

            # Strict: we cancel all the atoms in conflict in the mcs and 
            # delete all eventually non ring atoms that are left 
            def extend_conflict(mol, conflict):
                """
            
                This function check if rings have been broken during the MCS mapping
                deleting all the remaining atom rings. In strict mode all the 
                conflicting ring atoms are deleted. In loose mode only non planar
                atom rings are deleted
 

                Parameters
                ----------
                mol : RDKit molecule obj
                    the mcs molecule
                conflict : set
                    the set of atoms in Moli and Molj that are in conflict with 
                    the MCS molecule. A conflict is generated if the ring counter
                    between the MCS and Moli/Molj changes

                
                Returns
                -------
                mcs_mol : RDKit molecule obj
                    a copy of the edited mcs molecule
                       
                """
     
                mcs_conflict = list(conflict)
                mcs_conflict.sort(reverse=True)


                # Editing the mcs molecule deleting all the selected conficting atoms
                edit_mcs_mol = Chem.EditableMol(mol)

                # WARNING: atom indexes are changed
                for i in mcs_conflict:
                    edit_mcs_mol.RemoveAtom(i) 
                
                mcs_mol = edit_mcs_mol.GetMol()
              
                # The mcs molecule could be empty at this point
                if not mcs_mol.GetNumAtoms():
                    return mcs_mol
                
                # Deleting broken ring atoms if the atom rc > 0 and the atom is not 
                # in a ring anymore
                mcs_conflict = [ at.GetIdx()  for at in mcs_mol.GetAtoms() if int(at.GetProp('rc')) > 0 and not at.IsInRing()]
                
                mcs_conflict.sort(reverse=True)

                edit_mcs_mol = Chem.EditableMol(mcs_mol)
                
                # WARNING: atom indexes are changed
                for i in mcs_conflict:
                    edit_mcs_mol.RemoveAtom(i) 
                    
                mcs_mol = edit_mcs_mol.GetMol()

                # The mcs molecule could be empty at this point
                if not mcs_mol.GetNumAtoms():
                    return mcs_mol

                # Deleting eventually disconnected parts and keep the max fragment left
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

                # WARNING: atom indexes have changed
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
                
                # Moli atom is a ring atom (rc>0) and its rc is different from 
                # the corresponding mcs rc atom  
                if moli_idx_rc > 0 and (moli_idx_rc != int(at.GetProp('rc'))):
                    if strict_flag: # In strict mode we add the atom
                        mcs_conflict.add(at.GetIdx())
                    else: # In loose mode we add the atom if it is not an aromatic atom
                        if not at.GetIsAromatic():
                            mcs_conflict.add(at.GetIdx())
                        

                # Molj atom is a ring atom (rc>0) and its rc is different 
                # from the corresponding mcs rc atom 
                if molj_idx_rc > 0 and (molj_idx_rc  != int(at.GetProp('rc'))):
                    if strict_flag: # In strict mode we add the atom
                        mcs_conflict.add(at.GetIdx())
                    else: # In loose mode we add the atom if it is not an aromatic atom
                        if not at.GetIsAromatic():
                            mcs_conflict.add(at.GetIdx())

            mcs_mol = extend_conflict(self.mcs_mol, mcs_conflict)
            
            return mcs_mol


        mcs_mol_copy = Chem.Mol(self.mcs_mol)

        orig_nha_mcs_mol = mcs_mol_copy.GetNumHeavyAtoms() 


        # At this point the mcs_mol_copy has changed 
        mcs_mol_copy = delete_broken_ring()

        # The mcs molecule could be empty at this point
        if not mcs_mol_copy.GetNumAtoms():
            return math.exp(-2*beta*(orig_nha_mcs_mol))

        ### Chiral Atoms ###

        # Deleting Chiral Atoms
        mcs_ring_set = set()
        mcs_chiral_set = set()

        for atom in mcs_mol_copy.GetAtoms():
                if atom.IsInRing():
                    mcs_ring_set.add(atom.GetIdx())
                if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                    mcs_chiral_set.add(atom.GetIdx())
        
        
        # Loop over the mcs chirial atoms to check if they are also ring atoms
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
                # If the chiral atom is not a ring atom, we delete it
                delete_atoms.add(atom_idx)


        delete_atoms = list(delete_atoms)

        delete_atoms.sort(reverse=True)
    
        edit_mcs_mol = Chem.EditableMol(mcs_mol_copy)

        # WARNING atom indexes have changed
        for idx in delete_atoms:
            edit_mcs_mol.RemoveAtom(idx)

        mcs_mol_copy = edit_mcs_mol.GetMol()

        
        # The mcs molecule could be empty at this point
        if not mcs_mol_copy.GetNumAtoms():
            return math.exp(-2*beta*(orig_nha_mcs_mol))


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
        
        # The number of heavy atoms in the max fragment
        max_frag_num_heavy_atoms = 0
        for idx in max_frag:
            at = mcs_mol_copy.GetAtomWithIdx(idx)
            if at.GetAtomicNum() > 1:
                max_frag_num_heavy_atoms += 1   


        return math.exp(-2*beta*(orig_nha_mcs_mol - max_frag_num_heavy_atoms))
        

if ("__main__" == __name__) :

    mola = Chem.MolFromMol2File('test/basic/2-methylnaphthalene.mol2', sanitize=False, removeHs=False)    
    molb = Chem.MolFromMol2File('test/basic/2-naftanol.mol2', sanitize=False, removeHs=False)

    # MCS calculation
    try:
        MC = MCS(mola,molb)
    except Exception:
        raise ValueError('NO MCS FOUND......')
       
    # Mapping between the passed molecules    
    mcs_map = MC.getMap()
    
    print mcs_map

    # Draw the molecules andd their MCS
    MC.draw_mcs()
    
    # Rules calculations
    mcsr = MC.mcsr()
    mncar =  MC.mncar()
    
    strict = MC.tmcsr(strict_flag=True)
    loose = MC.tmcsr(strict_flag=False)

    print 'TMCRS STRICT = %f , TMCRS LOOSE = %f' % (strict, loose)
    print 'MCSR = ', mcsr
    print 'MNCAR = ', mncar
    
    tmp = mcsr * mncar 
    
    print 'Total Strict = %f , Total Loose = %f' % (tmp * strict, tmp * loose)  
