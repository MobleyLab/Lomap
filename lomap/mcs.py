# ******************
# MODULE DOCSTRING
# ******************

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

# *****************************************************************************
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
# *****************************************************************************


# ****************
# MODULE IMPORTS
# ****************


from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Geometry.rdGeometry import Point3D
import sys
import math
from rdkit import RDLogger
import logging
import argparse

# *******************************
# Maximum Common Subgraph Class
# *******************************


__all__ = ['MCS']


class MCS(object):
    """

    This class is used to compute the Maximum Common Subgraph (MCS) between two
    RDkit molecule objects and to score their similarity by using defined rules 
    
    """

    def __init__(self, moli, molj, options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=False)):
        """
        Initialization function
    
        Parameters
        ----------

        moli : RDKit molecule object 
            the first molecule used to perform the MCS calculation
        molj : RDKit molecule object 
            the second molecule used to perform the MCS calculation
        options : argparse python object 
            the list of user options 
       
        """

        def substructure_centre(mol, mol_sub):
            """

            This function takes a molecule and a list of atom indices
            in that molecule and returns an RDKit Point3D representing
            the geometric centre of the atoms in the list

            """

            sum = Point3D()
            for i in mol_sub:
                sum += mol.GetConformer().GetAtomPosition(i)
            return sum / len(mol_sub)


        def best_substruct_match(moli,molj,by_rmsd=True):
            """

            This function looks over all of the substructure matches and returns the one
            with the best 3D correspondence (if by_rmsd is true), or the fewest number
            of atomic number mismatches (if by_rmsd is false)

            Note that the 3D correspondence does a translational centreing (but
            does not rotate).

            """

            # Sanity checking
            if not moli.HasSubstructMatch(self.mcs_mol):
                raise ValueError('RDkit MCS Subgraph first molecule search failed')

            if not molj.HasSubstructMatch(self.mcs_mol):
                raise ValueError('RDkit MCS Subgraph second molecule search failed')

            moli_sub = moli.GetSubstructMatches(self.mcs_mol,uniquify=False)
            molj_sub = molj.GetSubstructMatches(self.mcs_mol,uniquify=False)
            best_rmsd=1e8
            for mapi in moli_sub:
                for mapj in molj_sub:
                    # Compute the translation to bring molj's centre over moli
                    coord_delta = (substructure_centre(moli,mapi)
                                 - substructure_centre(molj,mapj))
                    rmsd=0
                    for pair in zip(mapi,mapj):
                        if by_rmsd:
                            rmsd += (moli.GetConformer().GetAtomPosition(pair[0]) 
                                   - molj.GetConformer().GetAtomPosition(pair[1])
                                   - coord_delta).LengthSq()
                        elif (moli.GetAtomWithIdx(pair[0]).GetAtomicNum() != 
                              molj.GetAtomWithIdx(pair[1]).GetAtomicNum()):
                            rmsd+=1
                    if rmsd < best_rmsd:
                        besti=mapi
                        bestj=mapj
                        best_rmsd=rmsd

            return (besti,bestj)

        def trim_mcs_mol(max_deviation=2.0):
            """

            This function is used to trim the MCS molecule to remove mismatched atoms i.e atoms
            where the topological mapping does not work in 3D coordinates.

            The sets of mapped atoms are translated to bring their geometric centres
            into alignment before trimming
           
            Parameters
            ----------

            max_deviation : the maximum difference in Angstroms between mapped atoms to allow

            """

            while True:
                (mapi,mapj) = best_substruct_match(self.__moli_noh,self.__molj_noh,by_rmsd=True)
                # Compute the translation to bring molj's centre over moli
                coord_delta = (substructure_centre(self.__moli_noh,mapi)
                             - substructure_centre(self.__molj_noh,mapj))
                worstatomidx=-1
                worstdist=0
                atomidx=0
                for pair in zip(mapi,mapj):
                    dist = (self.__moli_noh.GetConformer().GetAtomPosition(pair[0])
                          - self.__molj_noh.GetConformer().GetAtomPosition(pair[1])
                          - coord_delta).Length()
                    if dist > worstdist:
                        worstdist=dist
                        worstatomidx=atomidx
                    atomidx=atomidx+1

                if worstdist > max_deviation:
                    # Remove the furthest-away atom and try again
                    rwm = Chem.RWMol(self.mcs_mol)
                    print("REMOVING ATOM",worstatomidx," with distance", worstdist)
                    rwm.RemoveAtom(worstatomidx)
                    self.mcs_mol=Chem.Mol(rwm)
                else:
                    break




        def map_mcs_mol():
            """

            This function is used to define a map between the generated mcs, the
            molecules and vice versa
           
            """

            # mcs indexes mapped back to the first molecule moli


            # Get self-mapping
            mcsi_sub = tuple(range(self.mcs_mol.GetNumAtoms()))

            (moli_sub,molj_sub) = best_substruct_match(self.__moli_noh,self.__molj_noh,by_rmsd=self.options.threed)

            # mcs to moli
            map_mcs_mol_to_moli_sub = zip(mcsi_sub, moli_sub)

            # An RDkit atomic property is defined to store the mapping to moli
            for idx in map_mcs_mol_to_moli_sub:
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp('to_moli', str(idx[1]))

            mcsj_sub = tuple(range(self.mcs_mol.GetNumAtoms()))

            # mcs to molj
            map_mcs_mol_to_molj_sub = zip(mcsj_sub, molj_sub)

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
            # property. This could be very useful in the code development
            # when deletion or atom insertions are performed
            for at in self.mcs_mol.GetAtoms():
                at.SetProp('org_idx', str(at.GetIdx()))

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
                at.SetProp('rc', '0')

            rginfo = mol.GetRingInfo()

            rgs = rginfo.AtomRings()

            # print rgs

            rgs_set = set([e for l in rgs for e in l])

            for idx in rgs_set:
                for r in rgs:
                    if idx in r:
                        val = int(mol.GetAtomWithIdx(idx).GetProp('rc'))
                        val = val + 1
                        mol.GetAtomWithIdx(idx).SetProp('rc', str(val))
            return

        # Set logging level and format
        logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.INFO)

        self.options=options

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

        try:  # Try to sanitize the MCS molecule
            Chem.SanitizeMol(self.mcs_mol)
        except Exception:  # if not, try to recover the atom aromaticity wich is
            # important for the ring counter
            sanitFail = Chem.SanitizeMol(self.mcs_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
                                         catchErrors=True)
            if sanitFail:  # if not, the MCS is skipped
                raise ValueError('Sanitization Failed...')

        # Trim the MCS to remove atoms with too-large real-space deviations
        try:
            trim_mcs_mol(max_deviation=self.options.max3d)
        except Exception as e:
            raise ValueError(str(e))

        # Mapping between the found MCS molecule and moli,  molj
        try:
            map_mcs_mol()
        except Exception as e:
            raise ValueError(str(e))

        # Set the ring counters for each molecule
        set_ring_counter(self.__moli_noh)
        set_ring_counter(self.__molj_noh)
        set_ring_counter(self.mcs_mol)

        # for at in self.mcs_mol.GetAtoms():
        #     print 'at = %d rc = %d' % (at.GetIdx(), int(at.GetProp('rc')))

        if not options.verbose == 'pedantic':
            lg.setLevel(RDLogger.WARNING)

        return

    def get_map(self):
        """

        This function is used to return a list of pairs of atom indexes generated
        by the mapping between the two molecules used to calculate the MCS. 
        The calculated mapping is performed without considering hydrogens 

        Returns
        -------
        pair of indexes related to the atom mapping 

        """

        return self.__map_moli_molj

    # Note - not used in the main LOMAP calculation - here primarily for testing?
    @staticmethod
    def getMapping(moli, molj, hydrogens=False, fname=None, time_out=150):

        """
        Compute the MCS between two passed molecules
    
        Parameters
        ----------

        moli : RDKit molecule object 
            the first molecule used to perform the MCS calculation
        molj : RDKit molecule object 
            the second molecule used to perform the MCS calculation
        hydrogens : bool 
            incluse or not the hydrogens in the MCS calculation

        fname : string 
            the filename used to output a png file depicting the MCS mapping 

        time_out: int
            the max time in seconds used to compute the MCS

        Returns:
        --------
        map_moli_molj: python list of tuple [...(i,j)...]
            the list of tuple which contains the atom mapping indexes between 
            the two molecules. The indexes (i,j) are resplectively related to 
            the first (moli) and the second (molj) passed molecules 
                 
        """

        # Molecule copies
        moli_c = Chem.Mol(moli)
        molj_c = Chem.Mol(molj)

        if not hydrogens:
            moli_c = AllChem.RemoveHs(moli_c)
            molj_c = AllChem.RemoveHs(molj_c)

            # MCS calculaton. In RDKit the MCS is a smart string. Ring atoms are
        # always mapped in ring atoms. 
        mcs = rdFMCS.FindMCS([moli_c, molj_c],
                             timeout=time_out,
                             atomCompare=rdFMCS.AtomCompare.CompareAny,
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             matchValences=False,
                             ringMatchesRingOnly=True,
                             completeRingsOnly=False,
                             matchChiralTag=False)

        # Checking
        if mcs.canceled:
            raise ValueError('Timeout! No MCS found between passed molecules')

        if mcs.numAtoms == 0:
            raise ValueError('No MCS was found between the molecules')

        # The found MCS pattern (smart strings) is converted to a RDKit molecule
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

        try:
            Chem.SanitizeMol(mcs_mol)
        except Exception:  # if not try to recover the atom aromaticity wich is
            # important for the ring counter
            sanitFail = Chem.SanitizeMol(mcs_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
                                         catchErrors=True)
            if sanitFail:  # if not the MCS is skipped
                raise ValueError('Sanitization Failed...')

        # mcs indexes mapped back to the first molecule moli
        if moli_c.HasSubstructMatch(mcs_mol):
            moli_sub = moli_c.GetSubstructMatch(mcs_mol)
        else:
            raise ValueError('RDkit MCS Subgraph first molecule search failed')
        # mcs indexes mapped back to the second molecule molj
        if molj_c.HasSubstructMatch(mcs_mol):
            molj_sub = molj_c.GetSubstructMatch(mcs_mol)
        else:
            raise ValueError('RDkit MCS Subgraph second molecule search failed')

        mcs_sub = tuple(range(mcs_mol.GetNumAtoms()))

        # Map between the two molecules
        map_moli_to_molj = zip(moli_sub, molj_sub)

        # depict the mapping by using a .png file
        if fname:
            AllChem.Compute2DCoords(moli_c)
            AllChem.Compute2DCoords(molj_c)
            AllChem.Compute2DCoords(mcs_mol)
            DrawingOptions.includeAtomNumbers = True
            moli_fname = 'Moli'
            molj_fname = 'Molj'
            mcs_fname = 'Mcs'

            img = Draw.MolsToGridImage([moli_c, molj_c, mcs_mol],
                                       molsPerRow=3, subImgSize=(400, 400),
                                       legends=[moli_fname, molj_fname, mcs_fname],
                                       highlightAtomLists=[moli_sub, molj_sub, mcs_sub])

            img.save(fname)

            DrawingOptions.includeAtomNumbers = False

        return map_moli_to_molj

        ############ MCS BASED RULES ############

    # MCSR Rule
    # the mtansr method is not used but be retained here in case need to use in the future.
    def mtansr(self, ):
        """
        This rule computes the structural similarity between the two passed molecules 
        using the tanimoto score. 
        Returns
        -------
        scr_tan : float
            the rule score
        """
        fps_moli = FingerprintMols.FingerprintMol(self.moli)
        fps_molj = FingerprintMols.FingerprintMol(self.molj)
        scr_tan = DataStructs.FingerprintSimilarity(fps_moli, fps_molj)
        return scr_tan

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
        scr_mcsr = math.exp(-beta * (nha_moli + nha_molj - 2 * nha_mcs_mol))

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

        scr_mncar = float((nha_mcs_mol >= ths) or (nha_moli < ths + 3) or (nha_molj < ths + 3))

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
                mcs_conflict = [at.GetIdx() for at in mcs_mol.GetAtoms() if
                                int(at.GetProp('rc')) > 0 and not at.IsInRing()]

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

                for idx in range(0, len(fragments)):
                    lgt = len(fragments[idx])
                    if lgt > lgt_max:
                        lgt_max = lgt
                        max_idx = idx

                max_frag = fragments[max_idx]
                mcs_conflict = [at.GetIdx() for at in mcs_mol.GetAtoms() if not at.GetIdx() in max_frag]
                mcs_conflict.sort(reverse=True)
                edit_mcs_mol = Chem.EditableMol(mcs_mol)

                # WARNING: atom indexes have changed
                for i in mcs_conflict:
                    edit_mcs_mol.RemoveAtom(i)
                mcs_mol = edit_mcs_mol.GetMol()

                return mcs_mol

            mcs_conflict = set()
            for at in self.mcs_mol.GetAtoms():

                moli_idx = int(at.GetProp('to_moli'))
                molj_idx = int(at.GetProp('to_molj'))

                moli_idx_rc = int(self.__moli_noh.GetAtomWithIdx(moli_idx).GetProp('rc'))
                molj_idx_rc = int(self.__molj_noh.GetAtomWithIdx(molj_idx).GetProp('rc'))

                # Moli atom is a ring atom (rc>0) and its rc is different from 
                # the corresponding mcs rc atom  
                if moli_idx_rc > 0 and (moli_idx_rc != int(at.GetProp('rc'))):
                    if strict_flag:  # In strict mode we add the atom
                        mcs_conflict.add(at.GetIdx())
                    else:  # In loose mode we add the atom if it is not an aromatic atom
                        if not at.GetIsAromatic():
                            mcs_conflict.add(at.GetIdx())

                # Molj atom is a ring atom (rc>0) and its rc is different 
                # from the corresponding mcs rc atom 
                if molj_idx_rc > 0 and (molj_idx_rc != int(at.GetProp('rc'))):
                    if strict_flag:  # In strict mode we add the atom
                        mcs_conflict.add(at.GetIdx())
                    else:  # In loose mode we add the atom if it is not an aromatic atom
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
            return math.exp(-2 * beta * orig_nha_mcs_mol)

        ### Chiral Atoms ###

        # Deleting Chiral Atoms
        mcs_ring_set = set()
        mcs_chiral_set = set()

        for atom in mcs_mol_copy.GetAtoms():
            if atom.IsInRing():
                mcs_ring_set.add(atom.GetIdx())
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                mcs_chiral_set.add(atom.GetIdx())

        # Loop over the mcs chiral atoms to check if they are also ring atoms
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
            return math.exp(-2 * beta * orig_nha_mcs_mol)

        # self.draw_molecule(mcs_mol_copy)
        fragments = Chem.rdmolops.GetMolFrags(mcs_mol_copy)
        max_idx = 0
        lgt_max = 0

        for idx in range(0, len(fragments)):
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

        return math.exp(-2 * beta * (orig_nha_mcs_mol - max_frag_num_heavy_atoms))

    # AtomicNumber rule 
    def atomic_number_rule(self,beta=0.1):

        """
        This rule checks how many elements have been changed in the MCS 
        and a score based on the fraction of MCS matches that are the same atomicnumber.
        When used with beta=0.1 and multiplied by mcsr, this is equivalent to counting
        mismatched atoms at only half weight.
             
        """
        natoms=0
        nmatch=0
        for at in self.mcs_mol.GetAtoms():
            natoms=natoms+1
            moli_idx = int(at.GetProp('to_moli'))
            molj_idx = int(at.GetProp('to_molj'))
            moli_a = self.moli.GetAtoms()[moli_idx]
            molj_a = self.molj.GetAtoms()[molj_idx]

            if moli_a.GetAtomicNum() == molj_a.GetAtomicNum():
                nmatch=nmatch+1

        return math.exp(-1 * beta * (natoms-nmatch))

    # Sulfonamides rule
    def sulfonamides_rule(self):

        """
        This rule checks to see if we are growing a complete sulfonamide, and 
        returns 0 if we are. This means that if this rule is used we effectively disallow
        this transition. Testing has shown that growing -SO2NH2 from scratch performs
        very badly.
             
        """

        def adds_sulfonamide(mol):
            """
            Returns true if the removal of the MCS from the provided molecule
            leaves a sulfonamide
            """

            if not mol.HasSubstructMatch(self.mcs_mol):
                raise ValueError('RDkit MCS Subgraph molecule search failed in sulfonamide check')

            
            rwm=rdmolops.DeleteSubstructs(mol, self.mcs_mol)
            return rwm.HasSubstructMatch(Chem.MolFromSmarts('S(=O)(=O)N'))

        retval = 0 if (adds_sulfonamide(self.__moli_noh)) else 1
        retval = 0 if (adds_sulfonamide(self.__molj_noh)) else retval
        return retval

    # Heterocycles rule
    def heterocycles_rule(self):

        """
        This rule checks to see if we are growing a heterocycle from a hydrogen, and 
        returns 0 if we are. This means that if this rule is used we effectively disallow
        this transition. Testing has shown that growing a pyridine or other heterocycle
        is unlikely to work (better to grow phenyl then mutate)
             
        """

        def adds_heterocycle(mol):
            """
            Returns true if the removal of the MCS from the provided molecule
            leaves a sulfonamide
            """

            if not mol.HasSubstructMatch(self.mcs_mol):
                raise ValueError('RDkit MCS Subgraph molecule search failed in sulfonamide check')

            
            rwm=rdmolops.DeleteSubstructs(mol, self.mcs_mol)
            # Only picking up N/C containing heterocycles - odd cases like pyran derivatives are not caught
            grow6mheterocycle =  rwm.HasSubstructMatch(Chem.MolFromSmarts('[n]1[c,n][c,n][c,n][c,n][c,n]1'))

            # Note that growing pyrrole, furan or thiophene is allowed
            grow5mheterocycle =  rwm.HasSubstructMatch(Chem.MolFromSmarts('[o,n&X3,s]1[n][c,n][c,n][c,n]1'))
            grow5mheterocycle |=  rwm.HasSubstructMatch(Chem.MolFromSmarts('[o,n&X3,s]1[c,n][n][c,n][c,n]1'))
            return (grow6mheterocycle | grow5mheterocycle)



        retval = 0 if (adds_heterocycle(self.__moli_noh)) else 1
        retval = 0 if (adds_heterocycle(self.__molj_noh)) else retval
        return retval

    def transmuting_methyl_into_ring_rule(self):

        """
         Rule to prevent turning a methyl into a ring atom and similar transformations
         (you can grow a ring, but you can't transmute into one)

        """
        moli=self.__moli_noh
        molj=self.__molj_noh

        # Get list of bonds in mol i and j that go from the MCS to a non-MCS atom,
        # arranged in tuples with the index of the MCS atom
        moli_sub = moli.GetSubstructMatch(self.mcs_mol)
        molj_sub = molj.GetSubstructMatch(self.mcs_mol)

        is_bad=False

        for i in range(0,len(moli_sub)):
            edge_bondsi = [ b.GetBeginAtomIdx() for b in moli.GetBonds() if (b.GetEndAtomIdx()==moli_sub[i] and not b.GetBeginAtomIdx() in moli_sub) ]
            edge_bondsi += [ b.GetEndAtomIdx() for b in moli.GetBonds() if (b.GetBeginAtomIdx()==moli_sub[i] and not b.GetEndAtomIdx() in moli_sub) ]
            edge_bondsj = [ b.GetBeginAtomIdx() for b in molj.GetBonds() if (b.GetEndAtomIdx()==molj_sub[i] and not b.GetBeginAtomIdx() in molj_sub) ]
            edge_bondsj += [ b.GetEndAtomIdx() for b in molj.GetBonds() if (b.GetBeginAtomIdx()==molj_sub[i] and not b.GetEndAtomIdx() in molj_sub) ]
            #print("Atom",i,"index",moli_sub[i],"edge atoms on mol 1 are",edge_bondsi);
            #print("Atom",i,"index",molj_sub[i],"edge atoms on mol 2 are",edge_bondsj);

            for edgeAtom_i in edge_bondsi:
                for edgeAtom_j in edge_bondsj:
                    if (moli.GetAtomWithIdx(edgeAtom_i).IsInRing() ^ molj.GetAtomWithIdx(edgeAtom_j).IsInRing()):
                        is_bad=True

        return 0 if is_bad else 1


if "__main__" == __name__:

    #mola = Chem.MolFromMol2File('../test/basic/2-methylnaphthalene.mol2', sanitize=False, removeHs=False)
    #molb = Chem.MolFromMol2File('../test/basic/2-naftanol.mol2', sanitize=False, removeHs=False)
    #mola = Chem.MolFromMolFile('../test/transforms/chlorotoluyl1.sdf', sanitize=False, removeHs=False)
    #molb = Chem.MolFromMolFile('../test/transforms/chlorotoluyl2.sdf', sanitize=False, removeHs=False)
    #mola = Chem.MolFromMolFile('../test/transforms/toluyl3.sdf', sanitize=False, removeHs=False)
    mola = Chem.MolFromMolFile('../test/transforms/chlorophenol.sdf', sanitize=False, removeHs=False)
    molb = Chem.MolFromMolFile('../test/transforms/phenylfuran.sdf', sanitize=False, removeHs=False)

    mp = MCS.getMapping(mola, molb, hydrogens=False, fname='mcs.png')

    print(mp)

    # MCS calculation
    try:
        MC = MCS(mola, molb)
    except Exception:
        raise ValueError('NO MCS FOUND......')

    # # Rules calculations
    mcsr = MC.mcsr()
    mncar = MC.mncar()
    atnum = MC.atomic_number_rule()

    strict = MC.tmcsr(strict_flag=True)
    loose = MC.tmcsr(strict_flag=False)

    print('TMCRS STRICT = %f , TMCRS LOOSE = %f' % (strict, loose))
    print('MCSR = ', mcsr)
    print('MNCAR = ', mncar)
    print('ATNUM = ', atnum)

    tmp = mcsr * mncar

    print('Total Strict = %f , Total Loose = %f' % (tmp * strict, tmp * loose))

    for at in MC.mcs_mol.GetAtoms():
        moli_idx = int(at.GetProp('to_moli'))
        molj_idx = int(at.GetProp('to_molj'))
        moli_a = mola.GetAtoms()[moli_idx]
        molj_a = molb.GetAtoms()[molj_idx]
        print("MCS match: ",moli_idx,moli_a.GetAtomicNum(),molj_idx,molj_a.GetAtomicNum())

    print("sulfonamides:",MC.sulfonamides_rule())
    print("heterocycles:",MC.heterocycles_rule())
    print("growring:",MC.transmuting_methyl_into_ring_rule())

