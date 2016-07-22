r"""
Various rules to compute the similarity score.
"""

import copy_reg
import types
from functools import partial

import numpy as np
import rdkit.Chem as rdchem
from rdkit.Chem.rdFMCS import FindMCS, AtomCompare, BondCompare

import mol
from logger import logger


def check_mcs(method):
    """
    Decorator to ensure that the MCS has been computed.
    """

    def wrapper(self, *args, **kwargs):
        if not self.mapping:
            self._compute_mcs()

        result = method(self)
        return result

    return wrapper

# FMCS parameters, mostly defaults
_fmcs_params = dict(maximizeBonds = False, threshold = 1.0,
                    verbose = False, matchValences = False,
                    ringMatchesRingOnly = True, completeRingsOnly = True,
                    atomCompare = AtomCompare.CompareIsotopes,
                    bondCompare = BondCompare.CompareAny)

class Rules(object):
    def __init__(self, min_mcs=4, maxtime=60):
        self.min_mcs = min_mcs
        self.maxtime = maxtime
        self.mapping = None

    def update_pair(self, morph_pair):
        self.morph_pair = morph_pair
        self.mol0 = morph_pair.mol0.molecule
        self.mol1 = morph_pair.mol1.molecule

    def same_charges(self):
        """
        Check if total molecule charge is the same.

        :param mol0: molecule 0
        :type mol0: mol.Molecule
        :param mol1:  molecule 1
        :type mol1: mol.Molecule
        """

        charge0 = 0.0

        for atom in self.mol0.GetAtoms():
            charge0 += atom.GetFormalCharge()

        charge1 = 0.0

        for atom in self.mol1.GetAtoms():
            charge1 += atom.GetFormalCharge()

        if abs(charge0 - charge1) < 1.0e-3:
            score = 1
        else:
            score = 0

        return score

    def _compute_mcs(self):
        """
        Compute the MCS for the morph pair using FMCS from RDKit
        """

        _fmcs_params.update(timeout=int(self.maxtime))

        # FIXME: heavyu atoms only
        mcs = FindMCS( (self.mol0, self.mol1), **_fmcs_params)

        smarts = mcs.smartsString
        completed = not mcs.canceled

        if not completed:
            logger.warn('MCSS timed out after %.2fs' % self.maxtime)

        if not smarts:
            # FIXME: more detailled error message
            raise ValueError('No MCS was found between the morph pair')

        pattern = rdchem.MolFromSmarts(smarts)

        # FIXME: how to deal with multiple MCSs
        match0 = self.mol0.GetSubstructMatch(pattern)
        match1 = self.mol1.GetSubstructMatch(pattern)

        self.mapping = zip(match0, match1)
        self.morph_pair.mapping = self.mapping

    @check_mcs
    def minimum_mcs(self):
        """
        Check for minimum size i.e. the number of atoms in the MCS
        """

        # FIXME: heavyu atoms only
        natoms_mcs = len(self.mapping)
        natoms_mol0 = self.mol0.GetNumAtoms()
        natoms_mol1 = self.mol1.GetNumAtoms()
        threshold = self.min_mcs + 3

        return float(natoms_mcs >= self.min_mcs or natoms_mol0 < threshold or
                     natoms_mol0 < threshold)


class Scorer(object):
    def __init__(self, rules):
        self.rules_class = rules

    def compute(self, mol0, mol1):
        all_scores = []
        self.rules = self.rules_class(min_mcs=4, maxtime=20)
        self.scorers = [self.rules.same_charges, self.rules.minimum_mcs]
        morph_pair = mol.MorphPair(mol0, mol1)
        self.rules.update_pair(morph_pair)

        # NOTE: In parallel mode no two processes must access the same
        #       morph pair!

        for scorer in self.scorers:
            part_score = scorer()

            # the total score is the product of all individual scores
            if part_score == 0:
                all_scores = []
                break

            all_scores.append(part_score)

        if all_scores:
            total_score = reduce(lambda x, y: x * y, all_scores)
        else:
            total_score = 0.0

        # FIXME: how to differentiate between strict and loose rule
        morph_pair.strict_score = total_score
        morph_pair.loose_score = total_score

        return morph_pair


# pickling support for multiprocessing
# https://stackoverflow.com/questions/25156768/cant-pickle-type-instancemethod-using-pythons-multiprocessing-pool-apply-a
def _pickle_method(m):
    if m.__self__ is None:
        return getattr, (m.__class__, m.im_func.func_name)
    else:
        return getattr, (m.__self__, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)

def compute_similarity_matrix(mols, maxtime, nproc):
    """
    """

    scorer = Scorer(Rules)

    scores = []
    N = len(mols)
    simmat = np.zeros(shape=(N,N), dtype=mol.MorphPair)

    if nproc > 1 or nproc <= 0:
        import multiprocessing as mp

        parallel = True
        maxproc = mp.cpu_count()

        if nproc > maxproc:
            logger.warn('limiting number of processors to %i' % maxproc)
            nproc = maxproc
        elif nproc <= 0:
            nproc = maxproc

        pool = mp.Pool(nproc)
        map_func = pool.imap
    else:
        parallel = False
        map_func = map

    for i in range(N-1):
        partial_func = partial(scorer.compute, mols[i])
        scores.append(map_func(partial_func, mols[i+1:N]) )

    for i, row in enumerate(scores):
        simmat[i][i+1:N] = [s for s in row]

    if parallel:
        pool.close()
        pool.join()

    return simmat
