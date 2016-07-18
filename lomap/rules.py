r"""
Various rules to compute the similarity score.
"""

import copy_reg
import types
from functools import partial

import numpy as np

import mol
from logger import logger


# https://stackoverflow.com/questions/25156768/cant-pickle-type-instancemethod-using-pythons-multiprocessing-pool-apply-a
def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)


def same_charges(mol0, mol1, *args):
    """
    Check if total molecule charge is the same.

    :param mol0: molecule 0
    :type mol0: mol.Molecule
    :param mol1:  molecule 1
    :type mol1: mol.Molecule
    """

    charge0 = 0.0

    for atom in mol0.molecule.GetAtoms():
        charge0 += atom.GetFormalCharge()

    charge1 = 0.0

    for atom in mol1.molecule.GetAtoms():
        charge1 += atom.GetFormalCharge()

    if abs(charge0 - charge1) < 1.0e-3:
        score = 1
    else:
        score = 0

    return score


class Scorer(object):
    def __init__(self, scorers = []):
        self.scorers = scorers

    def compute(self, *args):
        all_scores = []

        for scorer in self.scorers:
            part_score = scorer(*args)

            # the total score is the product of all individual scores
            if part_score == 0:
                return 0.0

            all_scores.append(part_score)

        return reduce(lambda x, y: x * y, all_scores)


def compute_similarity_matrix(mols, methods, nproc):
    """
    """

    scores = []
    scorer = Scorer(methods)
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
