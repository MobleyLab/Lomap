"""
Regression tests for the MCS class
"""
import pkg_resources
import pytest
from rdkit import Chem

from lomap import mcs

def _rf(fn):
    # get path to file from inside lomap installation
    f = pkg_resources.resource_filename('lomap', 'test/' + fn)
    return f.replace('/lomap/test', '/test')


@pytest.fixture
def toluene():
    f = _rf('basic/toluene.mol2')

    return Chem.MolFromMol2File(f)


@pytest.fixture
def methylcyclohexane():
    f = _rf('basic/methylcyclohexane.mol2')

    return Chem.MolFromMol2File(f)


def test_toluene_to_methylcyclohexane(toluene, methylcyclohexane):
    mapper = mcs.MCS(toluene, methylcyclohexane)

    mapping = mapper.heavy_atom_mcs_map()

    assert sorted(mapping) == [(0, 0), (1, 1), (2, 6), (3, 5), (4, 4),
                               (5, 3), (6, 2)]
