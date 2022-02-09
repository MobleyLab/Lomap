"""
Regression tests for the MCS class

These might start failing if RDKit changes canonical order of atoms.
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


@pytest.fixture
def methylnaphthalene():
    f = _rf('basic/2-methylnaphthalene.mol2')

    return Chem.MolFromMol2File(f)


@pytest.fixture
def trimethylnaphthalene():
    f = _rf('basic/1,3,7-trimethylnaphthalene.mol2')

    return Chem.MolFromMol2File(f)


@pytest.fixture
def butylmethylbenzene():
    f = _rf('basic/1-butyl-4-methylbenzene.mol2')

    return Chem.MolFromMol2File(f)


@pytest.fixture
def dimethylnaphthalene():
    f = _rf('basic/2,6-dimethylnaphthalene.mol2')

    return Chem.MolFromMol2File(f)


@pytest.fixture
def methylpropylnaphthalene():
    f = _rf('basic/2-methyl-6-propylnaphthalene.mol2')

    return Chem.MolFromMol2File(f)


@pytest.fixture
def naphthol():
    f = _rf('basic/2-naftanol.mol2')

    return Chem.MolFromMol2File(f)


def test_toluene_to_methylcyclohexane(toluene, methylcyclohexane):
    mapper = mcs.MCS(toluene, methylcyclohexane)

    mapping = mapper.heavy_atom_mcs_map()

    assert sorted(mapping) == [(0, 0), (1, 1), (2, 6), (3, 5), (4, 4), (5, 3),
                               (6, 2)]


def test_toluene_to_methylnaphthalene(toluene, methylnaphthalene):
    mapper = mcs.MCS(toluene, methylnaphthalene)

    mapping = mapper.heavy_atom_mcs_map()

    assert len(mapping) == 7
    assert sorted(mapping) == [(0, 0), (1, 1), (2, 2), (3, 3), (4, 8), (5, 9),
                               (6, 10)]


def test_toluene_to_trimethylnaphthalene(toluene, trimethylnaphthalene):
    mapper = mcs.MCS(toluene, trimethylnaphthalene)

    mapping = mapper.heavy_atom_mcs_map()

    assert sorted(mapping) == [(0, 0), (1, 1), (2, 10), (3, 9), (4, 8), (5, 3),
                               (6, 2)]


def test_toluene_to_naphthol(toluene, naphthol):
    mapper = mcs.MCS(toluene, naphthol)

    mapping = mapper.heavy_atom_mcs_map()

    assert sorted(mapping) == [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 9),
                               (6, 10)]


def test_toluene_to_butylmethylbenzene(toluene, butylmethylbenzene):
    mapper = mcs.MCS(toluene, butylmethylbenzene)

    mapping = mapper.heavy_atom_mcs_map()

    assert sorted(mapping) == [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8),
                               (6, 9)]


def test_toluene_dimethylnaphthalene(toluene, dimethylnaphthalene):
    mapper = mcs.MCS(toluene, dimethylnaphthalene)

    mapping = mapper.heavy_atom_mcs_map()

    assert sorted(mapping) == [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 10),
                               (6, 11)]


def test_toluene_to_methylpropylnaphthalene(toluene, methylpropylnaphthalene):
    mapper = mcs.MCS(toluene, methylpropylnaphthalene)

    mapping = mapper.heavy_atom_mcs_map()

    assert sorted(mapping) == [(0, 2), (1, 3), (2, 8), (3, 7), (4, 6), (5, 5),
                               (6, 4)]