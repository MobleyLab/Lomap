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


@pytest.fixture(scope='session')
def toluene():
    return Chem.MolFromMol2File(_rf('basic/toluene.mol2'))


@pytest.fixture(scope='session')
def methylcyclohexane():
    return Chem.MolFromMol2File(_rf('basic/methylcyclohexane.mol2'))


@pytest.fixture(scope='session')
def methylnaphthalene():
    return Chem.MolFromMol2File(_rf('basic/2-methylnaphthalene.mol2'))


@pytest.fixture(scope='session')
def trimethylnaphthalene():
    return Chem.MolFromMol2File(_rf('basic/1,3,7-trimethylnaphthalene.mol2'))


@pytest.fixture(scope='session')
def butylmethylbenzene():
    return Chem.MolFromMol2File(_rf('basic/1-butyl-4-methylbenzene.mol2'))


@pytest.fixture(scope='session')
def dimethylnaphthalene():
    return Chem.MolFromMol2File(_rf('basic/2,6-dimethylnaphthalene.mol2'))


@pytest.fixture(scope='session')
def methylpropylnaphthalene():
    return Chem.MolFromMol2File(_rf('basic/2-methyl-6-propylnaphthalene.mol2'))


@pytest.fixture(scope='session')
def naphthol():
    return Chem.MolFromMol2File(_rf('basic/2-naftanol.mol2'))


@pytest.fixture(scope='function', params=[
    'methylcyclohexane', 'methylnaphthalene', 'trimethylnaphthalene',
    'butylmethylbenzene', 'dimethylnaphthalene', 'methylpropylnaphthalene',
    'naphthol',
])
def regression_mappings(toluene, methylcyclohexane, methylnaphthalene,
                        trimethylnaphthalene, butylmethylbenzene,
                        dimethylnaphthalene, methylpropylnaphthalene,
                        naphthol, request):
    regression_mappings = {
        'methylcyclohexane': (methylcyclohexane,
                              [(0, 0), (1, 1), (2, 6), (3, 5), (4, 4), (5, 3),
                               (6, 2)]),
        'methylnaphthalene':  (methylnaphthalene,
                               [(0, 0), (1, 1), (2, 2), (3, 3), (4, 8), (5, 9),
                               (6, 10)]),
        'trimethylnaphthalene': (trimethylnaphthalene,
                                 [(0, 0), (1, 1), (2, 10), (3, 9), (4, 8),
                                  (5, 3), (6, 2)]),
        'butylmethylbenzene': (butylmethylbenzene,
                               [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8),
                                (6, 9)]),
        'dimethylnaphthalene': (dimethylnaphthalene,
                                [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4),
                                 (5, 10), (6, 11)]),
        'methylpropylnaphthalene': (methylpropylnaphthalene,
                                    [(0, 2), (1, 3), (2, 8), (3, 7), (4, 6),
                                     (5, 5), (6, 4)]),
        'naphthol': (naphthol,
                     [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 9), (6, 10)])
    }
    other, other_mapping = regression_mappings[request.param]
    # a bit paranoid, but return copies
    return Chem.Mol(toluene), Chem.Mol(other), other_mapping


def test_toluene_to_other(regression_mappings):
    toluene, other, ref_mapping = regression_mappings

    mapper = mcs.MCS(toluene, other)

    mapping = mapper.heavy_atom_mcs_map()

    assert sorted(mapping) == ref_mapping
