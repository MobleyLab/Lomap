import pytest
from lomap.mcs import MCS
from rdkit import RDLogger
import pickle
from lomap.dbmol import DBMolecules
import multiprocessing


def test_mcs():
    f = open('test/basic/MCS.pickle', 'rb')
    data = pickle.load(f)
    data_no_hydrogens = data[0]
    data_hydrogens = data[1]

    db = DBMolecules('test/basic/', parallel=1, verbose='off', output=False, time=20, ecrscore=0.0 ,name='out', display=False, max=6, cutoff=0.4, radial=False, hub=None)


    nohyds = {}
    hyds = {}

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    for i in range(0, db.nums()):
        for j in range(i + 1, db.nums()):
            MCS_no_hyds = MCS.getMapping(db[i].getMolecule(), db[j].getMolecule())
            MCS_hyds = MCS.getMapping(db[i].getMolecule(), db[j].getMolecule(), hydrogens=True)

            nohyds[(i, j)] = list(MCS_no_hyds)
            hyds[(i, j)] = list(MCS_hyds)

    assert(nohyds == data_no_hydrogens)
    assert(hyds == data_hydrogens)

# Check serial and parallel mode
def test_serial_parallel():
    db =  DBMolecules('test/basic')
    s_strict, s_loose = db.build_matrices()
    db.options.paralell = multiprocessing.cpu_count()
    p_strict, p_loose = db.build_matrices()
    
    assert (all(s_strict == p_strict))
    assert (all(s_loose == p_loose))

def test_read_mol2_files():
    db =  DBMolecules('test/basic')
    db.options.directory = 'test/'
    with pytest.raises(IOError):
        db.read_mol2_files()
