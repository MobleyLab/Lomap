import pytest
from lomap.mcs import MCS
from rdkit import RDLogger
import pickle
from lomap.dbmol import DBMolecules
import multiprocessing
import networkx as nx
import networkx.algorithms.isomorphism as iso
import subprocess


@pytest.fixture
def executable():
    return 'lomap'


def test_insufficient_arguments(executable):
    cmd = [executable]
    error_string = b'error: the following arguments are required: directory'
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    assert(error_string in stderr)


def test_mcs():
    f = open('test/basic/MCS.pickle', 'rb')
    data = pickle.load(f)
    data_no_hydrogens = data[0]
    data_hydrogens = data[1]

    db = DBMolecules('test/basic/', parallel=1, verbose='off', output=False, time=20, ecrscore=0.0, name='out',
                     display=False, max=6, cutoff=0.4, radial=False, hub=None)

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


# Check iter and next
def test_iter_next():
    inst = DBMolecules('test/basic/', parallel=1, verbose='off', output=False, time=20, ecrscore=0.0, name='out',
                       display=False, max=6, cutoff=0.4, radial=False, hub=None)
    for i in range(0, inst.nums()):
        inst.next()
    with pytest.raises(StopIteration):
        inst.next()


# Check get, set and add
def test_get_set_add():
    inst = DBMolecules('test/basic/', parallel=1, verbose='off', output=False, time=20, ecrscore=0.0, name='out',
                       display=False, max=6, cutoff=0.4, radial=False, hub=None)
    with pytest.raises(IndexError):
        inst.__getitem__(inst.nums()+1)
    with pytest.raises(IndexError):
        inst.__setitem__(inst.nums()+1,inst[1])
    with pytest.raises(ValueError):
        inst.__setitem__(0, 'no_mol_obj')
    with pytest.raises(ValueError):
        inst.__add__('no_mol_obj')


# Check serial and parallel mode
def test_serial_parallel():
    db = DBMolecules('test/basic')
    s_strict, s_loose = db.build_matrices()
    db.options.paralell = multiprocessing.cpu_count()
    p_strict, p_loose = db.build_matrices()
    
    assert (all(s_strict == p_strict))
    assert (all(s_loose == p_loose))


def test_read_mol2_files():
    db = DBMolecules('test/basic')
    db.options.directory = 'test/'
    with pytest.raises(IOError):
        db.read_mol2_files()


# Testing the graphs
# def test_graph():
#     db = DBMolecules('test/basic/', parallel=1, verbose='off', output=False, time=20, ecrscore=0.0, name='out',
#                      display=False, max=6, cutoff=0.4, radial=False, hub=None)
#
#     strict, loose = db.build_matrices()
#     graph = db.build_graph()
#
#     mol2_graph = nx.read_gpickle("test/basic/molecules.gpickle", 'rb')
#
#     dic1_nodes = graph.nodes()
#     dic1_edges = graph.edges()
#
#     dic2_nodes = mol2_graph.nodes()
#     dic2_edges = mol2_graph.edges()
#
#     assert(dic1_nodes.keys() == dic2_nodes.keys())
#     assert(dic1_edges.keys() == dic2_edges.keys())
#
#     nm = iso.categorical_node_match(['fname_comp', 'ID'], ['noname', -1])
#     em = iso.categorical_edge_match(['strict_flag', 'similarity'], [False, -1.0])
#
#     assert (nx.is_isomorphic(graph, mol2_graph, node_match=nm, edge_match=em))
#
#
# def test_graph_radial():
#     db = DBMolecules('test/radial/', radial=True)
#
#     strict, loose = db.build_matrices()
#     graph = db.build_graph()
#
#     mol2_graph = nx.read_gpickle("test/radial/radial.gpickle")
#
#     dic1_nodes = graph.nodes()
#     dic1_edges = graph.edges()
#
#     dic2_nodes = mol2_graph.nodes()
#     dic2_edges = mol2_graph.edges()
#
#     assert(dic1_nodes.keys() == dic2_nodes.keys())
#     assert(dic1_edges.keys() == dic2_edges.keys())
#
#     nm = iso.categorical_node_match(['fname_comp', 'ID'], ['noname', -1])
#     em = iso.categorical_edge_match(['strict_flag', 'similarity'], [False, -1.0])
#
#     assert(nx.is_isomorphic(graph, mol2_graph, node_match=nm, edge_match=em))
#
#
# def test_graph_radial_hub():
#
#     db = DBMolecules('test/radial/', radial = True, hub="ejm_46.mol2")
#
#     strict,loose = db.build_matrices()
#     graph = db.build_graph()
#
#     mol2_graph = nx.read_gpickle("test/radial/radial_hub.gpickle")
#
#     dic1_nodes = graph.nodes()
#     dic1_edges = graph.edges()
#
#     dic2_nodes = mol2_graph.nodes()
#     dic2_edges = mol2_graph.edges()
#
#     assert(dic1_nodes.keys() == dic2_nodes.keys())
#     assert(dic1_edges.keys() == dic2_edges.keys())
#
#     nm = iso.categorical_node_match(['fname_comp','ID'],['noname',-1])
#     em = iso.categorical_edge_match(['strict_flag','similarity'],[False,-1.0])
#
#     assert(nx.is_isomorphic(graph, mol2_graph, node_match=nm, edge_match=em))
#
#
# def test_graph_radial_hub_fingerprint():
#
#     db = DBMolecules('test/radial/', radial=True, fingerprint=True, hub="ejm_46.mol2")
#
#     strict,loose = db.build_matrices()
#     graph = db.build_graph()
#
#     mol2_graph = nx.read_gpickle("test/radial/radial_hub_fingerprint.gpickle")
#
#     dic1_nodes = graph.nodes()
#     dic1_edges = graph.edges()
#
#     dic2_nodes = mol2_graph.nodes()
#     dic2_edges = mol2_graph.edges()
#
#     assert(dic1_nodes.keys() == dic2_nodes.keys())
#     assert(dic1_edges.keys() == dic2_edges.keys())
#
#     nm = iso.categorical_node_match(['fname_comp','ID'],['noname',-1])
#     em = iso.categorical_edge_match(['strict_flag','similarity'],[False,-1.0])
#
#     assert(nx.is_isomorphic(graph, mol2_graph , node_match=nm, edge_match=em))
#
#
# def test_graph_radial_hub_fingerprint():
#
#     db = DBMolecules('test/radial/', radial=True, fingerprint=True, fast=True, hub="ejm_46.mol2")
#
#     strict,loose = db.build_matrices()
#     graph = db.build_graph()
#
#     mol2_graph = nx.read_gpickle("test/radial/radial_hub_fingerprint_fast.gpickle")
#
#     dic1_nodes = graph.nodes()
#     dic1_edges = graph.edges()
#
#     dic2_nodes = mol2_graph.nodes()
#     dic2_edges = mol2_graph.edges()
#
#     assert(dic1_nodes.keys() == dic2_nodes.keys())
#     assert(dic1_edges.keys() == dic2_edges.keys())
#
#     nm = iso.categorical_node_match(['fname_comp','ID'],['noname',-1])
#     em = iso.categorical_edge_match(['strict_flag','similarity'],[False,-1.0])
#
#     assert(nx.is_isomorphic(graph, mol2_graph , node_match=nm, edge_match=em))
