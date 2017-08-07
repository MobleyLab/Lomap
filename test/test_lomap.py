import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

import unittest
from unittest import skipIf
from lomap.dbmol import DBMolecules
from lomap.graphgen import GraphGen
from lomap.mcs import MCS
import argparse
import multiprocessing
import networkx as nx
import networkx.algorithms.isomorphism as iso
import pickle
from rdkit import RDLogger

# Python graph section must be update to fix a bug
py_ver = int(sys.version[0])

GR_COMP = False

if py_ver < 3:
    GR_COMP=True

    
class TestLomap(unittest.TestCase):
    
    def setUp(self):
        self.inst = DBMolecules('test/basic/', parallel=1, verbose='off', output=False, time=20, ecrscore=0.0 ,name='out', display=False, max=6, cutoff=0.4, radial=False, hub=None)
    
    # Test class Inizialitazion

    # Check wrong number of argumetns or wrong options
    def test_InsufficientArgs(self):
        self.assertRaises(TypeError, self.inst.__init__)
       
    # Check passed input Types
    def test_TypeArgs(self):
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/nodir/')
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/basic', time=-1)
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/basic', parallel=-1)
        self.assertRaises(SystemExit, self.inst.__init__, 'test/basic', parallel=-1.5)
        self.assertRaises(SystemExit, self.inst.__init__, 'test/basic', verbose='err_option')
        self.assertRaises(SystemExit, self.inst.__init__, 'test/basic', time=-1.5)
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/basic', ecrscore=-1.5)
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/basic', ecrscore=2.0)
        self.assertRaises(TypeError, self.inst.__init__, 'test/basic', output=-5.0)
        self.assertRaises(TypeError, self.inst.__init__, 'test/basic', radial=-5.0)
        self.assertRaises(TypeError, self.inst.__init__, 'test/basic', display=-5.0)
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/basic', max=-5)
        self.assertRaises(SystemExit, self.inst.__init__, 'test/basic', max=-5.0)
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/basic', cutoff=-5)
        self.assertRaises(SystemExit, self.inst.__init__, 'test/basic', max='string')

    # Test class methods
    
    # Check iter and next
    def test_iter_next(self):
        for i in range(0,self.inst.nums()):
            self.inst.next()        
        self.assertRaises(StopIteration, self.inst.next)

    # Check get, set and add
    def test_get_set_add(self):
        self.assertRaises(IndexError, self.inst.__getitem__, self.inst.nums()+1)
        self.assertRaises(IndexError, self.inst.__setitem__, self.inst.nums()+1, self.inst[1])
        self.assertRaises(ValueError, self.inst.__setitem__, 0, 'no_mol_obj')
        self.assertRaises(ValueError, self.inst.__add__, 'no_mol_obj')

    # Check read mol2 files
    def test_read_mol2_files(self):        
        db =  DBMolecules('test/basic')
        db.options.directory = 'test/'
        self.assertRaises(IOError, db.read_mol2_files)
        

    # Check serial and parallel mode
    def test_serial_parallel(self):
        db =  DBMolecules('test/basic')
        s_strict, s_loose = db.build_matrices()
        db.options.paralell = multiprocessing.cpu_count()
        p_strict, p_loose = db.build_matrices()
        
        self.assertEqual(True, all(s_strict == p_strict))
        self.assertEqual(True, all(s_loose == p_loose))
    
    # Check Graph
    @skipIf(not GR_COMP, 'The graph test has been skipped untill a bug in the graph generation between py2 and py3 will be fixed')
    def test_graph(self):

        db =  self.inst
        
        strict,loose = db.build_matrices()
        graph = db.build_graph()
        
        self.assertRaises(IOError, nx.nx_agraph.write_dot, graph, '/check.dot') 

        mol2_graph = nx.read_gpickle("test/basic/molecules.gpickle")

        dic1_nodes = graph.nodes(data=True)
        dic1_edges = graph.edges(data=True)

        dic2_nodes = mol2_graph.nodes(data=True)
        dic2_edges = mol2_graph.edges(data=True)

        self.assertEqual(True, dic1_nodes == dic2_nodes)
        self.assertEqual(True, dic2_edges == dic2_edges)
        
        nm = iso.categorical_node_match(['fname_comp','ID'],['noname',-1])
        em = iso.categorical_edge_match(['strict_flag','similarity'],[False,-1.0])
        
        self.assertEqual(True, nx.is_isomorphic(graph, mol2_graph , node_match=nm, edge_match=em))

    #check the radial graph option
    @skipIf(not GR_COMP, 'The graph test has been skipped untill a bug in the graph generation between py2 and py3 will be fixed')
    def test_graph_radial(self):

        db = DBMolecules('test/radial/', radial = True)
        
        strict,loose = db.build_matrices()
        graph = db.build_graph()
        
        mol2_graph = nx.read_gpickle("test/radial/radial.gpickle")

        dic1_nodes = graph.nodes(data=True)
        dic1_edges = graph.edges(data=True)

        dic2_nodes = mol2_graph.nodes(data=True)
        dic2_edges = mol2_graph.edges(data=True)

        self.assertEqual(True, dic1_nodes == dic2_nodes)
        self.assertEqual(True, dic2_edges == dic2_edges)
        
        nm = iso.categorical_node_match(['fname_comp','ID'],['noname',-1])
        em = iso.categorical_edge_match(['strict_flag','similarity'],[False,-1.0])

        self.assertEqual(True, nx.is_isomorphic(graph, mol2_graph , node_match=nm, edge_match=em))

    #check the radial hub graph option 
    @skipIf(not GR_COMP, 'The graph test has been skipped untill a bug in the graph generation between py2 and py3 will be fixed')
    def test_graph_radial(self):

        db = DBMolecules('test/radial/', radial = True, hub="ejm_46.mol2")
        
        strict,loose = db.build_matrices()
        graph = db.build_graph()
        
        mol2_graph = nx.read_gpickle("test/radial/radial_hub.gpickle")

        dic1_nodes = graph.nodes(data=True)
        dic1_edges = graph.edges(data=True)

        dic2_nodes = mol2_graph.nodes(data=True)
        dic2_edges = mol2_graph.edges(data=True)

        self.assertEqual(True, dic1_nodes == dic2_nodes)
        self.assertEqual(True, dic2_edges == dic2_edges)
        
        nm = iso.categorical_node_match(['fname_comp','ID'],['noname',-1])
        em = iso.categorical_edge_match(['strict_flag','similarity'],[False,-1.0])

        self.assertEqual(True, nx.is_isomorphic(graph, mol2_graph , node_match=nm, edge_match=em))
        
        
    def test_mcs(self):
        f = open('test/basic/MCS.pickle','rb')
        data = pickle.load(f)
        data_no_hydrogens = data[0]
        data_hydrogens = data[1]
        
        db = self.inst
        
        nohyds = {}
        hyds = {}
        
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        
        for i in range(0,db.nums()):
            for j in range(i+1,db.nums()):
                MCS_no_hyds = MCS.getMapping(db[i].getMolecule(), db[j].getMolecule())
                MCS_hyds = MCS.getMapping(db[i].getMolecule(), db[j].getMolecule(), hydrogens=True)
                if py_ver < 3:
                    nohyds[(i,j)] = MCS_no_hyds
                    hyds[(i,j)] = MCS_hyds
                else:
                    nohyds[(i,j)] = list(MCS_no_hyds)
                    hyds[(i,j)] = list(MCS_hyds)
                    
        self.assertEqual(True, nohyds == data_no_hydrogens)
        self.assertEqual(True, hyds  == data_hydrogens)

    

if __name__ =='__main__':
    unittest.main()
    
