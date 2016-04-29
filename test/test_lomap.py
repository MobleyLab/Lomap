import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

import unittest
from unittest import skipIf
from lomap.dbmol import DBMolecules
import argparse
import multiprocessing
import networkx as nx
import networkx.algorithms.isomorphism as iso

class TestLomap(unittest.TestCase):
    
    def setUp(self):
        self.inst = DBMolecules('test/basic/', time=20, parallel=1, verbose='off', output=False, name='out', display=False, max=6, cutoff=0.4)
    
    # Test class Inizialitazion

    # Check wrong number of argumetns or wrong options
    def test_InsufficientArgs(self):
        self.assertRaises(TypeError, self.inst.__init__)
       
    # Check passed input Types
    def test_TypeArgs(self):
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/nodir/')
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/basic', time=-1)
        self.assertRaises(SystemExit, self.inst.__init__, 'test/basic', time=-1.5)
        self.assertRaises(argparse.ArgumentTypeError, self.inst.__init__, 'test/basic', parallel=-1)
        self.assertRaises(SystemExit, self.inst.__init__, 'test/basic', parallel=-1.5)
        self.assertRaises(SystemExit, self.inst.__init__, 'test/basic', verbose='err_option')
        self.assertRaises(TypeError, self.inst.__init__, 'test/basic', output=-5.0)
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
    def test_graph(self):

        db =  DBMolecules('test/basic')
        
        strict,loose = db.build_matrices()
        graph = db.build_graph()
        
        mol2_graph = nx.read_gpickle("test/basic/molecules.gpikle")

        nm = iso.categorical_node_match(['fname_comp','ID'],['noname',-1])
        em = iso.categorical_edge_match(['strict_flag','similarity'],[False,-1.0])
        
        self.assertEqual(True, nx.is_isomorphic(graph, mol2_graph , node_match=nm, edge_match=em))
        
    


if __name__ =='__main__':
    unittest.main()
    
