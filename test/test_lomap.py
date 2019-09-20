import unittest
from lomap.mcs import MCS
from rdkit import RDLogger,Chem
from lomap.dbmol import DBMolecules
import multiprocessing
import subprocess
import math
import argparse


def executable():
    return '/home/mark/gui/trunk/buildDependencies/Python/linux-x86_64/bin/lomap'

def isclose(a,b):
    return (abs(a-b)<1e-5)

class TestLomap(unittest.TestCase):

    def test_insufficient_arguments(self):
        cmd = [executable()]
        error_string = b'error: the following arguments are required: directory'
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        assert(error_string in stderr)

    def test_mcsr(self):
        # MolA, molB, 3D?, max3d, mcsr, atomic_number_rule
        data=[ ('test/transforms/phenyl.sdf','test/transforms/toluyl.sdf', False, 1000, math.exp(-0.1 * (6 + 7 - 2*6)), 1) ,
               ('test/transforms/phenyl.sdf','test/transforms/chlorophenyl.sdf', False, 1000, math.exp(-0.1 * (6 + 7 - 2*6)), 1) ,
               ('test/transforms/toluyl.sdf','test/transforms/chlorophenyl.sdf', False, 1000, 1, math.exp(-0.1 * 0.5)) ,
               ('test/transforms/toluyl.sdf','test/transforms/chlorophenol.sdf', False, 1000, math.exp(-0.1 * (7 + 8 - 2*7)), math.exp(-0.1 * 0.5)),
               ('test/transforms/phenyl.sdf','test/transforms/napthyl.sdf', False, 1000, math.exp(-0.1 * (8 + 12 - 2*8)), 1),
               ('test/transforms/chlorophenyl.sdf','test/transforms/fluorophenyl.sdf', False, 1000, 1, math.exp(-0.1 * 0.5 )),
               ('test/transforms/chlorophenyl.sdf','test/transforms/bromophenyl.sdf', False, 1000, 1, math.exp(-0.1 * 0.15)),
               ('test/transforms/chlorophenyl.sdf','test/transforms/iodophenyl.sdf', False, 1000, 1, math.exp(-0.1 * 0.35)),

               # Compare with and without 3D pruning
               ('test/transforms/chlorophenyl.sdf','test/transforms/chlorophenyl2.sdf', False, 1000, 1, 1),
               ('test/transforms/chlorophenyl.sdf','test/transforms/chlorophenyl2.sdf', False, 2, math.exp(-0.1 * (7 + 7 - 2*6)), 1) ,

               # Compare with and without 3D matching
               ('test/transforms/chlorotoluyl1.sdf','test/transforms/chlorotoluyl2.sdf', False, 1000, 1, 1),
               ('test/transforms/chlorotoluyl1.sdf','test/transforms/chlorotoluyl2.sdf', True, 1000, 1, math.exp(-0.05 * 2)) 
            ]


        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)

        for d in data:
            mola = Chem.MolFromMolFile(d[0], sanitize=False, removeHs=False)
            molb = Chem.MolFromMolFile(d[1], sanitize=False, removeHs=False)
            MC = MCS(mola, molb, argparse.Namespace(time=20, verbose='info', max3d=d[3], threed=d[2]))
            mcsr = MC.mcsr()
            mncar = MC.mncar()
            atnum = MC.atomic_number_rule()
            strict = MC.tmcsr(strict_flag=True)
            loose = MC.tmcsr(strict_flag=False)

            assert(isclose(mcsr,d[4]))
            assert(isclose(atnum,d[5]))

    # Check iter and next
    def test_iter_next(self):
        inst = DBMolecules('test/basic/', parallel=1, verbose='off', output=False, time=20, ecrscore=0.0, name='out',
                           display=False, max=6, cutoff=0.4, radial=False, hub=None)
        for i in range(0, inst.nums()):
            inst.next()
        with self.assertRaises(StopIteration):
            inst.next()


    # Check get, set and add
    def test_get_set_add(self):
        inst = DBMolecules('test/basic/', parallel=1, verbose='off', output=False, time=20, ecrscore=0.0, name='out',
                           display=False, max=6, cutoff=0.4, radial=False, hub=None)
        with self.assertRaises(IndexError):
            inst.__getitem__(inst.nums()+1)
        with self.assertRaises(IndexError):
            inst.__setitem__(inst.nums()+1,inst[1])
        with self.assertRaises(ValueError):
            inst.__setitem__(0, 'no_mol_obj')
        with self.assertRaises(ValueError):
            inst.__add__('no_mol_obj')


    # Check serial and parallel mode
    def test_serial_parallel(self):
        db = DBMolecules('test/basic')
        s_strict, s_loose = db.build_matrices()
        db.options.paralell = multiprocessing.cpu_count()
        p_strict, p_loose = db.build_matrices()
        
        assert (all(s_strict == p_strict))
        assert (all(s_loose == p_loose))


    def test_read_mol2_files(self):
        db = DBMolecules('test/basic')
        db.options.directory = 'test/'
        with self.assertRaises(IOError):
            db.read_molecule_files()


    # Test which heterocycles I can grow (growing off a phenyl)
    # Test by Max indicates that growing complex heterocycles tends
    # to fail, so only allow growing phenyl, furan and pyrrole
    def test_heterocycle_scores(self):
        testdata=[('phenylfuran.sdf',1),
                 ('phenylimidazole.sdf',math.exp(-0.1*4)),
                 ('phenylisoxazole.sdf',math.exp(-0.1*4)),
                 ('phenyloxazole.sdf',math.exp(-0.1*4)),
                 ('phenylpyrazole.sdf',math.exp(-0.1*4)),
                 ('phenylpyridine1.sdf',math.exp(-0.1*4)),
                 ('phenylpyridine2.sdf',math.exp(-0.1*4)),
                 ('phenylpyrimidine.sdf',math.exp(-0.1*4)),
                 ('phenylpyrrole.sdf',1),
                 ('phenylphenyl.sdf',1)]
        parent=Chem.MolFromMolFile('test/transforms/phenyl.sdf',sanitize=False, removeHs=False)
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        for d in testdata:
            comp=Chem.MolFromMolFile('test/transforms/'+d[0],sanitize=False, removeHs=False)
            MC=MCS(parent,comp)
            self.assertEqual(MC.heterocycles_rule(penalty=4),d[1],'Fail on heterocycle check for '+d[0])

    # Tests by Max indicate that growing a sulfonamide all in one go is 
    # dodgy, so disallow it
    def test_sulfonamide_scores(self):
        testdata=[('cdk2_lig11.sdf',math.exp(-0.1*4)),
                 ('cdk2_lig1.sdf',1),
                 ('cdk2_lig2.sdf',math.exp(-0.1*4)),
                 ('cdk2_lig13.sdf',1),
                 ('cdk2_lig14.sdf',1),
                 ('cdk2_lig15.sdf',1) ]
        parent=Chem.MolFromMolFile('test/transforms/cdk2_lig16.sdf',sanitize=False, removeHs=False)
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        for d in testdata:
            comp=Chem.MolFromMolFile('test/transforms/'+d[0],sanitize=False, removeHs=False)
            MC=MCS(parent,comp)
            self.assertEqual(MC.sulfonamides_rule(penalty=4),d[1],'Fail on sulfonamide check for '+d[0])

    # Test to check symmetry equivalence by matching atomic numbers where possible
    def test_symmetry_matchheavies(self):
        mol1 = Chem.MolFromMolFile('test/transforms/chlorophenol.sdf',sanitize=False, removeHs=False)
        mol2 = Chem.MolFromMolFile('test/transforms/chlorophenyl.sdf',sanitize=False, removeHs=False)
        mol3 = Chem.MolFromMolFile('test/transforms/chlorophenyl2.sdf',sanitize=False, removeHs=False)
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        MCS1 = MCS(mol1,mol2)
        MCS2 = MCS(mol2,mol3)
        MCS3 = MCS(mol1,mol3)
        self.assertEqual(MCS1.mcs_mol.GetNumHeavyAtoms(),9)
        self.assertEqual([int(at.GetProp('to_moli')) for at in MCS1.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 9]);
        self.assertEqual([int(at.GetProp('to_molj')) for at in MCS1.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 8]);
        self.assertEqual([int(at.GetProp('to_moli')) for at in MCS2.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 8]);
        self.assertEqual([int(at.GetProp('to_molj')) for at in MCS2.mcs_mol.GetAtoms()],[4, 5, 0, 1, 2, 3, 7, 6, 8]);
        self.assertEqual([int(at.GetProp('to_moli')) for at in MCS3.mcs_mol.GetAtoms()],[4, 5, 0, 1, 2, 3, 7, 6, 9]);
        self.assertEqual([int(at.GetProp('to_molj')) for at in MCS3.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 8]);

    
    # Test to check symmetry equivalence by matching 3D coordinates rather than atomic numbers
    def test_symmetry_match3d(self):
        mol1 = Chem.MolFromMolFile('test/transforms/chlorophenol.sdf',sanitize=False, removeHs=False)
        mol2 = Chem.MolFromMolFile('test/transforms/chlorophenyl.sdf',sanitize=False, removeHs=False)
        mol3 = Chem.MolFromMolFile('test/transforms/chlorophenyl2.sdf',sanitize=False, removeHs=False)
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        MCS1 = MCS(mol1,mol2,options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=True))
        MCS2 = MCS(mol2,mol3,options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=True))
        MCS3 = MCS(mol1,mol3,options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=True))
        self.assertEqual(MCS1.mcs_mol.GetNumHeavyAtoms(),9)
        # MCS1 and MCS2 are the same as in the matchheavies case, but MCS3 gives a diffrent answer
        self.assertEqual([int(at.GetProp('to_moli')) for at in MCS1.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 9]);
        self.assertEqual([int(at.GetProp('to_molj')) for at in MCS1.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 8]);
        self.assertEqual([int(at.GetProp('to_moli')) for at in MCS2.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 8]);
        self.assertEqual([int(at.GetProp('to_molj')) for at in MCS2.mcs_mol.GetAtoms()],[4, 5, 0, 1, 2, 3, 7, 6, 8]);
        self.assertEqual([int(at.GetProp('to_moli')) for at in MCS3.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 8, 6, 9]);
        self.assertEqual([int(at.GetProp('to_molj')) for at in MCS3.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 8]);

    # Test to check removing atoms from the MCS when the 3D coords are too far apart
    def test_clip_on_3d(self):
        mol1 = Chem.MolFromMolFile('test/transforms/chlorophenyl.sdf',sanitize=False, removeHs=False)
        mol2 = Chem.MolFromMolFile('test/transforms/chlorophenyl2.sdf',sanitize=False, removeHs=False)
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        MCS1 = MCS(mol1,mol2,options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=True))
        MCS2 = MCS(mol1,mol2,options=argparse.Namespace(time=20, verbose='info', max3d=2, threed=True))
        self.assertEqual(MCS1.mcs_mol.GetNumHeavyAtoms(),9)
        self.assertEqual(MCS2.mcs_mol.GetNumHeavyAtoms(),8)

    # Test disallowing turning a methyl group (or larger) into a ring atom
    def test_transmuting_methyl_into_ring_rule(self):
        testdata=[('phenyl.sdf','toluyl3.sdf',1),
                 ('toluyl3.sdf','chlorotoluyl1.sdf',1),
                 ('toluyl3.sdf','phenylfuran.sdf',math.exp(-0.1*4)),
                 ('toluyl3.sdf','phenylpyridine1.sdf',math.exp(-0.1*4)),
                 ('phenyl.sdf','phenylfuran.sdf',1),
                 ('phenyl.sdf','phenylpyridine1.sdf',1),
                 ('chlorophenol.sdf','phenylfuran.sdf',1)
                 ]
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        for d in testdata:
            parent=Chem.MolFromMolFile('test/transforms/'+d[0],sanitize=False, removeHs=False)
            comp=Chem.MolFromMolFile('test/transforms/'+d[1],sanitize=False, removeHs=False)
            MC=MCS(parent,comp)
            self.assertEqual(MC.transmuting_methyl_into_ring_rule(penalty=4),d[2],'Fail on transmuting-methyl-to-ring check for '+d[0]+' '+d[1])

    # Test penalising hybridization changes
    def test_hybridization_rule(self):
        testdata=[('napthyl.sdf','tetrahydronaphthyl.sdf',math.exp(-0.1 * 4))
                 ]
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        for d in testdata:
            parent=Chem.MolFromMolFile('test/transforms/'+d[0],sanitize=False, removeHs=False)
            comp=Chem.MolFromMolFile('test/transforms/'+d[1],sanitize=False, removeHs=False)
            MC=MCS(parent,comp)
            assert(isclose(MC.hybridization_rule(),d[2]))

    # Test disallowing turning a ring into an incompatible ring
    def test_transmuting_ring_sizes_rule(self):
        testdata=[('phenyl.sdf','phenylcyclopropyl.sdf',1),
                 ('toluyl.sdf','phenylcyclopropyl.sdf',1),  # Disallowed by test_transmuting_methyl_into_ring_rule instead
                 ('phenylcyclopropyl.sdf','phenylcyclobutyl.sdf',0),
                 ('phenylcyclopropyl.sdf','phenylcyclopentyl.sdf',0),
                 ('phenylcyclopropyl.sdf','phenylcyclononyl.sdf',0),
                 ('phenylcyclobutyl.sdf','phenylcyclopentyl.sdf',0),
                 ('phenylcyclobutyl.sdf','phenylcyclononyl.sdf',0),
                 ('phenylcyclopentyl.sdf','phenylcyclononyl.sdf',1)
                 ]
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        for d in testdata:
            parent=Chem.MolFromMolFile('test/transforms/'+d[0],sanitize=False, removeHs=False)
            comp=Chem.MolFromMolFile('test/transforms/'+d[1],sanitize=False, removeHs=False)
            MC=MCS(parent,comp)
            self.assertEqual(MC.transmuting_ring_sizes_rule(),d[2],'Fail on transmuting-ring-size check for '+d[0]+' '+d[1])

    # Test getting the mapping string out of the MCS
    def test_mapping_string_heavy(self):
        testdata=[('phenyl.sdf','toluyl3.sdf',"0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7"),
                 ('toluyl2.sdf','chlorotoluyl1.sdf',"0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:8,8:9"),
                 ('toluyl3.sdf','phenylfuran.sdf',"0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7")
                 ]
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        for d in testdata:
            parent=Chem.MolFromMolFile('test/transforms/'+d[0],sanitize=False, removeHs=False)
            comp=Chem.MolFromMolFile('test/transforms/'+d[1],sanitize=False, removeHs=False)
            MC=MCS(parent,comp)
            self.assertEqual(MC.heavy_atom_match_list(), d[2], 'Fail on heavy atom match list for '+d[0]+' '+d[1])

    # Test getting the mapping string including hydrogens out of the MCS
    def test_mapping_string_hydrogen(self):
        testdata=[('phenyl.sdf','toluyl3.sdf',"0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:9,9:10,10:8,11:11,12:17,13:12,14:13,15:14,16:15,17:16"),
                 ('toluyl2.sdf','chlorotoluyl1.sdf',"0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:8,8:9,9:10,10:7,11:11,12:12,13:13,14:14,15:15,16:16,17:17,18:18,19:19,20:20"),
                 ('toluyl3.sdf','phenylfuran.sdf',"0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7,9:13,10:14,11:15,12:17,13:18,14:19,15:20,16:21,17:16"),
                 ('toluyl.sdf','phenylmethylamino.sdf',"0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:8,9:10,10:11,11:12,12:13,13:14,14:15,15:16,16:17,17:18,18:9,19:19,20:20")
                 ]
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        for d in testdata:
            parent=Chem.MolFromMolFile('test/transforms/'+d[0],sanitize=False, removeHs=False)
            comp=Chem.MolFromMolFile('test/transforms/'+d[1],sanitize=False, removeHs=False)
            MC=MCS(parent,comp)
            self.assertEqual(MC.all_atom_match_list(), d[2], 'Fail on all-atom match list for '+d[0]+' '+d[1])

    # Test to check correct handling of chirality
    def test_chirality_handling(self):
        testdata=[('Chiral1R.sdf','Chiral1S.sdf',6),
                  ('Chiral1R.sdf','Chiral2R.sdf',7),
                  ('Chiral1S.sdf','Chiral2R.sdf',6),
                  ('Chiral3RS.sdf','Chiral3SS.sdf',11),
                  ('Chiral3SR.sdf','Chiral3SS.sdf',10),
                  ('Chiral3SR.sdf','Chiral3RS.sdf',9),
                  ('Chiral4RR.sdf','Chiral4RS.sdf',5),
                  ('RingChiralR.sdf','RingChiralS.sdf',6),
                  ('SpiroR.sdf','SpiroS.sdf',6)
                ]
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.INFO)
        for d in testdata:
            parent=Chem.MolFromMolFile('test/chiral/'+d[0],sanitize=False, removeHs=False)
            comp=Chem.MolFromMolFile('test/chiral/'+d[1],sanitize=False, removeHs=False)
            MC=MCS(parent,comp, argparse.Namespace(time=20, verbose='info', max3d=5, threed=True))
            self.assertEqual(MC.mcs_mol.GetNumHeavyAtoms(), d[2], 'Fail on chiral MCS size for '+d[0]+' '+d[1])

    # Test to check correct trimming of rings when 3D coordinate matching is used
    def test_ring_trimming_on_3d_match(self):
        parent=Chem.MolFromMolFile('test/transforms/phenylcyclopentylmethyl1.sdf',sanitize=False, removeHs=False)
        comp=Chem.MolFromMolFile('test/transforms/phenylcyclopentylmethyl2.sdf',sanitize=False, removeHs=False)
        MC=MCS(parent,comp, argparse.Namespace(time=20, verbose='info', max3d=2, threed=True))
        self.assertEqual(MC.mcs_mol.GetNumHeavyAtoms(), 9, 'Fail on ring trim on 3D match')

if __name__ == '__main__':
    unittest.main()
            
