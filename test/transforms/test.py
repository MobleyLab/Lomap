import unittest
import sys
import argparse
sys.path.insert(0,'../../')
from lomap import mcs 
from rdkit import Chem

class TestLomap(unittest.TestCase):

    # Test which heterocycles I can grow (growing off a phenyl)
    # Test by Max indicates that growing complex heterocycles tends
    # to fail, so only allow growing phenyl, furan and pyrrole
    def test_heterocycle_scores(self):
        testdata=[('phenylfuran.sdf',1),
                 ('phenylimidazole.sdf',0),
                 ('phenylisoxazole.sdf',0),
                 ('phenyloxazole.sdf',0),
                 ('phenylpyrazole.sdf',0),
                 ('phenylpyridine1.sdf',0),
                 ('phenylpyridine2.sdf',0),
                 ('phenylpyrimidine.sdf',0),
                 ('phenylpyrrole.sdf',1),
                 ('phenylphenyl.sdf',1)]
        parent=Chem.MolFromMolFile('phenyl.sdf',sanitize=False, removeHs=False)
        for d in testdata:
            comp=Chem.MolFromMolFile(d[0],sanitize=False, removeHs=False)
            MC=mcs.MCS(parent,comp)
            self.assertEqual(MC.heterocycles_rule(),d[1],'Fail on heterocycle check for '+d[0])

    # Tests by Max indicate that growing a sulfonamide all in one go is 
    # dodgy, so disallow it
    def test_sulfonamide_scores(self):
        testdata=[('cdk2_lig11.sdf',0),
                 ('cdk2_lig1.sdf',1),
                 ('cdk2_lig2.sdf',0),
                 ('cdk2_lig13.sdf',1),
                 ('cdk2_lig14.sdf',1),
                 ('cdk2_lig15.sdf',1) ]
        parent=Chem.MolFromMolFile('cdk2_lig16.sdf',sanitize=False, removeHs=False)
        for d in testdata:
            comp=Chem.MolFromMolFile(d[0],sanitize=False, removeHs=False)
            MC=mcs.MCS(parent,comp)
            self.assertEqual(MC.sulfonamides_rule(),d[1],'Fail on sulfonamide check for '+d[0])

    # Test to check symmetry equivalence by matching atomic numbers where possible
    def test_symmetry_matchheavies(self):
        mol1 = Chem.MolFromMolFile('chlorophenol.sdf',sanitize=False, removeHs=False)
        mol2 = Chem.MolFromMolFile('chlorophenyl.sdf',sanitize=False, removeHs=False)
        mol3 = Chem.MolFromMolFile('chlorophenyl2.sdf',sanitize=False, removeHs=False)
        MCS1 = mcs.MCS(mol1,mol2)
        MCS2 = mcs.MCS(mol2,mol3)
        MCS3 = mcs.MCS(mol1,mol3)
        self.assertEqual(MCS1.mcs_mol.GetNumHeavyAtoms(),9)
        self.assertEqual([int(at.GetProp('to_moli')) for at in MCS1.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 9]);
        self.assertEqual([int(at.GetProp('to_molj')) for at in MCS1.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 8]);
        self.assertEqual([int(at.GetProp('to_moli')) for at in MCS2.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 8]);
        self.assertEqual([int(at.GetProp('to_molj')) for at in MCS2.mcs_mol.GetAtoms()],[4, 5, 0, 1, 2, 3, 7, 6, 8]);
        self.assertEqual([int(at.GetProp('to_moli')) for at in MCS3.mcs_mol.GetAtoms()],[4, 5, 0, 1, 2, 3, 7, 6, 9]);
        self.assertEqual([int(at.GetProp('to_molj')) for at in MCS3.mcs_mol.GetAtoms()],[0, 5, 4, 3, 2, 1, 7, 6, 8]);

    
    # Test to check symmetry equivalence by matching 3D coordinates rather than atomic numbers
    def test_symmetry_match3d(self):
        mol1 = Chem.MolFromMolFile('chlorophenol.sdf',sanitize=False, removeHs=False)
        mol2 = Chem.MolFromMolFile('chlorophenyl.sdf',sanitize=False, removeHs=False)
        mol3 = Chem.MolFromMolFile('chlorophenyl2.sdf',sanitize=False, removeHs=False)
        MCS1 = mcs.MCS(mol1,mol2,options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=True))
        MCS2 = mcs.MCS(mol2,mol3,options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=True))
        MCS3 = mcs.MCS(mol1,mol3,options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=True))
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
        mol1 = Chem.MolFromMolFile('chlorophenyl.sdf',sanitize=False, removeHs=False)
        mol2 = Chem.MolFromMolFile('chlorophenyl2.sdf',sanitize=False, removeHs=False)
        MCS1 = mcs.MCS(mol1,mol2,options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=True))
        MCS2 = mcs.MCS(mol1,mol2,options=argparse.Namespace(time=20, verbose='info', max3d=2, threed=True))
        self.assertEqual(MCS1.mcs_mol.GetNumHeavyAtoms(),9)
        self.assertEqual(MCS2.mcs_mol.GetNumHeavyAtoms(),8)

    # Test disallowing turning a methyl group (or larger) into a ring atom
    def test_transmuting_methyl_into_ring_rule(self):
        testdata=[('phenyl.sdf','toluyl3.sdf',1),
                 ('toluyl3.sdf','chlorotoluyl1.sdf',1),
                 ('toluyl3.sdf','phenylfuran.sdf',0),
                 ('toluyl3.sdf','phenylpyridine1.sdf',0),
                 ('phenyl.sdf','phenylfuran.sdf',1),
                 ('phenyl.sdf','phenylpyridine1.sdf',1),
                 ('chlorophenol.sdf','phenylfuran.sdf',1)
                 ]
        for d in testdata:
            parent=Chem.MolFromMolFile(d[0],sanitize=False, removeHs=False)
            comp=Chem.MolFromMolFile(d[1],sanitize=False, removeHs=False)
            MC=mcs.MCS(parent,comp)
            self.assertEqual(MC.transmuting_methyl_into_ring_rule(),d[2],'Fail on transmuting-methyl-to-ring check for '+d[0]+' '+d[1])

if __name__ == '__main__':
    unittest.main()
            
