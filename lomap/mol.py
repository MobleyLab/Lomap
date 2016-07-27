r"""
Read molecular structure information with either RDKit or OpenBabel and store
this.
"""


import os
from logger import logger

import rdkit.Chem as rdchem
from rdkit import rdBase



class Molecule(object):
    """
    Simple container to keep basic information about a molecule.  The
    molecule itself is a rdkit.Chem.rdchem.Mol.

    Seems mostly redundant at the moment.
    """

    __slots__ = ['molecule', 'name', 'ID']

    # FIXME: we need to "know" that molecule is a rdkit.Chem.rdchem.Mol
    def __init__(self, rdmol, molname, molid=''):
        """
        :param rdmol: RDKit molecule
        :type rdmol: rdkit.Chem.rdchem.Mol
        :param molname: molecule name, can be the filename or less useful
                        strings, empty, etc.
        :type molname: str
        :param molid: the ID for the molecule
        :type molid: str
        """

        self.molecule = rdmol
        self.name = molname
        self.ID = molid


class MorphPair(object):
    """
    Simple container to keep basic information about a morph pair.
    """

    #__slots__ = []

    def __init__(self, mol0, mol1):
        """
        :param mol0: molecule for state 0
        :type mol0: Molecule
        :param mol1: molecule for state 1
        :type mol1: Molecule
        """

        self.mol0 = mol0
        self.mol1 = mol1
        self.mapping = None
        self.strict_score = None
        self.loose_score = None


def sdf_supplier(filename, *args, **kwargs):
    """
    Wrapper for RDKit's SDMolSupplier.

    :param filename: name of file with molecules
    :type filename: str
    :returns: all molecules found
    :rtype: list of rdkit.Chem.rdchem.Mol
    """

    mols = []
    supplier = rdchem.SDMolSupplier(filename, *args, **kwargs)

    for mol in supplier:
        if mol:
            mols.append(mol)

    return mols

tipos_mol = '@<TRIPOS>MOLECULE'

def itermol2(filename):
    """
    Mol2 file generator to return each individual Tripos molecule.

    :param filename: name of file with molecules
    :type filename: str
    :returns: next molecule string
    :rtype: str
    """
    
    lines = []
    first = True

    with open(filename, 'rb') as mol2_file:
        for line in mol2_file:
            if line.startswith(tipos_mol):
                if not first:
                    yield ''.join(lines)
                    lines = []
                else:
                    first = False

            lines.append(line)

    yield ''.join(lines)

def mol2_supplier(filename, *args, **kwargs):
    """
    Simplistic mol2 supplier for RDKit.

    :param filename: name of file with molecules
    :type filename: str
    :returns: all molecules found
    :rtype: list of rdkit.Chem.rdchem.Mol
    """

    mols = []

    for mol in itermol2(filename):
        mols.append(rdchem.MolFromMol2Block(mol, *args, **kwargs))

    return mols

def fake_supplier(func):
    """
    A fake supplier for molecule file formats containing only single molecules.
    Implemented as closure.

    :param func: name of file with molecules
    :type func: function
    :returns: fake generator that yields only on molecule
    :rtype: function
    """

    def iterfile(filename, *args, **kwargs):
        """Fake generator."""
        yield func(filename, *args, **kwargs)

    return iterfile

class RDKitMolReader(object):
    """
    Read molecular structure information from files with RDKit.
    """

    mol_readers = {'sdf': sdf_supplier,
                   'mol2': mol2_supplier,
                   'mol': fake_supplier(rdchem.MolFromMolFile),
                   'pdb': fake_supplier(rdchem.MolFromPDBFile),
                   'tpl': fake_supplier(rdchem.MolFromTPLFile)}

    def read_molecules(self, filename):
        """
        Read molecules from a file with RDKit and convert each structure
        to a RDKit molecule.

        :param filename: name of file that contains molecule(s)
        :type filename: str or None when failure
        """

        file_ext = os.path.splitext(filename)[1][1:]

        try:
            mol_reader = self.mol_readers[file_ext]
        except KeyError:
            logger.warn('cannot guess file format of %s' % filename)
            return None
        else:
            # FIXME: error handling
            mols = mol_reader(filename, sanitize=False, removeHs=False)
            all_mols = []

            for mol in mols:
                #rdchem.SanitizeMol(mol,
                #             rdchem.SANITIZE_ALL^rdchem.SANITIZE_KEKULIZE)
                all_mols.append(Molecule(mol, '', ''))

        return all_mols

class OBMolReader(object):
    """
    Read molecular structure information from files with OpenBabel.
    """

    ob = __import__('openbabel')

    def read_molecules(self, filename):
        """
        Read molecules from a file with Openbabel and convert each structure
        to a RDKit molecule.  This is currently done by RDKit reparsing the
        string output from a Openbabel molecule in MDL mol format.

        :param filename: name of file that contains molecule(s)
        :type filename: str or None when failure
        """

        ob = self.ob
        conv = ob.OBConversion()
        fmt = conv.FormatFromExt(filename)

        if not fmt:
            logger.warn('cannot guess file format of %s' % filename)
            return None

        conv.SetInAndOutFormats(fmt.GetID(), 'mol2')

        # FIXME: openbabel may produce unwanted warnings, do we want these?
        errlev = ob.obErrorLog.GetOutputLevel()
        ob.obErrorLog.SetOutputLevel(0)

        obmol = ob.OBMol()
        ok = conv.ReadFile(obmol, filename)

        if not ok:
            logger.warn('cannot read molecule data from file %s' % filename)
            return None

        obmols = []
        obdata = []
        inchi = ob.OBConversion()        # RDKit may not have INCHI compiled in
        inchi.SetOutFormat('inchikey')

        while ok:
            obmols.append(obmol)
            obdata.append((obmol.GetTitle(),
                            inchi.WriteString(obmol).rstrip()))
            obmol = ob.OBMol()
            ok = conv.Read(obmol)

        ob.obErrorLog.SetOutputLevel(errlev)

        rdmols = []

        # FIXME: do we really want this?
        rdBase.DisableLog('rdApp.warning')

        # FIXME: build RDKit Mol from scratch?
        for obmol in obmols:
            mol_str = conv.WriteString(obmol)
            rdmol = rdchem.MolFromMol2Block(mol_str, sanitize=False,
                                            removeHs=False)
            rdmols.append(rdmol)

            rdBase.EnableLog('rdApp.warning')

        mols = []

        for rdmol, data in zip(rdmols, obdata):
            rdchem.SanitizeMol(rdmol, rdchem.SANITIZE_ALL^rdchem.SANITIZE_KEKULIZE)
            mols.append(Molecule(rdmol, *data))

        return mols


if __name__ == '__main__':
    import sys
    import logging

    from logger import LogFormatter


    log = True

    if log:
        hdlr = logging.StreamHandler(sys.stdout)
        hdlr.setFormatter(LogFormatter())
        logger.addHandler(hdlr)

    all_mols = []
    mol_reader = RDKitMolReader()

    for filename in sys.argv[1:]:
        mols = mol_reader.read_molecules(filename)

        if mols:
            all_mols.extend(mols)

    if not all_mols:
        logger.error('no molecular structures found')
        sys.exit(1)

    for mol in all_mols:
        print mol.molecule.GetNumAtoms(), mol.name, mol.ID
