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

    # FIXME: we ned to "know" that molecule is a rdkit.Chem.rdchem.Mol
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


class RDKitMolReader(object):
    """
    Read molecular structure information from files with RDKit.
    """

    def __init__(self):
        self.mol_readers = {'sdf': self._sdf_supplier,
                            'mol': rdchem.MolFromMolFile,
                            'mol2': rdchem.MolFromMol2File,
                            'pdb': rdchem.MolFromPDBFile,
                            'tpl': rdchem.MolFromTPLFile}

    def _sdf_supplier(self, *args, **kwargs):
        mols = []

        supplier = rdchem.SDMolSupplier(*args, **kwargs)

        for mol in supplier:
            if mol:
                mols.append(mol)

        return mols

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
            mols = mol_reader(filename, sanitize=False, removeHs=False)
            all_mols = []

            # FIXME: this is ugly
            if hasattr(mols, '__iter__'):
                for mol in mols:
                    #rdchem.SanitizeMol(mol,
                    #             rdchem.SANITIZE_ALL^rdchem.SANITIZE_KEKULIZE)
                    all_mols.append(Molecule(mol, '', ''))
            else:
                all_mols = [Molecule(mols, '', '')]

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

    for filename in sys.argv[1:]:
        mols = read_molecules(filename)

        if mols:
            all_mols.extend(mols)

    if not all_mols:
        logger.error('no molecular structures found')
        sys.exit(1)

    for mol in all_mols:
        print mol.molecule.GetNumAtoms(), mol.name, mol.ID
