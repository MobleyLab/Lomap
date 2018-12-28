import copy


class Molecule(object):
    """
    This Class stores the molecule object,
    its identification and the molecule name.
    Based on the selected toolkit the molecule
    object can be an OE molecule or a Rdkit
    molecule
    """

    __counter = 0

    def __init__(self, molecule, molname):
        """
        Initialization class function

        Parameters
        ----------
        molecule :molecule object
           the molecule

        molname : str
           the molecule name

        """
        # The variable __molecule saves the current molecule object
        self.__molecule = molecule

        # The variable __ID saves the molecule identification number
        self.__ID = Molecule.__counter

        Molecule.__counter += 1

        # The variable __name saves the molecule identification name
        self.__name = molname

    def __del__(self):
        Molecule.__counter -= 1

    @property
    def id(self):
        """
        Get the molecule ID number


        Returns
        -------
           : int
           the molecule ID number

        """
        return self.__ID

    @property
    def molecule(self):
        """
        Get the molecule object

        Returns
        -------
        mol_copy :molecule object
           The copy of the molecule

        """
        mol_copy = copy.deepcopy(self.__molecule)
        return mol_copy

    @property
    def name(self):
        """
        Get the molecule file name

        Returns
        -------
           : str
           the molecule string file name

        """
        return self.__name

