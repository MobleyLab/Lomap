# *****************************************************************************
# Lomap: A toolkit to plan alchemical relative binding affinity calculations
# Copyright 2015 - 2019  UC Irvine and the Authors
#
# Authors: Dr Gaetano Calabro' and Dr David Mobley
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the MIT License along with this library;
# if not, see https://opensource.org/licenses/MIT
# *****************************************************************************

import importlib
import lomap

try:
    import openeye
    if not openeye.OEChemIsLicensed():
        raise ImportError
    from openeye import oechem
    HAS_OE = True
except ImportError:
    HAS_OE = False

try:
    from rdkit import Chem
    HAS_RDK = True
except ImportError:
    HAS_RDK = False

if HAS_RDK:
    DEFAULT = "RDK"
    from .rdk import read_molecules
elif HAS_OE:
    DEFAULT = "OE"
    from .oe import read_molecules
else:
    raise ImportError("No Cheminformatics toolkit has been found")


def set_default(toolkit):
    """
    Set the default toolkit

    Parameters
    ----------
    toolkit: str
        the toolkit string name. Currently supported toolkit names are:
        OE and RDK

    Returns
    -------
    DEFAULT : str
        the new default toolkit
    """

    if toolkit == "OE":
        if HAS_OE:
            lomap.toolkits.DEFAULT = "OE"
            mod = importlib.import_module("lomap.toolkits.oe")
        else:
            raise ValueError("OE toolkit is not installed")

    elif toolkit == "RDK":
        if toolkit == "RDK":
            if HAS_RDK:
                lomap.toolkits.DEFAULT = "RDK"
                mod = importlib.import_module("lomap.toolkits.rdk")
            else:
                raise ValueError("RDK toolkit is not installed")
    else:
        raise ValueError("The selected toolkit is not supported: {} \n"
                         "OE and RDK are currently supported".format(toolkit))

    lomap.toolkits.read_molecules = mod.__dict__['read_molecules']

    return lomap.toolkits.DEFAULT
