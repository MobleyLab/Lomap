"""
Lomap2
======

Alchemical free energy calculations hold increasing promise as an aid to drug 
discovery efforts. However, applications of these techniques in discovery 
projects have been relatively few, partly because of the difficulty of planning 
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an 
automated algorithm to plan efficient relative free energy calculations between 
potential ligands within a substantial of compounds.

Authors: Gaetano Calabro' <gcalabro@uci.edu> 
         David Mobley     <dmobley@uci.edu>


Licence: LGPL

URL: https://github.com/nividic/Lomap


Using
-----
      Just write in Python

      # Import lomap
      import lomap

      # Create the molecule database by using .mol2 files
      # The DBMolecule class must be created with a valid
      # directory name
    
      db_mol = lomap.DBMolecules('lomap/test/basic/')

      # Generate the strict and loose syimmetric similarity 
      # score matrices
      
      strict, loose = db_mol.build_matrices()

      # Convert the matrices in standard numpy matrices
      
      strict_numpy = strict.to_numpy_2D_array()
      loose_numpy = loose.to_numpy_2D_array()

      # Networkx graph generation based on the similarity 
      # score matrices
      
      nx_graph = db_mol.build_graph() 
            
      # Calculate the Maximum Common Subgraph (MCS) between 
      # the first two molecules in the molecule database

      MC = lomap.MCS(db_mol[0].getMolecule(), db_mol[1].getMolecule())

      # Get the map between atom indices 
      mcs_map = MC.getMap()

      # Output the MCS in a .png file

      MC.draw_mcs()
"""



from dbmol import DBMolecules
from dbmol import SMatrix
from dbmol import Molecule
from mcs import MCS

del dbmol
del graphgen
del mcs
