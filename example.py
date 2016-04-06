# Import lomap
import lomap
    
# Create the molecule database by using .mol2 files
# The DBMolecule class must be created with a valid
# directory name
        
db_mol = lomap.DBMolecules('lomap/test/basic')
    
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

# Output the MCS in a .png file

MC.draw_mcs()
