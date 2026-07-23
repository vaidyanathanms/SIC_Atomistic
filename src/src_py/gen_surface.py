from ase.build import bcc110, bcc100, bcc111
from ase.io import read, write
from ase.build import surface
import numpy as np

# ===============================
# USER INPUTS:
# ===============================

# Target lateral dimensions (in Angstroms)
#a_li     = -1 # -1: default, or give exact input
a_li     = 3.4824 # average of lx and ly divided by number of x and y
# repeats for large_bare_surf simulations
target_x = 81.2662 # Obtain from averaged bulk simulations 
target_y = 82.1931

nx = 1; ny = 2 # number of atoms in surface unit cells in X,Y
nz_repeats = 4 # repetitions along surface normal


baseunits = f'{nx}{ny}{nz_repeats}'

# ===============================
# STEP 1: Build Lithium (110) Surface using ASE
# ===============================
# Create (110) surface: size = (nx, ny, nz)
# Note ny should be even
if a_li == -1:
    unit_slab = bcc110('Li', size=(nx, ny, nz_repeats),vacuum=0.0, orthogonal=True)
else:
    unit_slab = bcc110('Li', size=(nx, ny, nz_repeats), a=a_li, vacuum=0.0, orthogonal=True)
write(f'../../InputStructures/inpcoord_files/li_surface/Li_{baseunits}_surface.pdb',unit_slab)
write(f'../../InputStructures/inpcoord_files/li_surface/Li_{baseunits}_surface.xyz',unit_slab)


# Get unit cell dimensions
cell_x, cell_y, cell_z = unit_slab.get_cell().lengths()
print(cell_x, cell_y, cell_z)

# Calculate replication factors
repeat_x = int(target_x / cell_x) + 1
repeat_y = int(target_y / cell_y) + 1

print(f"Repeating Lithium slab {repeat_x} times in X and {repeat_y} times in Y.")

total_li = nx*ny*repeat_x*repeat_y*nz_repeats 
print(f"Total lithium atoms is {total_li}")

# ===============================
# STEP 2: Replicate slab
# ===============================
li_large = unit_slab.repeat((repeat_x, repeat_y, 1))

# Optional: Add some vacuum along Z if needed
vacuum_z = 0.0  # No vacuum for the slab itself
li_large.center(vacuum=vacuum_z, axis=2)

# ===============================
# STEP 3: Save the expanded slab
# ===============================
output_file = f'../../InputStructures/inpcoord_files/li_surface/lithium_large_{repeat_x}_{repeat_y}_{nz_repeats}.pdb'
write(output_file, li_large)

print(f"Expanded Lithium slab written to {output_file}.")
print(f"New cell dimensions: {li_large.get_cell().lengths()}")
