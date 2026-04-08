from ase.build import bcc110, bcc100, bcc111
from ase.io import read, write
from ase.build import surface

# ===============================
# USER INPUTS:
# ===============================

# Target lateral dimensions (in Angstroms)
target_x = 74.0 # Obtain from excel sheet for x-y dimensions
target_y = 74.0

# ===============================
# STEP 1: Build Lithium (110) Surface using ASE
# ===============================
# Create (110) surface: size = (nx, ny, nz)
# Note ny should be even
<<<<<<< HEAD
unit_slab = bcc110('Li', size=(1, 2, 12), vacuum=0.0, orthogonal=True)
=======
unit_slab = bcc110('Li', size=(1, 2, 16), vacuum=0.0, orthogonal=True)
>>>>>>> origin/main
write('Li_124_surface.pdb',unit_slab)
write('Li_124_surface.xyz',unit_slab)


# Get unit cell dimensions
cell_x, cell_y, cell_z = unit_slab.get_cell().lengths()

# Calculate replication factors
repeat_x = int(target_x / cell_x) + 1
repeat_y = int(target_y / cell_y) + 1

print(f"Repeating Lithium slab {repeat_x} times in X and {repeat_y} times in Y.")

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
output_file = '../../InputStructures/inp_coordfiles/lithium_large.pdb'
write(output_file, li_large)

print(f"Expanded Lithium slab written to {output_file}.")
print(f"New cell dimensions: {li_large.get_cell().lengths()}")
