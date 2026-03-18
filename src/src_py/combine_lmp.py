#!/usr/bin/env python3

from lmp_define import LammpsData
import aux_functs as af

# Main file for combining data files
# one template file per species
file_specs = [
    ("../../InputStructures/inp_coordfiles/lithium_large_edited.data", 1),
    ("../../InputStructures/inp_coordfiles/60PEO_Optimized_CH3terminated_editedterminal_edited.data", 26),
    ("../../InputStructures/inp_coordfiles/Li_Atom_edited.data", 1284),
    ("../../InputStructures/inp_coordfiles/co3_2minus_edited.data", 642),
    ("../../InputStructures/inp_coordfiles/V27M5_150T_edited.data", 120),
    ("../../InputStructures/inp_coordfiles/Li_Atom_edited.data",600)
]

system = af.combine_lammps_system(file_specs, spacing=5.0)

af.write_lammps_data(filename="../../InputStructures/combined/combined_system.data", system=system,\
                     box=((0.0, 89.0), (0.0, 89.7), (0.0, 178.5)),\
                     write_atoms='../../InputStructures/combined/all_atoms.data')

print("Wrote combined_system.data")


