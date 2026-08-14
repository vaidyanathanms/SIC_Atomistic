#!/usr/bin/env python3

from lmp_define import LammpsData
import aux_functs as af


# Main file for combining data files
# one template file per species
# Use the atoms from equilibrated PDB file
# NOTE VERY IMPORTANT: THE ORDER IN PDB FILE SHOULD BE SAME AS THAT IN
# FILE_SPECS ARRAY

head_dir =
'/home/vaidyams/all_codes/files_interface/InputStructures/inpcoord_files'

xmin = -0.5; xmax = 89.0
ymin = -0.5; ymax = 89.7
zmin =  0.0; zmax = 178.5

file_specs = [
    (head_dir+"/li_surface/lithium_large_edited.data", 1),
    (head_dir+"/peo_polymer/PEO_Optimized_CH3terminated_editedterminal_edited.data", 26),
    (head_dir+"/li_monomer/Li_Atom_edited.data", 1284),
    (head_dir+"/co32m_monomer/co3_2minus_edited.data", 642),
    (head_dir+"/vecmtfsi_polymer/InputStructures/inp_coordfiles/V27M5_150T_edited.data", 120),
    (head_dir+"/li_monomer/InputStructures/inp_coordfiles/Li_Atom_edited.data",600)
]


system = af.combine_lammps_system(file_specs, equil_pdb_file)

af.write_lammps_data(filename="combined_system.data",system=system,\
                     box=((xmin,xmax),(ymin,ymax),(zmin,zmax)),\
                     write_atoms='all_atoms.data')

print("Wrote combined_system.data")






