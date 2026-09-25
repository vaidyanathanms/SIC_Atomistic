#!/usr/bin/env python3

from lmp_define import LammpsData
import aux_functs as af
import os

# Main file for combining data files
# one template file per species
# Use the atoms from equilibrated PDB file

# NOTE VERY IMPORTANT: THE ORDER IN PDB FILE SHOULD BE SAME AS THAT IN
# FILE_SPECS ARRAY

#Inputs
equil_pdb_file = 'li_co32m_peo_vecmtfsi_fan_0pt06.pdb_FORCED'
if not os.path.exists(equil_pdb_file):
    raise RuntimeError(f"{equil_pdb_file} not found in {os.getcwd}")

head_dir = '/home/vaidyams/all_codes/files_interface/InputStructures/inpcoord_files'

xmin = -0.05; xmax = 83.5776000
ymin = -0.05; ymax = 83.7225700
zmin = -0.60; zmax = 169.355206

# File-specification inputs - Maintain the order as that of PACKMOL input
file_specs = [
    (head_dir+"/li_surface/lithium_large_24_17_6_edited.data", 1),
    (head_dir+"/peo_polymer/60PEO_Optimized_CH3terminated_editedterminal_edited.data", 20),
    (head_dir+"/li_monomer/Li_Atom_edited.data", 1000),
    (head_dir+"/co32m_monomer/co3_2minus_edited.data", 500),
    (head_dir+"/vecmtfsi_polymer/V30M2_150T_edited.data", 120),
    (head_dir+"/li_monomer/Li_Atom_edited.data",240)
]


# Main analysis
system = af.combine_lammps_system(file_specs, equil_pdb_file)
af.write_lammps_data(filename="combined_system.data",system=system,\
                     box=((xmin,xmax),(ymin,ymax),(zmin,zmax)),\
                     write_atoms='all_atoms.data')

print("Wrote combined_system.data")






