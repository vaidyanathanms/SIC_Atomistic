# SIC_Atomistic
Atomistic simulations of single-ion conductors

Assumptions
-----------
Impropers are not added (can be added later)
Have not added restraints to CO32- to mimic semi-solid
Li-surface motions are restricted to the x-y plane
NpT of bulk and crystal have not been done. Using the standard cutoff for now. This can be changed though.

Methodology
------------
1. Check how many lithium ions, CO_3^(2-) ions, VEC-MTFSI chains and the number of layers necessary.
2. Get the charge data are obtained from Gaussian in the Excel sheet.
3. Use Materials studio to generate a single molecule of each type and export as mol2 file, car and mdf files
4. Use src_tcl_inp/mol2topdb.tcl to convert to pdb files
5. Use src_tcl_inp/change_atq_vecmtfsi.tcl to change the resnames/atomtypes and create datafiles
5. Copy all needed PDB/datafiles to InputStructures/inp_coordfiles directory
6. Use PACKMOL to pack the necessary number of atoms/chains.
7. Use VMD/topotools and combine_lmp.py to generate the initial file. 
