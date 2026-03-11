# To convert mol2 to PDB/PSF and change the res names for MTFSI
# Rename residue names. Will change atomtype and charges using change_atq.tcl
# Example: set sel [atomselect $molID "index 0"] ; $sel set type "opls_135"
package require topotools

# Set the input name without mol2 extension
set molname ../MTFSI_Monomer/MTFSI_Monomer_O2Changed

# Load the MOL2 file
mol new ${molname}.mol2

# Get the molecule ID
set molID [molinfo top]

# Select all atoms
set sel [atomselect $molID "all"]
set natoms [$sel num]

#----------------Main manipulations----------------------------------

#---Residue name change
# Change identity of anion-group
set C8H13NO [atomselect $molID "resname C8H13NO"]
$C8H13NO set resname "ANI"

# Do topo operations to guess bonds/angles/dihedrals
topo numatomtypes
topo retypebonds
topo guessangles
topo guessdihedrals

# Set box sizes
molinfo $molID set a 100
molinfo $molID set b 100
molinfo $molID set c 100
molinfo $molID set alpha 90
molinfo $molID set beta  90
molinfo $molID set gamma 90


# Write PDB/PSF files
topo writelammpsdata ${molname}.data full
animate write psf ${molname}.psf
animate write pdb ${molname}.pdb
