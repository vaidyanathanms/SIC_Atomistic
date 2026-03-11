# To convert mol2 to PDB/PSF

package require topotools

# Set the input name without mol2 extension
set molname ../VEC_monomer/VEC_monomer_Optim_O2Changed

# Load the MOL2 file
mol new ${molname}.mol2

# Get the molecule ID
set molID [molinfo top]

# Select all atoms
set sel [atomselect $molID "all"]
set natoms [$sel num]

#----------------Main manipulations----------------------------------

# Example: set sel [atomselect $molID "index 0"] ; $sel set type "opls_135"
# Rename residue names. Will change atomtype and charges using change_atq.tcl


#---Residue name change
# Change identity of main group
set oxyethy [atomselect $molID "resname oxyethy"]
$oxyethy set resname "EO"

# Change identity of anion-group
#set C8H11NO [atomselect $molID "resname C8H11NO"]
#$C8H11NO set resname "ANI"


# Change identity of terminal-group-1
set methyl [atomselect $molID "resname methyl"]
$methyl set resname "CTR"

# Change identity of terminal-group-2
set hydroxy [atomselect $molID "resname hydroxy"]
$hydroxy set resname "OTR"


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
