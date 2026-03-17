# Change all necessary atomtypes/residues and charges

package require topotools

set molname "../inpcoord_files/lithium_large" ;# without the extension
set outname "${molname}_edited"

# Load the PDB file - Should have the necessary PDB files
mol new ${molname}.pdb type pdb

# Get the molecule ID
set molID [molinfo top]

# Select all atoms
set sel [atomselect $molID "all"]

#----------------Main manipulations----------------------------------

# Change identity of each group
# set oldresname [atomselect $molID "resname oldresname"]
# $oldresname set resname "newresname"

# change atomtype, charge, atomname
# Change Oxygens in COI
set sel2 [atomselect $molID "resname COI and name O00"]
set n [$sel2 num]
if {$n == 0} {
    puts "No atoms matched selection: $sel2"
    $sel2 delete
} else {
    set newResname "MOL"
    set newQ 0.0000
    set newType "Li"

    $sel2 set type      $newType
    $sel2 set charge    $newQ
    $sel2 set resname   $newResname
    $sel2 delete
}


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



# Optional: write out modified structure as MOL2 or PDB+PSF
# Requires TopoTools to write PSF
topo writelammpsdata ${outname}.data full
animate write psf ${outname}.psf
animate write pdb ${outname}.pdb

quit
