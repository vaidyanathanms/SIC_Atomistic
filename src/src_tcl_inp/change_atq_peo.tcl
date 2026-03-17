# Change all necessary atomtypes/residues and charges
# Version: Feb-16-2026, VMS
# See All_Charges.xlsx for the details

package require topotools

# Main procedure to change so that I can avoid writing the set variable 3 times
proc change_at_type_name_q {oldName newName newType {newCharge "0.0"} {newResname "BAD"} {scopeSel "all"} {molid top}} {
    set seltext "1"
    if {$scopeSel ne ""} {
	set seltext "$scopeSel"
    }
    
    puts "oldName $oldName newResname $newResname"

    # DO NOT FORGET THE '' around the atomtypes/names if there is a dot in the type or name
    set seltext "$seltext and name $oldName"
    set sel [atomselect $molid $seltext]

    set n [$sel num]

    if {$n == 0} {

	puts "No atoms matched selection: $seltext"
	$sel delete
	return

    }

    # Grab current properties for reporting
    set idxs   [$sel get index]
    set resn   [$sel get resname]
    set resid  [$sel get resid]
    set aname  [$sel get name]
    set atype  [$sel get type]
    set chgs   [$sel get charge]

    puts "Matched $n atoms for: $seltext"
    puts "Index  Resname Resid  Name   Type   Charge"

    for {set i 0} {$i < $n} {incr i} {
        puts [format "%5s  %6s %5s  %6s %6s %8.5f" \
            [lindex $idxs  $i] \
            [lindex $resn  $i] \
            [lindex $resid $i] \
            [lindex $aname $i] \
            [lindex $atype $i] \
            [lindex $chgs  $i]]
    }

    puts "Will set: type=$newType, name=$newName, charge=$newCharge, newResname=$resn"
    # Apply updates
    $sel set type      $newType
    $sel set name      $newName
    $sel set charge    $newCharge
    $sel set resname   $newResname

    $sel delete
    puts "Done."
}

#-----------------------------------Main function--------------------------------------
# Add the file prefix
set fname "../SEI_PEO_Polymer/60PEO_CH3terminated/60PEO_Optimized_CH3terminated_editedterminal"
set outname "${fname}_edited"

puts "${fname}.car"
# Load the car file
mol new ${fname}.car type car
mol addfile ${fname}.mdf type mdf

#----------------Main manipulations----------------------------------

# Get the molecule ID and the residue ID details
set molID [molinfo top]
set all_resids [lsort -unique [[atomselect top "all"] get resid]]
set max_resid [tcl::mathfunc::max {*}$all_resids]

# Loop through ResIDs
foreach ref_resid $all_resids {

    # NOTE: loop optimization by filtering common atoms did NOT work
    # So explicitly adding the definition for each atom
    # See https://pmc.ncbi.nlm.nih.gov/articles/PMC8037826/ for charges
    # Select all atoms in each resID
    set selmain [atomselect $molID "resid $ref_resid"]
    set natoms [$selmain num]
    set ref_aname [$selmain get name]
    set ref_resname "PEO" ;# set default value
    puts "$natoms, $ref_resid"

    # Set resname
    # Set charges according to https://pmc.ncbi.nlm.nih.gov/articles/PMC8037826/

    if { $ref_resid == $max_resid } {  
	set ref_resname "EOT" ;# Terminal residues from Materials Studio
	change_at_type_name_q	"H11"	"H11"	"HP"	"0.029"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H12"	"H12"	"HP"	"0.029"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H13"	"H13"	"HP"	"0.029"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"CT1"	"CT1"	"C1P"	"0.161" "$ref_resname"	"resid $ref_resid"	"top"
    } elseif { $ref_resid == $max_resid-1 } {  
	set ref_resname "EOT" ;# Terminal residues from Materials Studio
	change_at_type_name_q	"H21"	"H21"	"HP"	"0.029"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H22"	"H22"	"HP"	"0.029"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H23"	"H23"	"HP"	"0.029"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"CT2"	"CT2"	"C1P"	"0.161" "$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O"	"O"	"O1P"	"-0.472" "$ref_resname"	"resid $ref_resid"	"top"
    }  elseif { $ref_resid == 1 } {
	set ref_resname "EOI" ;
	# Note the charges of C2,H21 and H22 and O are changed since they are
	# part of "end monomer"
	change_at_type_name_q	"H11"	"H11"	"HP"	"-0.002" "$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H12"	"H12"	"HP"	"-0.002" "$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C1"	"C1"	"C1P"	"0.309"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H21"	"H21"	"HP"	"0.028" "EOT"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H22"	"H22"	"HP"	"0.028" "EOT"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C2"	"C2"	"C1P"	"0.168"	"EOT"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O"	"O"	"O1P"	"-0.472" "EOT"	"resid $ref_resid"	"top"
    } elseif { $ref_resid == $max_resid-2 } {
	set ref_resname "EOI" ;
	# Note the charges of C1,H11 and H12 and O are changed since they are
	# part of "end monomer"
	change_at_type_name_q	"H11"	"H11"	"HP"	"0.028" "EOT"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H12"	"H12"	"HP"	"0.028" "EOT"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C1"	"C1"	"C1P"	"0.168"	"EOT"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H21"	"H21"	"HP"	"-0.002" "$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H22"	"H22"	"HP"	"-0.002" "$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C2"	"C2"	"C1P"	"0.309"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O"	"O"	"O1P"	"-0.610" "$ref_resname"	"resid $ref_resid"	"top"
    }  else {
	set ref_resname "EOI" ; #Only for oxygen in resid=1
	change_at_type_name_q	"H11"	"H11"	"HP"	"-0.002" "$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H12"	"H12"	"HP"	"-0.002" "$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C1"	"C1"	"C1P"	"0.309"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H21"	"H21"	"HP"	"-0.002" "$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H22"	"H22"	"HP"	"-0.002" "$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C2"	"C2"	"C1P"	"0.309"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O"	"O"	"O1P"	"-0.610" "$ref_resname"	"resid $ref_resid"	"top"
    }  
}


#----------------Set atomnames-----------------------------------------
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

# Write out modified structure
# Requires TopoTools to write PSF
topo writelammpsdata ${outname}.data full
animate write psf ${outname}.psf
animate write pdb ${outname}.pdb


# Write top file
#mol load psf ${outname}.psf pdb ${outname}.pdb waitfor all
topo writegmxtop ${outname}.top 

quit
