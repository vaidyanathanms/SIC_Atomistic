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
set fname "../Copolymer_SingleChain/V27M5_150T"
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
    # See All_Charges.xlsx for the details
    # Select all atoms in each resID
    set selmain [atomselect $molID "resid $ref_resid"]
    set natoms [$selmain num]
    set ref_resname "BAD" ;# set default value
    puts "$natoms, $ref_resid, $ref_resname"
    
    # Set resname
    if { ($natoms == 15) && ($ref_resid == 1) } {  
	set ref_resname "VTR" ;# First VTR

	change_at_type_name_q	"O1"	"O1"	"OR"	"-0.4053"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C3"	"C3"	"CR"	"0.09452"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C4"	"C4"	"CR"	"0.31929"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C2"	"C2"	"C2"	"-0.16611"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C1"	"C1"	"C1"	"-0.28627"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O2"	"O2"	"OR"	"-0.44778"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C5"	"C5"	"CO"	"0.91462"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O3"	"O3"	"OC"	"-0.57795"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H7"	"H7"	"HR2"	"0.04915"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H6"	"H6"	"HR2"	"0.04915"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H5"	"H5"	"HR1"	"0.0386"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H4"	"H4"	"H2"	"0.06922"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H3"	"H3"	"H2"	"0.06922"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H2"	"H2"	"H1"	"0.13982"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H1"	"H1"	"H1"	"0.13982"	"$ref_resname"	"resid $ref_resid"	"top"
	
    } elseif { ($natoms == 15) && ($ref_resid == $max_resid) } {
	set ref_resname "VTR" ;# Last VTR

	change_at_type_name_q	"O1"	"O1"	"OR"	"-0.4053"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C3"	"C3"	"CR"	"0.09453"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C4"	"C4"	"CR"	"0.31929"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C2"	"C2"	"C2"	"-0.16611"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C1"	"C1"	"C1"	"-0.28627"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O2"	"O2"	"OR"	"-0.44778"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C5"	"C5"	"CO"	"0.91462"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O3"	"O3"	"OC"	"-0.57795"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H7"	"H7"	"HR2"	"0.04915"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H6"	"H6"	"HR2"	"0.04915"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H5"	"H5"	"HR1"	"0.0386"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H4"	"H4"	"H2"	"0.13844"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H3"	"H3"	"H1"	"0.09321"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H2"	"H2"	"H1"	"0.09321"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H1"	"H1"	"H1"	"0.09321"	"$ref_resname"	"resid $ref_resid"	"top"

    }  elseif { $natoms == 14 } {
	set ref_resname "VEC" ;# VEC
	
	change_at_type_name_q	"O1"	"O1"	"OR"	"-0.4053"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C3"	"C3"	"CR"	"0.09452"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C4"	"C4"	"CR"	"0.31929"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C2"	"C2"	"C2"	"-0.16611"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C1"	"C1"	"C1"	"-0.28627"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O2"	"O2"	"OR"	"-0.44778"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C5"	"C5"	"CO"	"0.91462"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O3"	"O3"	"OC"	"-0.57795"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H6"	"H6"	"HR2"	"0.04915"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H5"	"H5"	"HR2"	"0.04915"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H4"	"H4"	"HR1"	"0.0386"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H3"	"H3"	"H2"	"0.13844"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H2"	"H2"	"H1"	"0.13982"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H1"	"H1"	"H1"	"0.13982"	"$ref_resname"	"resid $ref_resid"	"top"

    }  elseif { $natoms == 31 } {
	set ref_resname "MT2" ;# VEC
	change_at_type_name_q	"H3"	"H3"	"HCT3"	"0.01929"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C3"	"C3"	"CT3"	"-0.11836"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C2"	"C2"	"CT2"	"0.07952"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H4"	"H4"	"HCT3"	"0.01929"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H5"	"H5"	"HCT3"	"0.01929"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C1"	"C1"	"CT3"	"-0.04352"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C4"	"C4"	"CMA"	"0.72866"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H00"	"H00"	"HCT2"	"0.00000"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H0"	"H0"	"HCT3"	"0.00000"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H1"	"H1"	"HCT3"	"0.0101"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H2"	"H2"	"HCT3"	"0.0101"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O1"	"O1"	"OM1"	"-0.43316"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O2"	"O2"	"OM2"	"-0.61245"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C5"	"C5"	"C2O"	"0.0934"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C6"	"C6"	"C2C"	"0.01849"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H6"	"H6"	"H1CO"	"0.0725"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H7"	"H7"	"H2CO"	"0.0725"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C7"	"C7"	"C2S"	"-0.18375"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H8"	"H8"	"HCT2"	"0.03809"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H9"	"H9"	"HCT2"	"0.03809"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"S1"	"S1"	"STF"	"1.24001"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H10"	"H10"	"H1CS"	"0.05658"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"H11"	"H11"	"H2CS"	"0.05658"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"N1"	"N1"	"NTF"	"-0.75545"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O3"	"O3"	"OS"	"-0.62063"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O4"	"O4"	"OS"	"-0.62063"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"S2"	"S2"	"STF"	"1.15213"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"C8"	"C8"	"CF"	"0.36973"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O5"	"O5"	"OS"	"-0.58437"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"O6"	"O6"	"OS"	"-0.58437"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"F1"	"F1"	"FC"	"-0.17922"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"F2"	"F2"	"FC"	"-0.17922"	"$ref_resname"	"resid $ref_resid"	"top"
	change_at_type_name_q	"F3"	"F3"	"FC"	"-0.17922"	"$ref_resname"	"resid $ref_resid"	"top"
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
