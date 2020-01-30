#################################################
# Dakshitha Abegunawrdana and Sanjaya Viraj		#
#   13.12.2019									#
# script to select random phenol molecules      #
#################################################


set sel_all [atomselect top all]
 #give how many molecules you want to select
set goal 47
#give the selecton condition. within or all
set sel_dn [atomselect top "(same residue as (within 20 of resname CORE) or same residue as (within 7 of segname D1)) and resname PHNL" frame last]
set id_array [$sel_dn get index]
set idlength [$sel_all num]
#puts "$id_arry"

set random_list {};
set i 0;
while {1} {

	set random [expr {int(rand()*$idlength)}]
	
	#set random [expr { int (( expr {rand()} ) * ({llength $id_arry} ))} ]
	if {[lsearch $random_list $random]== -1 && [lsearch $id_array $random] != -1} {
		lappend random_list $random
		#puts "After append : $random"
		
		#find bonded atoms to this index 
		set sel_bonded [atomselect top "index $random"]
		set bonded_list [$sel_bonded getbonds]
		puts "$bonded_list"


		set vari [lindex $bonded_list 0 0]
		set sel_bonded_level2 [atomselect top "index $vari"]
		set bonded_list2 [$sel_bonded_level2 getbonds]

		#if your molecle contain only one selection skip this. if it is like a molecule with four beads use this
		#getbonds returns list of a lists (2tiered list) 
			lappend random_list [lindex $bonded_list 0 0]
			lappend random_list [lindex $bonded_list2 0 0]
			lappend random_list [lindex $bonded_list 0 1]
			lappend random_list [lindex $bonded_list2 0 1]
			puts "$i"

		incr i
	}
	
	if {$i == $goal} {

		break

	}

}

set final_phnl [atomselect top "index $random_list"]
$final_phnl writepdb random_slct.pdb

	

