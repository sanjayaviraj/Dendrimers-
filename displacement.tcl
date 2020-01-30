
#################################
# sanjaya viraj bandara 2019-Dec#
#  sanjayavirajbandara@gmail.com#
#################################
#find average displacement of phenol beads within r=20 of center of mass of dendrimer

package require pbctools
pbc wrap -centersel "segname D1" -center com -compound residue -all

set nf [molinfo top get numframes]
set outfile [open "displacements_within20ofDen_COM.dat" w]


for {set f 1} {$f < [expr $nf-1] } {incr f} {

#select all phenol beads in two frames close to eachother
	set sel_frame1 [atomselect top "resname PHNL" frame $f]
	set sel_frame2 [atomselect top "resname PHNL" frame [expr $f+1]]
	
	#get the cordinates (a list) from two frames seperately 
	set cordf1 [$sel_frame1 get {x y z}]
	set cordf2 [$sel_frame2 get {x y z}]
	
	#sub the x y z coordinates and get the veclength for a atom in two frames = distance traveled
	set Tvec_list {}
	set j 0
	
	#cal separatly and append to a new list | can sub bcz x y z are in accending order of indexes**
	while {1} {
		set Tvec [veclength  [vecsub [lindex $cordf1 $j] [lindex $cordf2 $j]]]
		lappend Tvec_list $Tvec
		incr j
		if {$j == [llength $cordf1]} {
			break
		}
	}

	#calculate the displacement of each atom and append to a new list
	set displacement_list {}
	
	#displacement A' per ps
	set time 20 
	foreach var $Tvec_list {
		set displacement [expr $var/$time]
		lappend displacement_list $displacement	
	}
	


	set den_sel [atomselect top "segname D1" frame [expr $f+1]] 
	set sel_center [measure center $den_sel weight mass]
	set rx [lindex $sel_center 0]
	set ry [lindex $sel_center 1]
	set rz [lindex $sel_center 2]
	#get atom indexes within r10 of CORE
	#call those indexes in displacement_arry and append to a new list r10_velo
	set r10 [atomselect top "(($rx-x)*($rx-x)+($ry-y)*($ry-y)+($rz-z)*($rz-z) < 400 ) and resname PHNL" frame [expr $f+1]]
	#if no beads selected ,skip this frame 
	set num_check [$r10 num]
		
	if {$num_check != 0 } {
		set sel_all [atomselect top all]
		set sel_all_list [$sel_all list]
		set r10indexlist [$r10 list]
		set allphnlindexlist [$sel_frame1 list]
		set matchindex_list {}
		#find the list indexes that matches the r10 selection from all selection 

		set k 0
		
		while {1} {			
			if {[lsearch $r10indexlist $k] != -1} {
				if {[lsearch $allphnlindexlist $k] !=-1} {
					set matchindex [lsearch $allphnlindexlist $k]
					lappend matchindex_list $matchindex			
				}
			}
			
			incr k
			if {$k == [llength $sel_all_list]} {
				break
			}
		}
			
		set r10_velo_list {}
		
		#get velocities from main list
		for {set i 0} {$i < [llength $matchindex_list]} {incr i} {
			set r10_velo [lindex $displacement_list	[lindex $matchindex_list $i]]
			lappend r10_velo_list $r10_velo
		}

		set a [vecsum $r10_velo_list]
		set r10_average_velo [expr $a/[llength $r10_velo_list]]

		puts $outfile "[expr $f+0.5] $r10_average_velo" ; # fame+0.5 , displacement in Angstrom per picosecound (A'/ps)
		
		if {[expr $f%20]==0} {
			puts "frame: $f"	
		}
	}	
}