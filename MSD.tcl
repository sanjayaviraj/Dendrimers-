######################################################
# find Diffusion coefficient of one type of selected #
# molecule within user defined radius.				 #
# 	sanjaya viraj bandara 2019-Dec					 #			
# 	sanjayavirajbandara@gmail.com					 #
######################################################

package require pbctools
#dont wrap when finding diffusion parameters

#set nf [molinfo top get numframes]
set outfile [open "MSD.dat" w]
#find center of mass of dendrimer(replace with your selection, reference mol)
set den_sel [atomselect top "segname D1" frame 0] 
set sel_center [measure center $den_sel weight mass]
set rx [lindex $sel_center 0]
set ry [lindex $sel_center 1]
set rz [lindex $sel_center 2]
#selection within 20 of center of mass
set sel_frame0 [atomselect top "(($rx-x)*($rx-x)+($ry-y)*($ry-y)+($rz-z)*($rz-z) < 400 ) and resname PHNL" frame 0] ; # replace Resname PHNL with atom type you want to find D
set atom_indexes [$sel_frame0 list]
set nf 50; #we are tracking entraped molecules; after 100 frames all initally traped moleculs dispers to water in worse case(z1);aire
set Max_tau 50
#bin is the time interval or Tau value
for {set bin 1} {$bin < $Max_tau} {incr bin} {
	
	set allMSD {}
	for {set f 0} {$f < $nf} {incr f} {
		#we wand to see how initialy traped molecules diffuse ,so selection should initialy traped molecules
		set sel1 [atomselect top "index $atom_indexes" frame $f]
		set sel2 [atomselect top "index $atom_indexes" frame [expr $f+$bin]]

		set sel1_cord [$sel1 get {x y z}]
		set sel2_cord [$sel2 get {x y z}]

		#cal separatly and append to a new list | can sub bcz x y z are in accending order of indexes**
		set veclength_list {}
		set j 0
		while {1} {
			set vec [veclength  [vecsub [lindex $sel1_cord $j] [lindex $sel2_cord $j]]]
			lappend veclength_list $vec
			incr j

			#puts "$vec"
			if {$j == [llength $sel1_cord]} {
				break
			}
		}	

		#get square of length vector
		set k 0
		set sq_veclength {}
		while {1} {
			set variable [lindex $veclength_list $k]
			set square_value [expr $variable*$variable]
			lappend sq_veclength $square_value
			#puts "$square_value"
			incr k
			if {$k==[llength $veclength_list]} {
				break
			}
		}

		set sq_veclength_sum [vecsum $sq_veclength]
		#puts "[llength $sq_veclength]"
		set nmol [$sel_frame0 num]
		#puts "$nmol"
		set MSD [expr $sq_veclength_sum / $nmol]
		lappend allMSD $MSD
		#puts "$MSD"

		#if {[expr $f%100]==0} {
		#	puts "frame: $f"
		#}
	}

	set sum_MSD [vecsum $allMSD]
	set av_MSD [expr $sum_MSD / ($nf+1) ]
	puts "$av_MSD"
	puts $outfile "[expr $bin*20] [expr $av_MSD/100]" ; #bin = Tau value in picosecounds , av_MSD in square nanometers
}
close $outfile
puts "Done"