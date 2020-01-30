
#################################################
# Find Root mean square dispalcement/seperation #
#  sanjaya viraj bandara 2019-Dec		#
#  sanjayavirajbandara@gmail.com 		#
#################################################
package require pbctools
pbc wrap -centersel "segname D1" -center com -compound residue -all
set nf [molinfo top get numframes]
set output [open "RMSS.dat" w]

for {set f 500} {$f < $nf} {incr f} {
	
	set sel [atomselect top "index 5 6 9 8 28 29 25 26 40 41 11 12 19 20 44 43" frame $f]
	set cord [$sel get {x y z}]

	set distance_vec_list {}

	for {set i 0} {$i < [llength $cord]} {incr i} {

		for {set j 0} {$j < [llength $cord]} {incr j} {
			set distance_vec_sqr [ veclength2 [vecsub [lindex $cord $i] [lindex $cord $j]]]; #veclength2 returns the square of the scalar length of v
			if {[lsearch $distance_vec_list $distance_vec_sqr]==-1} {
				lappend distance_vec_list $distance_vec_sqr
				
			}
		}
		
	}

	#remove 0 that getting by substracting same indexes
	set in [lsearch $distance_vec_list 0.0]
	set distance_vec_list [lreplace $distance_vec_list $in $in]
	set sum_distance_vec_sqr [vecsum $distance_vec_list]
	set MSS [expr $sum_distance_vec_sqr/[llength $distance_vec_list]]
	set RMSS [expr { sqrt ($MSS)}]
	puts $output "$f $RMSS" ; # frame and RMMS in Angstrom
}

puts "DONE"
close $output


