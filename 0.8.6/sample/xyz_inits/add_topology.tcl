
set N_ATS_PER_NT 5
set S_IND 3
set P_IND 4
set RCUT 10

set MCFILE [open "native.xyz" r]
set file_data [read $MCFILE]
set data [split $file_data "\n"]
#set N_NT [exec head -1 native.xyz]
set N_NT [expr [lindex $file_data 0]/$N_ATS_PER_NT]

puts [lindex $data 4]
for { set i 0 } {$i<[expr $N_NT -1] } {incr i} {
    #puts [expr $i*$N_ATS_PER_NT + $S_IND]
    set SUG1 [expr $i*$N_ATS_PER_NT+$S_IND+2]
    set SUG2 [expr ($i+1)*$N_ATS_PER_NT+$S_IND+2]
    #puts [lindex $data $SUG1]
    set x1 [lindex [lindex $data $SUG1] 1]
    set y1 [lindex [lindex $data $SUG1] 2]
    set z1 [lindex [lindex $data $SUG1] 3]

    set x2 [lindex [lindex $data $SUG2] 1]
    set y2 [lindex [lindex $data $SUG2] 2]
    set z2 [lindex [lindex $data $SUG2] 3]
    
    set DIST [expr sqrt(($x1-$x2)*($x1-$x2)+($y1-$y2)*($y1-$y2)+($z1-$z2)*($z1-$z2))]
    if  {$DIST < $RCUT } {
	topo addbond [expr $i*$N_ATS_PER_NT + $S_IND] [expr ($i+1)*$N_ATS_PER_NT + $P_IND]
    }
    
}


for { set i 0 } {$i<[expr $N_NT] } {incr i} {
    #puts [expr $i*$N_ATS_PER_NT + $S_IND]
    topo addbond [expr $i*$N_ATS_PER_NT + $S_IND] [expr $i*$N_ATS_PER_NT ]
    topo addbond [expr $i*$N_ATS_PER_NT + $S_IND] [expr $i*$N_ATS_PER_NT + $P_IND]
    topo addbond [expr $i*$N_ATS_PER_NT + 1] [expr $i*$N_ATS_PER_NT ]
    topo addbond [expr $i*$N_ATS_PER_NT + 2] [expr $i*$N_ATS_PER_NT ]
    topo addbond [expr $i*$N_ATS_PER_NT + 1] [expr $i*$N_ATS_PER_NT +2]
}
