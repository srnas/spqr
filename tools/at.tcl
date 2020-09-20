set N_ATS_PER_NT 5
set S_IND 3
set P_IND 4


set tempsel [atomselect top all]
set N_ATS [$tempsel num]

set N_NT [expr $N_ATS/$N_ATS_PER_NT]

for { set i 0 } {$i<[expr $N_NT -1] } {incr i} {
    set ch1 [lindex [$tempsel get chain] [expr $i*$N_ATS_PER_NT]]
    set ch2 [lindex [$tempsel get chain] [expr ($i+1)*$N_ATS_PER_NT]]
    set re1 [lindex [$tempsel get resid] [expr $i*$N_ATS_PER_NT]]
    set re2 [lindex [$tempsel get resid] [expr ($i+1)*$N_ATS_PER_NT]]
    
    if  { $ch1 == $ch2 && ( ( [expr $re1 - $re2] == 1 ) || ( [expr $re1 - $re2] == -1 ) ) }  {
	topo addbond [expr $i*$N_ATS_PER_NT + $S_IND] [expr ($i+1)*$N_ATS_PER_NT + $P_IND]
    }
    
}

#internally

for { set i 0 } {$i<[expr $N_NT] } {incr i} {
    #puts [expr $i*$N_ATS_PER_NT + $S_IND]
    topo addbond [expr $i*$N_ATS_PER_NT + $S_IND] [expr $i*$N_ATS_PER_NT ]
    topo addbond [expr $i*$N_ATS_PER_NT + $S_IND] [expr $i*$N_ATS_PER_NT + $P_IND]
    topo addbond [expr $i*$N_ATS_PER_NT + 1] [expr $i*$N_ATS_PER_NT ]
    topo addbond [expr $i*$N_ATS_PER_NT + 2] [expr $i*$N_ATS_PER_NT ]
    topo addbond [expr $i*$N_ATS_PER_NT + 1] [expr $i*$N_ATS_PER_NT +2]
}
mol modcolor 0 top Index
mol modstyle 0 top Licorice 0.300000 12.000000 12.000000
