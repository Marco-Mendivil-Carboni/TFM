set res 32
set r_LJ 5.0
if {$argc==2} {
    package require topotools
    set sim_dir [lindex $argv 0]
    set sim_idx [lindex $argv 1]
    set gro_file [format "%s/initial-configuration-%03d.gro" $sim_dir $sim_idx]
    set pattern [format "%s/trajectory-positions-%03d-*.trr" $sim_dir $sim_idx]
    set trr_files [lsort [glob $pattern]]
    color Display Background white
    display cuedensity 0.2
    mol new $gro_file autobonds off
    set N [molinfo top get numatoms]
    for {set i 0} {$i<[expr $N-1]} {incr i} { topo addbond $i [expr $i+1]}
    set sel [atomselect top "type X"]
    $sel set radius $r_LJ
    mol modstyle 0 top VDW 1.0 $res
    # mol modstyle 0 top licorice $r_LJ $res $res
    # mol modstyle 0 top CPK 4.0 [expr 4.0*$r_LJ/2.0] $res $res
    foreach trr_file $trr_files { mol addfile $trr_file}
    set param_file [format "%s/adjustable-parameters.dat" $sim_dir]
    set param_fp [open $param_file "r"]
    set param_lines [split [read $param_fp] "\n"]
    set R [expr 10.0*[scan [lindex $param_lines 2] "R\t%f"]]
    draw material Transparent
    draw sphere {0 0 0} radius $R resolution $res
} else {
    puts "You forgot the input."
    exit
}
