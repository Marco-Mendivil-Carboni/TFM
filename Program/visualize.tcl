#parameters
set res 32
set r_LJ 5.0

if {$argc==2} {
    #load initial condition
    package require topotools
    set sim_dir [lindex $argv 0]
    set sim_idx [lindex $argv 1]
    set gro_file [format "%s/initial-condition-%03d.gro" $sim_dir $sim_idx]
    color Display Background white
    display cuedensity 0.2
    display shadows on
    mol new $gro_file autobonds off
    set N [molinfo top get numatoms]
    for {set i 0} {$i<[expr $N-1]} {incr i} { topo addbond $i [expr $i+1]}
    set sel [atomselect top "name A"]
    $sel set radius $r_LJ
    color Name "A" 0
    set sel [atomselect top "name B"]
    $sel set radius $r_LJ
    color Name "B" 1
    mol modstyle 0 top CPK 4.0 [expr 4.0*$r_LJ/2.0] $res $res

    #draw sphere
    set param_file [format "%s/adjustable-parameters.dat" $sim_dir]
    set param_fp [open $param_file "r"]
    set param_lines [split [read $param_fp] "\n"]
    set R [expr 10.0*[scan [lindex $param_lines 1] "confinement_radius %f"]]
    draw material Transparent
    draw sphere {0 0 0} radius $R resolution $res

    #load trajectory files
    set pattern [format "%s/trajectory-%03d-*.trr" $sim_dir $sim_idx]
    set trr_files [lsort [glob $pattern]]
    foreach trr_file $trr_files { mol addfile $trr_file}
} else {
    puts "You forgot the input."
    exit
}
