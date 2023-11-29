if {$argc==2} {
    #set sim_dir and sim_idx
    set sim_dir [lindex $argv 0]
    set sim_idx [lindex $argv 1]

    #set visualization preferences
    color Display Background white
    display resize 640 960
    display cuedensity 0.2
    display shadows on
    display ambientocclusion on
    display aoambient 0.8
    display aodirect 0.8
    axes location Off
    set res 32
    set p_rad 5.0

    #load lamina binding sites
    set gro_file [format "%s/lamina-binding-sites-%03d.gro" $sim_dir $sim_idx]
    mol new $gro_file autobonds off
    set sel [atomselect top "name C"]
    $sel set radius $p_rad
    color Name "C" 2
    mol modstyle 0 top CPK 4.0 [expr 4.0*$p_rad/2.0] $res $res

    #load initial condition
    set gro_file [format "%s/initial-condition-%03d.gro" $sim_dir $sim_idx]
    mol new $gro_file autobonds off
    package require topotools
    set N [molinfo top get numatoms]
    for {set i 0} {$i<[expr $N-1]} {incr i} { topo addbond $i [expr $i+1]}
    set sel [atomselect top "name A"]
    $sel set radius $p_rad
    color Name "A" 0
    set sel [atomselect top "name B"]
    $sel set radius $p_rad
    color Name "B" 1
    mol modstyle 0 top CPK 4.0 [expr 4.0*$p_rad/2.0] $res $res

    #draw nucleus
    set param_file [format "%s/adjustable-parameters.dat" $sim_dir]
    set param_fp [open $param_file "r"]
    set param_lines [split [read $param_fp] "\n"]
    set R_n [expr 10.0*[scan [lindex $param_lines 1] "nucleus_radius %f"]]
    draw material Transparent
    draw sphere {0 0 0} radius $R_n resolution $res

    #change viewpoint
    scale by 0.8

    #load trajectory files
    set pattern [format "%s/trajectory-%03d-*.trr" $sim_dir $sim_idx]
    set trr_files [lsort [glob $pattern]]
    foreach trr_file $trr_files { mol addfile $trr_file}
} else {
    puts "You forgot the input."
    exit
}

# display resize 1200 1800
