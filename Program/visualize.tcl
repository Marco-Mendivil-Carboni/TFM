#parameters
set res 32
set p_rad 5.0

set ::pi 3.1415926535897931
set ::ex {1.0 0.0 0.0}
set ::ey {0.0 1.0 0.0}
set ::ez {0.0 0.0 1.0}

proc draw_ring { res z_ctr_a z_ctr_b r_ctr_a r_ctr_b theta_a theta_b} {
    for {set j 0} {$j < $res} {incr j} {
        set theta_c [expr 2.0*$::pi*($j+0.0)/$res]
        set theta_d [expr 2.0*$::pi*($j+1.0)/$res]
        set x_v_0 [vecscale $::ex [expr $r_ctr_a*cos($theta_c)]]
        set x_v_1 [vecscale $::ex [expr $r_ctr_a*cos($theta_d)]]
        set x_v_2 [vecscale $::ex [expr $r_ctr_b*cos($theta_c)]]
        set x_v_3 [vecscale $::ex [expr $r_ctr_b*cos($theta_d)]]
        set y_v_0 [vecscale $::ey [expr $r_ctr_a*sin($theta_c)]]
        set y_v_1 [vecscale $::ey [expr $r_ctr_a*sin($theta_d)]]
        set y_v_2 [vecscale $::ey [expr $r_ctr_b*sin($theta_c)]]
        set y_v_3 [vecscale $::ey [expr $r_ctr_b*sin($theta_d)]]
        set z_v_0 [vecscale $::ez [expr $z_ctr_a]]
        set z_v_1 [vecscale $::ez [expr $z_ctr_a]]
        set z_v_2 [vecscale $::ez [expr $z_ctr_b]]
        set z_v_3 [vecscale $::ez [expr $z_ctr_b]]
        set v_0 [vecadd $x_v_0 $y_v_0 $z_v_0]
        set v_1 [vecadd $x_v_1 $y_v_1 $z_v_1]
        set v_2 [vecadd $x_v_2 $y_v_2 $z_v_2]
        set v_3 [vecadd $x_v_3 $y_v_3 $z_v_3]
        set x_n_0 [vecscale $::ex [expr cos($theta_c)*cos($theta_a)]]
        set x_n_1 [vecscale $::ex [expr cos($theta_d)*cos($theta_a)]]
        set x_n_2 [vecscale $::ex [expr cos($theta_c)*cos($theta_b)]]
        set x_n_3 [vecscale $::ex [expr cos($theta_d)*cos($theta_b)]]
        set y_n_0 [vecscale $::ey [expr sin($theta_c)*cos($theta_a)]]
        set y_n_1 [vecscale $::ey [expr sin($theta_d)*cos($theta_a)]]
        set y_n_2 [vecscale $::ey [expr sin($theta_c)*cos($theta_b)]]
        set y_n_3 [vecscale $::ey [expr sin($theta_d)*cos($theta_b)]]
        set z_n_0 [vecscale $::ez [expr sin($theta_a)]]
        set z_n_1 [vecscale $::ez [expr sin($theta_a)]]
        set z_n_2 [vecscale $::ez [expr sin($theta_b)]]
        set z_n_3 [vecscale $::ez [expr sin($theta_b)]]
        set n_0 [vecadd $x_n_0 $y_n_0 $z_n_0]
        set n_1 [vecadd $x_n_1 $y_n_1 $z_n_1]
        set n_2 [vecadd $x_n_2 $y_n_2 $z_n_2]
        set n_3 [vecadd $x_n_3 $y_n_3 $z_n_3]
        draw trinorm $v_0 $v_1 $v_2 $n_0 $n_1 $n_2
        draw trinorm $v_3 $v_2 $v_1 $n_3 $n_2 $n_1
    }
}

if {$argc==2} {
    #set sim_dir and sim_idx
    set sim_dir [lindex $argv 0]
    set sim_idx [lindex $argv 1]

    #set visualization preferences
    color Display Background white
    display resize 720 1280
    display cuedensity 0.2
    display shadows on
    display ambientocclusion on
    display aoambient 0.8
    display aodirect 0.6
    axes location Off

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
    set R_o [expr 10.0*[scan [lindex $param_lines 2] "opening_radius %f"]]
    set R_b [expr 10.0*[scan [lindex $param_lines 3] "bleb_radius %f"]]
    set d_b [expr sqrt($R_n*$R_n-$R_o*$R_o)+sqrt($R_b*$R_b-$R_o*$R_o)]
    set noa [expr 0.5*$::pi-asin($R_o/$R_n)]
    set boa [expr asin($R_o/$R_b)-0.5*$::pi]
    draw material Transparent
    for {set i 0} {$i < $res} {incr i} {
        set theta_a [expr ($noa+$::pi*0.5)*($i+0.0)/$res-$::pi*0.5]
        set theta_b [expr ($noa+$::pi*0.5)*($i+1.0)/$res-$::pi*0.5]
        set z_ctr_a [expr $R_n*sin($theta_a)]
        set z_ctr_b [expr $R_n*sin($theta_b)]
        set r_ctr_a [expr $R_n*cos($theta_a)]
        set r_ctr_b [expr $R_n*cos($theta_b)]
        draw_ring $res $z_ctr_a $z_ctr_b $r_ctr_a $r_ctr_b $theta_a $theta_b
        set theta_a [expr ($::pi*0.5-$boa)*($i+0.0)/$res+$boa]
        set theta_b [expr ($::pi*0.5-$boa)*($i+1.0)/$res+$boa]
        set z_ctr_a [expr $d_b+$R_b*sin($theta_a)]
        set z_ctr_b [expr $d_b+$R_b*sin($theta_b)]
        set r_ctr_a [expr $R_b*cos($theta_a)]
        set r_ctr_b [expr $R_b*cos($theta_b)]
        draw_ring $res $z_ctr_a $z_ctr_b $r_ctr_a $r_ctr_b $theta_a $theta_b
    }

    #change viewpoint
    rotate x to 270.0
    translate to 0.0 -0.4 0.0
    scale by 0.8

    #load trajectory files
    set pattern [format "%s/trajectory-%03d-*.trr" $sim_dir $sim_idx]
    set trr_files [lsort [glob $pattern]]
    foreach trr_file $trr_files { mol addfile $trr_file}
} else {
    puts "You forgot the input."
    exit
}
