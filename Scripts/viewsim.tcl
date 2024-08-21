#VMD visualization script

set ::pi 3.1415926535897931

proc draw_ring_quad {res z_a z_b r_a r_b a_a a_b a_c a_d} {
    set v_0 [list [expr $r_a*cos($a_c)] [expr $r_a*sin($a_c)] [expr $z_a]]
    set v_1 [list [expr $r_a*cos($a_d)] [expr $r_a*sin($a_d)] [expr $z_a]]
    set v_2 [list [expr $r_b*cos($a_c)] [expr $r_b*sin($a_c)] [expr $z_b]]
    set v_3 [list [expr $r_b*cos($a_d)] [expr $r_b*sin($a_d)] [expr $z_b]]
    set n_0 [list [expr cos($a_c)*cos($a_a)] [expr sin($a_c)*cos($a_a)] [expr sin($a_a)]]
    set n_1 [list [expr cos($a_d)*cos($a_a)] [expr sin($a_d)*cos($a_a)] [expr sin($a_a)]]
    set n_2 [list [expr cos($a_c)*cos($a_b)] [expr sin($a_c)*cos($a_b)] [expr sin($a_b)]]
    set n_3 [list [expr cos($a_d)*cos($a_b)] [expr sin($a_d)*cos($a_b)] [expr sin($a_b)]]
    graphics top trinorm $v_0 $v_1 $v_2 $n_0 $n_1 $n_2
    graphics top trinorm $v_3 $v_2 $v_1 $n_3 $n_2 $n_1
}

proc draw_ring {res z_a z_b r_a r_b a_a a_b} {
    for {set j 0} {$j < $res} {incr j} {
        set a_c [expr 2.0*$::pi*($j+0.0)/$res]
        set a_d [expr 2.0*$::pi*($j+1.0)/$res]
        draw_ring_quad $res $z_a $z_b $r_a $r_b $a_a $a_b $a_c $a_d
    }
}

proc draw_nucleus_without_bleb {res R_n} {
    for {set i 0} {$i < $res} {incr i} {
        set a_a [expr $::pi*($i+0.0)/$res-$::pi*0.5]
        set a_b [expr $::pi*($i+1.0)/$res-$::pi*0.5]
        set z_a [expr $R_n*sin($a_a)]
        set z_b [expr $R_n*sin($a_b)]
        set r_a [expr $R_n*cos($a_a)]
        set r_b [expr $R_n*cos($a_b)]
        draw_ring $res $z_a $z_b $r_a $r_b $a_a $a_b
    }
}

proc draw_nucleus_with_bleb {res R_n R_o R_b} {
    set d_b [expr sqrt($R_n*$R_n-$R_o*$R_o)+sqrt($R_b*$R_b-$R_o*$R_o)]
    set noa [expr 0.5*$::pi-asin($R_o/$R_n)]
    set boa [expr asin($R_o/$R_b)-0.5*$::pi]
    for {set i 0} {$i < $res} {incr i} {
        set a_a [expr ($noa+$::pi*0.5)*($i+0.0)/$res-$::pi*0.5]
        set a_b [expr ($noa+$::pi*0.5)*($i+1.0)/$res-$::pi*0.5]
        set z_a [expr $R_n*sin($a_a)]
        set z_b [expr $R_n*sin($a_b)]
        set r_a [expr $R_n*cos($a_a)]
        set r_b [expr $R_n*cos($a_b)]
        draw_ring $res $z_a $z_b $r_a $r_b $a_a $a_b
        set a_a [expr ($::pi*0.5-$boa)*($i+0.0)/$res+$boa]
        set a_b [expr ($::pi*0.5-$boa)*($i+1.0)/$res+$boa]
        set z_a [expr $d_b+$R_b*sin($a_a)]
        set z_b [expr $d_b+$R_b*sin($a_b)]
        set r_a [expr $R_b*cos($a_a)]
        set r_b [expr $R_b*cos($a_b)]
        draw_ring $res $z_a $z_b $r_a $r_b $a_a $a_b
    }
}

if {$argc>=2} {
    set sim_dir [lindex $argv 0]
    set sim_idx [lindex $argv 1]
    if {$argc>=3} {
        set t_f_idx [lindex $argv 2]
    }

    color Display Background white
    display resize 960 960
    display cuedensity 0.2
    display shadows on
    display ambientocclusion on
    display aoambient 0.8
    display aodirect 0.8
    axes location Off

    set res 32
    set r_p 5.0
    set n_l 0
    set R_o 0.0
    set param_file [format "%s/adjustable-parameters.dat" $sim_dir]
    set param_fp [open $param_file "r"]
    while {[gets $param_fp line] != -1} {
        if {[regexp "nucleus_radius\\s+(.*)" $line all value]} {
            set R_n [expr 10.0*$value]
        }
        if {[regexp "opening_radius\\s+(.*)" $line all value]} {
            set R_o [expr 10.0*$value]
        }
        if {[regexp "bleb_radius\\s+(.*)" $line all value]} {
            set R_b [expr 10.0*$value]
        }
        if {[regexp "number_of_lbs\\s+(.*)" $line all value]} {
            set n_l $value
        }
    }
    close $param_fp

    mol new
    graphics top material Transparent
    if {$R_o==0.0} {
        draw_nucleus_without_bleb $res $R_n
    } else {
        draw_nucleus_with_bleb $res $R_n $R_o $R_b
    }
    mol rename top {nucleus}

    set gro_file [format "%s/initial-condition-%03d.gro" $sim_dir $sim_idx]
    mol new $gro_file autobonds off
    package require topotools
    set N [molinfo top get numatoms]
    set chr_lim [list 0 2284 4599 6687 9338 12508 15193 17477 19792 21880 24531 27701 30386]
    set n_chr [expr [llength $chr_lim]-1]
    if {$N==[lindex $chr_lim $n_chr]} {
        for {set i_c 0} {$i_c<$n_chr} {incr i_c} {
            set i_s [lindex $chr_lim $i_c]
            set i_e [expr [lindex $chr_lim [expr $i_c+1]]-1]
            for {set i $i_s} {$i<$i_e} {incr i} { topo addbond $i [expr $i+1]}
        }
    } else {
        for {set i 0} {$i<[expr $N-1]} {incr i} { topo addbond $i [expr $i+1]}
    }
    set sel [atomselect top "name A"]
    $sel set radius $r_p
    color Name "A" 1
    set sel [atomselect top "name B"]
    $sel set radius $r_p
    color Name "B" 0
    mol modstyle 0 top CPK 4.0 [expr 4.0*$r_p/2.0] $res $res
    mol rename top {chromatin}

    if {$argc>=3} {
        set trr_file [format "%s/trajectory-%03d-%03d.trr" $sim_dir $sim_idx $t_f_idx]
        mol addfile $trr_file
    } else {
        set pattern [format "%s/trajectory-%03d-*.trr" $sim_dir $sim_idx]
        set trr_files [lsort [glob $pattern]]
        foreach trr_file $trr_files { mol addfile $trr_file}
    }

    set gro_file [format "%s/lamina-binding-sites-%03d.gro" $sim_dir $sim_idx]
    if {$n_l!=0} {
        mol new $gro_file autobonds off
        set sel [atomselect top "name C"]
        $sel set radius $r_p
        color Name "C" 2
        mol modstyle 0 top CPK 4.0 [expr 4.0*$r_p/2.0] $res $res
        mol rename top {lamina binding sites}
    }

    mol top 0
    display resetview

    if {$R_o!=0.0} {
        display resize 640 960
        rotate x to -90.0
    }
} else {
    puts "You forgot the input."
    exit
}
