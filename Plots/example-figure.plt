reset

set terminal tikz color size 12cm,8cm
set output "example-figure.tex"

set style line 01 lw 2 lt 1 lc rgb '#d37b81'
set style line 02 lw 2 lt 1 lc rgb '#b96bcd'
set style line 03 lw 2 lt 1 lc rgb '#5a7ec8'
set style line 04 lw 2 lt 1 lc rgb '#4ac28c'

set style line 11 lw 2 lt 1 lc rgb '#d81e2c'
set style line 12 lw 2 lt 1 lc rgb '#a31cc5'
set style line 13 lw 2 lt 1 lc rgb '#194bb2'
set style line 14 lw 2 lt 1 lc rgb '#169f62'

set style line 21 lw 2 pt 1 lc rgb '#9f000b'
set style line 22 lw 2 pt 1 lc rgb '#6a0084'
set style line 23 lw 2 pt 1 lc rgb '#002269'
set style line 24 lw 2 pt 1 lc rgb '#004f2b'

set xrange [0:4*pi]
set yrange [-1:+1]
set xlabel "$x$"
set ylabel "$J_n(x)$"
set key bottom right
set samples 64
plot "+" u ($1):(besjn(1,$1)):(0) w filledcurves t "" fc '#d37b81',\
     "+" u ($1):(besjn(2,$1)):(0) w filledcurves t "" fc '#b96bcd',\
     "+" u ($1):(besjn(3,$1)):(0) w filledcurves t "" fc '#5a7ec8',\
     "+" u ($1):(besjn(4,$1)):(0) w filledcurves t "" fc '#4ac28c',\
     besjn(1,x) t "n=1" ls 11,\
     besjn(2,x) t "n=2" ls 12,\
     besjn(3,x) t "n=3" ls 13,\
     besjn(4,x) t "n=4" ls 14,\
     "+" u ($1):(besjn(1,$1)) w points t "" ls 21,\
     "+" u ($1):(besjn(2,$1)) w points t "" ls 22,\
     "+" u ($1):(besjn(3,$1)) w points t "" ls 23,\
     "+" u ($1):(besjn(4,$1)) w points t "" ls 24

unset output
