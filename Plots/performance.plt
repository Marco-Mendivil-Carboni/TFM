reset

set terminal tikz color size 12cm,8cm
set output "performance.tex"

set style line 11 lw 2 lt 1 pt 7 lc rgb '#d81e2c'
set style line 12 lw 2 lt 1 pt 7 lc rgb '#a31cc5'
set style line 13 lw 2 lt 1 pt 7 lc rgb '#194bb2'
set style line 14 lw 2 lt 1 pt 7 lc rgb '#169f62'

set logscale x
set logscale y
set xrange [128:65536]
set yrange [0.1:100.0]
set xlabel "$N$"
set ylabel "$t_e$ (ms)"
set key bottom right
plot "GeForce-920M.dat" w linespoints t "GeForce-920M" ls 11,\
     "GeForce-RTX-3050.dat" w linespoints t "GeForce-RTX-3050" ls 12,\
     "RTX-A4000.dat" w linespoints t "RTX-A4000" ls 13,\

unset output
