reset

set terminal tikz color size 12cm,12cm
set output "analysis-summary.tex"

set style line 01 lw 2 lt 1 lc rgb '#d37b81'
set style line 02 lw 2 lt 1 lc rgb '#b96bcd'
set style line 03 lw 2 lt 1 lc rgb '#5a7ec8'
set style line 04 lw 2 lt 1 lc rgb '#4ac28c'

set style line 11 lw 2 lt 1 pt 7 lc rgb '#d81e2c'
set style line 12 lw 2 lt 1 pt 7 lc rgb '#a31cc5'
set style line 13 lw 2 lt 1 pt 7 lc rgb '#194bb2'
set style line 14 lw 2 lt 1 pt 7 lc rgb '#169f62'

set style fill solid 1.0

set multiplot layout 2,1

data = "../Simulations/example/analysis-fin.dat"
unset key

x_p = 0
set xlabel "$r$"
set ylabel "$\\rho(r)$"
set format y "$ %.1t\\cdot10^{%T} $"
plot data i 1 u 1:2:(x_p):(x_p=$1):($2-$4):($2+$4) w boxxyerror ls 01,\
     data i 1 w fsteps ls 11
set format y " %g "

set logscale x
set logscale y
set xlabel "$s$"
set ylabel "$d(s)$"
set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax
plot data i 2 u 1:($2-$4):($2+$4) w filledcurves ls 02,\
     data i 2 w lines ls 12

unset multiplot

unset output
