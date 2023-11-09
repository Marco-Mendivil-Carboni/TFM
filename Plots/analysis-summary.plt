reset

set terminal tikz standalone color size 26cm,26cm

set output "analysis-summary.tex"

set style line 01 lw 1 lt 1 lc rgb '#d37b81'
set style line 02 lw 1 lt 1 lc rgb '#b96bcd'
set style line 03 lw 1 lt 1 lc rgb '#5a7ec8'
set style line 04 lw 1 lt 1 lc rgb '#4ac28c'

set style line 11 lw 1 lt 1 pt 1 lc rgb '#d81e2c'
set style line 12 lw 1 lt 1 pt 1 lc rgb '#a31cc5'
set style line 13 lw 1 lt 1 pt 1 lc rgb '#194bb2'
set style line 14 lw 1 lt 1 pt 1 lc rgb '#169f62'

set style fill solid 1.0

N=16384
cvf=0.300
laf=0.200

N(n)=4096*2**n
cvf(n)=(n+1)*0.1
laf(n)=n*0.15+0.05

analysis_summary="../Simulations/analysis-summary.dat"

key(n)=sprintf("$\\varphi=%.3f$",cvf(n))
file(n)=sprintf("../Simulations/%5d-%.3f-%.3f/analysis-fin.dat",N,cvf(n),laf)

set multiplot layout 4,2 scale 0.98,0.98

set xrange [0.05:0.45]
set yrange [0.0:4.0]
set xlabel "$\\varphi$"
set ylabel "$ d_{cm} $"
set format y " %g "
plot analysis_summary i 0 u 2:($1==N&&$3==laf?$04:1/0):(1/0):6\
         w xyerrorbars ls 11 notitle
unset xrange
unset yrange

set xrange [0.05:0.45]
set yrange [100.0:600.0]
set xlabel "$\\varphi$"
set ylabel "$ R_g^2 $"
set format y " %g "
plot analysis_summary i 0 u 2:($1==N&&$3==laf?$07:1/0):(1/0):9\
         w xyerrorbars ls 11 notitle
unset xrange
unset yrange

set xrange [0.05:0.45]
set yrange [-0.1:0.0]
set xlabel "$\\varphi$"
set ylabel "$ S $"
set format y " %g "
plot analysis_summary i 0 u 2:($1==N&&$3==laf?$10:1/0):(1/0):12\
         w xyerrorbars ls 11 notitle
unset xrange
unset yrange

set xrange [0:]
set yrange [0:]
set xlabel "$r$"
set ylabel "$\\rho(r)$"
set format y "$ %.1t\\cdot10^{%T} $"
plot for [n=0:3] file(n) i 1+(x_p=0) u 1:2:(x_p):(x_p=$1):($2-$4):($2+$4)\
         w boxxyerror ls 1+n%4 notitle,\
     for [n=0:3] file(n) i 1 w fsteps ls 11+n%4 t key(n)
plot for [n=0:3] file(n) i 2+(x_p=0) u 1:2:(x_p):(x_p=$1):($2-$4):($2+$4)\
         w boxxyerror ls 1+n%4 notitle,\
     for [n=0:3] file(n) i 2 w fsteps ls 11+n%4 t key(n)
plot for [n=0:3] file(n) i 3+(x_p=0) u 1:2:(x_p):(x_p=$1):($2-$4):($2+$4)\
         w boxxyerror ls 1+n%4 notitle,\
     for [n=0:3] file(n) i 3 w fsteps ls 11+n%4 t key(n)
unset xrange
unset yrange

set key bottom right
set logscale x
set logscale y
set xlabel "$s$"
set ylabel "$d(s)$"
set format y " %g "
set xrange [1:8192]
set autoscale yfixmin
set autoscale yfixmax
plot for [n=0:3] file(n) i 4 u 1:($2-$4):($2+$4)\
         w filledcurves ls 1+n%4 notitle,\
     for [n=0:3] file(n) i 4 w lines ls 11+n%4 t key(n)

unset multiplot

unset output
