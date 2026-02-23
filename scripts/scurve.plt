# usage:
#        gnuplot scurve.plt
#
set encoding utf8
set term pdfcairo transparent enhanced font "Helvetica,12" fontscale 1.0 size 5.00in, 5.00in
INFILE = "scurve.dat"
OUTFILE1 = "scurve1.pdf"
OUTFILE2 = "scurve2.pdf"

### x-axis
stats INFILE using 1 nooutput
xmin = 10.0**floor(log10(STATS_min))
xmax = 10.0**ceil(log10(STATS_max))
set xrange [xmin:xmax]
set logscale x
set format x "10^{%L}"
set mxtics 10
set xlabel "{/Symbol S}"

### y-axis
stats INFILE using 2 nooutput
ymin = 10.0**floor(log10(STATS_min))
ymax = 10.0**ceil(log10(STATS_max))
set yrange [ymin: ymax]
set logscale y
set format y "10^{%L}"
set mytics 10
set ylabel "Tmid"

# S-curve: Sigma vs T
set output OUTFILE1
plot "scurve.dat" u 1:2 with lines title "Sigma vs T"

unset output

# 不安定 branch のハイライト（stability 列が -1 の点）

x0 = 9.71684E+01 
y0 = 1.08438E+04

set output OUTFILE2
plot "scurve.dat" u 1:2:8 with lines notitle, \
     '+' u (x0):(y0) with points pt 7 notitle
#plot "scurve.dat" u 1:2:8 with lines, \
#     "scurve.dat" u ($8<0 ? $2:1/0):1 with points pt 7 title "unstable"
