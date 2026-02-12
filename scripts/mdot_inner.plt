# usage:
#        gnuplot mdot_inner.plt
#
set encoding utf8
set term pdfcairo transparent enhanced font "Helvetica,12" fontscale 1.0 size 5.00in, 4.00in
#set term pdfcairo transparent enhanced fontscale 1.0 size 5.00in, 4.00in 
INFILE = "mdot_inner.dat"

### Labels
set label 1 left at graph 0.05, 0.92 "Opacity: {/Symbol k}_{es}+{/Symbol k}_{ff}" font ",10"

### x-axis
stats INFILE using 2 nooutput
xmin = 10.0**(int(log10(STATS_min)))
xmax = 10.0**(int(log10(STATS_max)) + 1)
set xrange [xmin:xmax]
set logscale x
set format x "10^{%L}"
set xlabel "t [s]"


### y-axis
stats INFILE using 4 nooutput
ymax = 10.0**(int(log10(STATS_max)) + 1)
ymin = ymax * 1e-4
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set ylabel "Ṁ(r_{in}) [g/s]"

### legends
set key font ",10"
set key right top

OUTFILE = "mdot_inner.pdf"
set output OUTFILE
plot INFILE u 2:4 w l lw 3 lc "black" notitle

unset output
