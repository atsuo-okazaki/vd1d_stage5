# usage:
#        gnuplot Mdisk.plt
#
set term pdfcairo transparent enhanced fontscale 1.0 size 5.00in, 4.00in 
INFILE = "r_edge_vs_t.dat"

### Labels
set label 1 left at graph 0.05, 0.92 "Opacity: OP+S03" font ",10"

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
#ymin_log = int(log10(STATS_min))
ymin_log = 0
ymax_log = int(log10(STATS_max)) + 1
ymin = 10.0**ymin_log
ymax = 10.0**ymax_log
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
unset ytics
set ytics add
do for [i = ymin_log:ymax_log] {
    set ytics add sprintf("10^{%d}", i) 10**i
}
set mytics 10
#unset logscale y
#unset format y
#ymin = 0
#ymax = 10
#set yrange [ymin:ymax]
set ylabel "Truncation radius [R_*]"

### legends
set key font ",10"
set key right top

OUTFILE = "r_edge.pdf"
set output OUTFILE
plot INFILE u 2:4 w l lw 3 lc "black" notitle

unset output
