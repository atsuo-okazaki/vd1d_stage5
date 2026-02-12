# usage:
#        gnuplot Mdisk.plt
#
set term pdfcairo transparent enhanced fontscale 1.0 size 5.00in, 4.00in 
INFILE = "Mdisk_vs_t.dat"

### Labels
set label 1 left at graph 0.05, 0.92 "Opacity: OP+AES+S03" font ",9"

#stats INFILE every ::0::0 using 1 nooutput
#x0 = STATS_mean

### x-axis
stats INFILE using 2 nooutput
xmin = 10.0**(int(log10(STATS_min)))
xmax = 10.0**(int(log10(STATS_max)) + 1)
set xrange [xmin:xmax]
set logscale x
set format x "10^{%L}"
set xlabel "t [s]"


### y-axis
stats INFILE using 3 nooutput
ymin = 10.0**(int(log10(STATS_min)))
ymax = 10.0**(int(log10(STATS_max)) + 1)
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set ylabel "Disk mass [g]"

### legends
set key font ",10"
set key right top

OUTFILE = "Mdisk.pdf"
set output OUTFILE
plot INFILE u 2:3 w l lw 3 lc "black" notitle

unset output
