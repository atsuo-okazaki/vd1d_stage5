# usage:
#        gnuplot sigma.plt
#
set term pdfcairo transparent enhanced fontscale 1.0 size 5.00in, 5.00in 
INFILE = "disk_evolution.dat"
OUTFILE = "disk_evolution.pdf"

# Get r0 in the 2nd
r0 = real(system( \
    "grep '^# r_min' ".INFILE." | sed 's/.*r0 *= *//'"))

stats INFILE using 0 nooutput
ncol = STATS_columns

### Labels
set label 1 right at graph 0.95, 0.93 "Opacity: OP+S03" font ",9"
#set label 2 right at graph 0.95, 0.85 "w/ wind-driven ablation" font ",9"
set label 3 left at graph 0.05, 0.06 sprintf("R_* = %.2e cm", r0) font ",9"

set output OUTFILE

### x-axis
stats INFILE using 1 nooutput
xmin = 10.0**(int(log10(STATS_min/r0)))
xmax = 10.0**(int(log10(STATS_max/r0)) + 1)
set xrange [xmin:xmax]
set logscale x
set format x "10^{%L}"
set mxtics 10
#set xlabel "r [cm]"
set xlabel "r/R_*"
#set xlabel sprintf("r/R_* (R_* = %.2e cm)", r0)

### y-axis
stats INFILE using ncol nooutput
ymin = 1e-9
ymax = 10.0**(int(log10(STATS_max)) + 2)
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set mytics 10
set ylabel "log {/Symbol S} (g/cm^2)"

### Labels
#set label 1 left at graph 0.05, 0.9 "Opacity: OP+S03" font ",10"

### legends
set key font ",8"
set key right top

plot for [c=1:ncol:5] INFILE u ($1/r0):c w l lw 1 lc "black" notitle #, \
#     ymax*(x/xmin)**(-2) w l lw 2 lc 'black'
#plot for [c=2:15] INFILE u 1:c w p pt 7 ps 0.1 notitle

unset output
