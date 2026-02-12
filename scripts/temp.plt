# usage:
#        gnuplot sigma.plt
#
set term pdfcairo transparent enhanced fontscale 1.0 size 5.00in, 4.00in 

files = system("ls profiles/disk_t*.dat")
nfiles = words(files)

### x-axis
xmin = 1e10; xmax = 1e13
set xrange [xmin:xmax]
set logscale x
set format x "10^{%L}"
set xlabel "r (cm)"

### Labels
set label 1 left at graph 0.05, 0.9 "Opacity: {/Symbol k}_{es}+{/Symbol k}_{ff}" font ",10"


#################
# Plot 1: Tmid  #
#################
### y-axis
set logscale y
ymin = 1e2; ymax = 1e6
set yrange [ymin:ymax]
set format y "10^{%L}"
set ylabel "T_{mid} [K]"

### legends
set key font ",10"
set key right top

OUTFILE = "Tmid_evolution.pdf"
set output OUTFILE

plot for [i=1:nfiles:2] word(files, i) u 1:6 w l lc "black" notitle

unset output


#################
# Plot 1: Tsurf #
#################
### y-axis
set logscale y
ymin = 1e2; ymax = 1e6
set yrange [ymin:ymax]
set format y "10^{%L}"
set ylabel "T_{surf} [K]"

### legends
set key font ",10"
set key right top

OUTFILE = "Tsurf_evolution.pdf"
set output OUTFILE

plot for [i=1:nfiles:2] word(files, i) u 1:7 w l lc "blue" notitle

unset output
