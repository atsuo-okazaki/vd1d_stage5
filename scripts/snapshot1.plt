# usage:
#        gnuplot -e 'IT="002970001"' scripts/profiles.plt
#
set term pdfcairo transparent enhanced fontscale 0.8 size 12.00in, 6.00in 
#INFILE = "profiles/disk_t".IT.".dat"
INFILE = "disk_t".IT.".dat"
OUTFILE = "disk_t".IT.".pdf"
set output OUTFILE

# Get time from INFILE
raw = system(sprintf("grep '^#' %s | head -1 | sed 's/.*t =  \\([0-9.E+-]*\\).*/\\1/'", INFILE))
raw_clean = sprintf("%s", raw)
t = real(raw_clean)
t_exp = floor(log10(t))
t_base = t / (10.0**t_exp)
t_fmt = sprintf("%.3g×10^{%d}", t_base, t_exp)

# -----------------------------
# Column definitions
# -----------------------------
c_r     = 2        # r (or x) column used for x-axis normalization and plot x
c_sigma = 3
c_h     = 4        # H
c_y     = 5        # y(=H/r)
c_dydxi = 6        # dy/dxi
c_rho   = 7
c_nu    = 8
c_qvis  = 9
c_qirr  = 10
c_qrad  = 11
c_tmid  = 12
c_tsurf = 13
c_kappa = 14
c_tau   = 15

### Labels
#set label 1 left at graph 0.05, 0.93 "Opacity: OP+S03" font ",10"
set label 2 right at graph 0.96, 0.93 "Time = ".t_fmt." s" font ",9"

stats INFILE every ::0::0 using c_r nooutput
x0 = STATS_mean

### x-axis
stats INFILE using c_r nooutput
xmin_log = floor(log10(STATS_min/x0))
xmax_log = ceil(log10(STATS_max/x0))
xmin = 10.0**xmin_log
xmax = 10.0**xmax_log
set xrange [xmin:xmax]
set logscale x
set format x "10^{%L}"
set xlabel "r/r_{in}"

set multiplot layout 2,4

#################
# Plot 1: H/r   #
#################
# ratio = H / r  (assuming col1=r, col3=H)
### y-axis
stats INFILE using ( (column(c_r) != 0 && column(c_h)/column(c_r) > 0) ? (column(c_h)/column(c_r)) : 1/0 ) nooutput

ymin_log = floor(log10(STATS_min))
ymax_log = ceil(log10(STATS_max))

# ymin_log = ymin_log - 1
# ymax_log = ymax_log + 0

ymin = 10.0**ymin_log
ymax = 10.0**ymax_log

set yrange [ymin:ymax]

set logscale y
set format y "10^{%L}"
#unset ytics
#set ytics add
#do for [i = ymin_log:ymax_log] {
#    set ytics add sprintf("10^{%d}", i) 10.0**i
#}
#set mytics 10
set ylabel "H/r"

### legends
set key font ",10"
set key right top

plot INFILE using (column(c_r)/x0):(column(c_h)/column(c_r)) \
     with lines lw 2 lc rgb "black" notitle

# --- Reset ticks to default for subsequent plots
#set ytics autofreq
#unset mytics

set yrange [*:*]
unset logscale y


#########################
# Plot 2: Tmid vs Tsurf #
#########################
### y-axis
stats INFILE using (column(c_tmid)  > 0 ? column(c_tmid)  : 1/0) nooutput
tmin_log_9  = floor(log10(STATS_min))
tmax_log_9  = ceil(log10(STATS_max))

stats INFILE using (column(c_tsurf) > 0 ? column(c_tsurf) : 1/0) nooutput
tmin_log_10 = floor(log10(STATS_min))
tmax_log_10 = ceil(log10(STATS_max))

ymin_log = (tmin_log_9  < tmin_log_10) ? tmin_log_9  : tmin_log_10
ymax_log = (tmax_log_9  > tmax_log_10) ? tmax_log_9  : tmax_log_10

ymin = 10.0**ymin_log
ymax = 10.0**ymax_log
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set ylabel "Temperature [K]"

### legends
set key font ",10"
set key left bottom

plot INFILE using (column(c_r)/x0):(column(c_tmid)) \
     with lines lw 2 lc rgb "black" title "T_{mid}", \
     INFILE using (column(c_r)/x0):(column(c_tsurf)) \
     with lines lw 2 lc rgb "blue" title "T_{surf}"

set yrange [*:*]
unset logscale y


#################
# Plot 3: rho   #
#################
### y-axis
stats INFILE using (column(c_rho) > 0 ? column(c_rho) : 1/0) nooutput
ymin_log = floor(log10(STATS_min))
ymax_log = ceil(log10(STATS_max))
ymin = 10.0**ymin_log
ymax = 10.0**ymax_log
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set ylabel "{/Symbol r} [g/cm^3]"

### legends
set key font ",10"
set key right top

plot INFILE using (column(c_r)/x0):(column(c_rho)) \
     with lines lw 2 lc rgb "black" notitle

set yrange [*:*]
unset logscale y


#################
# Plot 4: Sigma #
#################
stats INFILE using (column(c_sigma) > 0 ? column(c_sigma) : 1/0) nooutput
ymin_log = floor(log10(STATS_min))
ymax_log = ceil(log10(STATS_max))

# ymin_log = ymin_log - 1
# ymax_log = ymax_log + 1

ymin = 10.0**ymin_log
ymax = 10.0**ymax_log
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set mytics 10
set ylabel "{/Symbol S} [g/cm^2]"

### legends
set key font ",10"
set key right top

plot INFILE using (column(c_r)/x0):(column(c_sigma)) \
     with lines lw 2 lc rgb "black"notitle

set yrange [*:*]
unset logscale y
unset format y
unset mytics


#################
# Plot 5: nu    #
#################
### y-axis
stats INFILE using (column(c_nu) > 0 ? column(c_nu) : 1/0) nooutput
ymin_log = floor(log10(STATS_min))
ymax_log = ceil(log10(STATS_max))
ymin = 10.0**ymin_log
ymax = 10.0**ymax_log
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set ylabel "{/Symbol n} [cm^2/s]"

### legends
set key font ",10"
set key right top

plot INFILE using (column(c_r)/x0):(column(c_nu)) \
     with lines lw 2 lc rgb "black" notitle

set yrange [*:*]
unset logscale y


##################
# Plot 6: kappaR #
##################
### y-axis
stats INFILE using (column(c_kappa) > 0 ? column(c_kappa) : 1/0) nooutput
ymin_log = floor(log10(STATS_min) - 1)
ymax_log = ceil(log10(STATS_max) + 1)
ymin = 10.0**ymin_log
ymax = 10.0**ymax_log
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set ylabel "{/Symbol k} [cm^2/g]"

### legends
set key font ",10"
set key right top

plot INFILE using (column(c_r)/x0):(column(c_kappa)) \
     with lines lw 2 lc rgb "black" notitle

set yrange [*:*]
unset logscale y


##################
# Plot 7: tauR   #
##################
### y-axis
stats INFILE using (column(c_tau) > 0 ? column(c_tau) : 1/0) nooutput
ymin_log = floor(log10(STATS_min))
ymax_log = ceil(log10(STATS_max))
ymin = 10.0**ymin_log
ymax = 10.0**ymax_log
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set ylabel "{/Symbol t}"

### legends
set key font ",10"
set key right top

plot INFILE using (column(c_r)/x0):(column(c_tau)) \
     with lines lw 2 lc rgb "black" notitle

set yrange [*:*]
unset logscale y


##############################
# Plot 8: Qvis, Qirr, Qrad   #
##############################
### y-axis
stats INFILE using (column(c_qrad) > 0 ? column(c_qrad) : 1/0) nooutput
#ymin_log = floor(log10(STATS_min))
ymax_log = ceil(log10(STATS_max))

stats INFILE using (column(c_qvis) > 0 ? column(c_qvis) : 1/0) nooutput
ymin_log = floor(log10(STATS_min))
#ymax_log = ceil(log10(STATS_max))

ymin = 10.0**ymin_log
ymax = 10.0**ymax_log
set yrange [ymin:ymax]
set logscale y
set format y "10^{%L}"
set ylabel "Q_{vis}, Q_{irr}, Q_{rad} [erg/s/cm^2]"

### legends
set key font ",10"
set key left bottom

plot INFILE using (column(c_r)/x0):(column(c_qrad)) \
     with lines lw 2 lc rgb "black" title "Q_{rad}", \
	 INFILE using (column(c_r)/x0):(column(c_qvis)) \
     with lines lw 2 lc rgb "red" title "Q_{vis}", \
	 INFILE using (column(c_r)/x0):(column(c_qirr)) \
     with lines lw 2 lc rgb "blue" title "Q_{irr}"

unset multiplot
unset output
