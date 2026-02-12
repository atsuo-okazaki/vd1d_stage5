# usage:
#        gnuplot -e 'IT="002970001"' scripts/profiles.plt
#
set term pdfcairo transparent enhanced fontscale 1.0 size 5.00in, 4.00in 
#INFILE = "profiles/disk_t".IT.".dat"
INFILE = "disk_t".IT.".dat"

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
m_col = 2          # r (or x) column used for x-axis normalization and plot x
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

m = m_col
stats INFILE every ::0::0 using m nooutput
x0 = STATS_mean

### x-axis
stats INFILE using m_col nooutput
xmin = 10.0**(int(log10(STATS_min/x0)))
xmax = 10.0**(int(log10(STATS_max/x0) + 0.99999))
set xrange [xmin:xmax]
set logscale x
set format x "10^{%L}"
set xlabel "r/r_{in}"


#################
# Plot 1: Sigma #
#################
n = c_sigma
stats INFILE using (column(n) > 0 ? column(n) : 1/0) nooutput
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

OUTFILE = "Sigma_t".IT.".pdf"
set output OUTFILE

plot INFILE using (column(m)/x0):(column(n)) \
     with lines lw 2 lc rgb "black"notitle

unset output
set yrange [*:*]
unset logscale y
unset format y

#################
# Plot 2: H/r   #
#################
# ratio = H/r and its derivative
### y-axis
n = c_y
#stats INFILE using ( (column(m) != 0 && column(n)/column(m) > 0) ? (column(n)/column(m)) : 1/0 ) nooutput
stats INFILE using (column(n) > 0 ? column(n) : 1/0) nooutput
ymin_log1 = floor(log10(STATS_min))
ymax_log1 = ceil(log10(STATS_max))

n = c_dydxi
stats INFILE using (column(n) > 0 ? column(n) : 1/0) nooutput
ymin_log2 = floor(log10(STATS_min))
ymax_log2 = ceil(log10(STATS_max))

if (ymin_log1 < ymin_log2) {
   ymin_log = ymin_log1
} else {
   ymin_log = ymin_log2
}
if (ymax_log1 > ymax_log2) {
   ymax_log = ymax_log1
} else {
   ymax_log = ymax_log2
}
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
set ylabel "Y(=H/r), dY/d{/Symbol x}"

### legends
set key font ",10"
set key left bottom

OUTFILE = "hor_t".IT.".pdf"
set output OUTFILE
plot INFILE using (column(m)/x0):(column(c_y)) \
     with lines lw 2 lc rgb "black" title "Y", \
	 INFILE using (column(m)/x0):(column(c_dydxi)) \
     with lines lw 2 lc rgb "blue" title "dY/d{/Symbol x}"

unset output
# --- Reset ticks to default for subsequent plots
#set ytics autofreq
#unset mytics

set yrange [*:*]
unset logscale y


#################
# Plot 3: rho   #
#################
### y-axis
n = c_rho
stats INFILE using (column(n) > 0 ? column(n) : 1/0) nooutput
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

OUTFILE = "rho_t".IT.".pdf"
set output OUTFILE
plot INFILE using (column(m)/x0):(column(n)) \
     with lines lw 2 lc rgb "black" notitle

unset output
set yrange [*:*]
unset logscale y

#################
# Plot 4: nu    #
#################
### y-axis
n = c_nu
stats INFILE using (column(n) > 0 ? column(n) : 1/0) nooutput
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

OUTFILE = "nu_t".IT.".pdf"
set output OUTFILE
plot INFILE using (column(m)/x0):(column(n)) \
     with lines lw 2 lc rgb "black" notitle

unset output
set yrange [*:*]
unset logscale y


##############################
# Plot 5: Qvis, Qirr, Qrad   #
##############################
### y-axis
n = c_qrad
stats INFILE using (column(n) > 0 ? column(n) : 1/0) nooutput
#ymin_log = floor(log10(STATS_min))
ymax_log = ceil(log10(STATS_max))

n = c_qvis
stats INFILE using (column(n) > 0 ? column(n) : 1/0) nooutput
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

OUTFILE = "E-balance_t".IT.".pdf"
set output OUTFILE

plot INFILE using (column(m)/x0):(column(c_qrad)) \
     with lines lw 2 lc rgb "black" title "Q_{rad}", \
	 INFILE using (column(m)/x0):(column(c_qvis)) \
     with lines lw 2 lc rgb "red" title "Q_{vis}", \
	 INFILE using (column(m)/x0):(column(c_qirr)) \
     with lines lw 2 lc rgb "blue" title "Q_{irr}"

unset output
set yrange [*:*]
unset logscale y


#########################
# Plot 6: Tmid vs Tsurf #
#########################
### y-axis
n1 = c_tmid
stats INFILE using (column(n1)  > 0 ? column(n1)  : 1/0) nooutput
tmin_log_9  = floor(log10(STATS_min))
tmax_log_9  = ceil(log10(STATS_max))

n2 = c_tsurf
stats INFILE using (column(n2) > 0 ? column(n2) : 1/0) nooutput
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

OUTFILE = "temp_t".IT.".pdf"
set output OUTFILE

plot INFILE using (column(m)/x0):(column(n1)) \
     with lines lw 2 lc rgb "black" title "T_{mid}", \
     INFILE using (column(m)/x0):(column(n2)) \
     with lines lw 2 lc rgb "blue" title "T_{surf}"

unset output
set yrange [*:*]
unset logscale y


##################
# Plot 7: kappaR #
##################
### y-axis
n = c_kappa
stats INFILE using (column(n) > 0 ? column(n) : 1/0) nooutput
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

OUTFILE = "kappaR_t".IT.".pdf"
set output OUTFILE

plot INFILE using (column(m)/x0):(column(n)) \
     with lines lw 2 lc rgb "black" notitle

unset output
set yrange [*:*]
unset logscale y


##################
# Plot 8: tauR   #
##################
### y-axis
n = c_tau
stats INFILE using (column(n) > 0 ? column(n) : 1/0) nooutput
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

OUTFILE = "tauR_t".IT.".pdf"
set output OUTFILE

plot INFILE using (column(m)/x0):(column(n)) \
     with lines lw 2 lc rgb "black" notitle

unset output
