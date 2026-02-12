# usage:
#        gnuplot scripts/snapshots.plt
#
set term pdfcairo transparent enhanced fontscale 0.8 size 12.00in, 6.00in 
set encoding utf8

# -----------------------------
# Input files
# -----------------------------
files = system("ls disk_t????????.dat")
N = words(files)
if (N <= 0) {
    print "ERROR: no input files matched disk_t*.dat"
    exit
}

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

# -----------------------------
# Plotting frequenct
# -----------------------------
INC = 1

# -----------------------------
# X-axis normalization: x0[i] = inner radius in each file
# (You used STATS_mean of the first row; keep as-is)
# -----------------------------
array x0[N]
do for [i=1:N] {
    f = word(files, i)
    stats f every ::0::0 using c_r name "X" nooutput
    x0[i] = X_mean
}

# -----------------------------
# X-axis setup
# -----------------------------
stats word(files, 1) using c_r nooutput
xmin_log = floor(log10(STATS_min/x0[1]))
xmax_log = ceil(log10(STATS_max/x0[1]))
xmin = 10.0**xmin_log
xmax = 10.0**xmax_log
set xrange [xmin:xmax]
set logscale x
set mxtics 10
set xlabel "r/r_{in}"
set format x "10^{%L}"

### legends
set key font ",10"
set key left bottom

do for [i=1:N] {
    f = word(files, i)
    # Get time from INFILE
    raw = system(sprintf("grep '^#' %s | head -1 | sed 's/.*t =  \\([0-9.E+-]*\\).*/\\1/'", f))
    raw_clean = sprintf("%s", raw)
    t = real(raw_clean)
    t_exp = floor(log10(t))
    t_base = t / (10.0**t_exp)
    t_fmt = sprintf("%.3g×10^{%d}", t_base, t_exp)
    set label right at graph 0.96, 0.93 "Time = ".t_fmt." s" font ",9"

	set output "disk_".substr(f, 6, 14).".pdf"

    set multiplot layout 2,4

    #################
    # Plot 1: H/r   #
    #################
    ymin = 1.0e-4
    ymax = 1.0e-1
    set yrange [ymin:ymax]
    set logscale y
    set format y "10^{%L}"
    set ylabel "H/r"

    plot f using (column(c_r)/x0[1]):(column(c_h)/column(c_r)) \
         with lines lw 2 lc rgb "black" notitle

    set yrange [*:*]
    unset logscale y

    #########################
    # Plot 2: Tmid vs Tsurf #
    #########################
    ymin = 1.0e1
    ymax = 1.0e7
    set yrange [ymin:ymax]
    set logscale y
    set format y "10^{%L}"
    set mytics 10
    set ylabel "Temperature [K]"

    plot f using (column(c_r)/x0[1]):(column(c_tmid)) \
         with lines lw 2 lc rgb "black" title "T_{mid}", \
         f using (column(c_r)/x0[1]):(column(c_tsurf)) \
         with lines lw 2 lc rgb "blue" title "T_{surf}"

    set yrange [*:*]
    unset logscale y


    #################
    # Plot 3: rho   #
    #################
    ymin = 1e-12
    ymax = 1e-5
    set yrange [ymin:ymax]
    set logscale y
    set format y "10^{%L}"
    set ylabel "{/Symbol r} [g/cm^3]"

    plot f using (column(c_r)/x0[1]):(column(c_rho)) \
         with lines lw 2 lc rgb "black" notitle

    set yrange [*:*]
    unset logscale y


    #################
    # Plot 4: Sigma #
    #################
    ymin = 1.0e-2
    ymax = 1.0e4
    set yrange [ymin:ymax]
    set logscale y
    set format y "10^{%L}"
    set mytics 10
    set ylabel "{/Symbol S} [g/cm^2]"

    plot f using (column(c_r)/x0[1]):(column(c_sigma)) \
         with lines lw 2 lc rgb "black"notitle

    set yrange [*:*]
    unset logscale y
    unset format y
    unset mytics


    #################
    # Plot 5: nu    #
    #################
    ymin = 1e12
    ymax = 1e17
    set yrange [ymin:ymax]
    set logscale y
    set format y "10^{%L}"
    set mytics 10
    set ylabel "{/Symbol n} [cm^2/s]"

    plot f using (column(c_r)/x0[1]):(column(c_nu)) \
         with lines lw 2 lc rgb "black" notitle

    set yrange [*:*]
    unset logscale y
    unset mytics

    ##################
    # Plot 6: kappaR #
    ##################
    ymin = 1e-4
    ymax = 1e2
    set yrange [ymin:ymax]
    set logscale y
    set format y "10^{%L}"
    set mytics 10
    set ylabel "{/Symbol k} [cm^2/g]"

    plot f using (column(c_r)/x0[1]):(column(c_kappa)) \
         with lines lw 2 lc rgb "black" notitle

    set yrange [*:*]
    unset logscale y
    unset mytics

    ##################
    # Plot 7: tauR   #
    ##################
    ymin = 1.0e-2
    ymax = 1e4
    set yrange [ymin:ymax]
    set logscale y
    set format y "10^{%L}"
    set mytics 10
    set ylabel "{/Symbol t}"

    plot f using (column(c_r)/x0[1]):(column(c_tau)) \
         with lines lw 2 lc rgb "black" notitle

    set yrange [*:*]
    unset logscale y
    unset mytics


    ##############################
    # Plot 8: Qvis, Qirr, Qrad   #
    ##############################
    ymin = 1.0e4
    ymax = 1.0e14
    set yrange [ymin:ymax]
    set logscale y
    set format y "10^{%L}"
    set mytics 10
    set ylabel "Q_{vis}, Q_{irr}, Q_{rad} [erg/s/cm^2]"

    plot f using (column(c_r)/x0[1]):(column(c_qrad)) \
         with lines lw 2 lc rgb "black" title "Q_{rad}", \
    	 f using (column(c_r)/x0[1]):(column(c_qvis)) \
         with lines lw 2 lc rgb "red" title "Q_{vis}", \
    	 f using (column(c_r)/x0[1]):(column(c_qirr)) \
         with lines lw 2 lc rgb "blue" title "Q_{irr}"

    unset multiplot
    unset label
    unset output
}
