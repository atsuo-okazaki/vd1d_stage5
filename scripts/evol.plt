# ============================================================
# evol.plt
# usage: gnuplot scripts/evol.plt
# ============================================================

reset
set term pdfcairo transparent enhanced fontscale 1.0 size 5.00in, 4.00in
set encoding utf8

# -----------------------------
# Input files
# -----------------------------
files = system("ls disk_t*.dat")
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
c_h     = 4        # H in H/r
c_rho   = 7
c_nu    = 8
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
# Labels / annotations
# -----------------------------
set label 1 right at graph 0.96, 0.93 "Opacity: OP+AES+S03" font ",10"

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
set format x "10^{%L}"    # optional

# -----------------------------
# Helpers (manual conventions)
# -----------------------------
# Safety large/small initial values for global min/max
BIG  = 1e300
SMALL= 1e-300

# Fallback ranges (used if no valid positive data exist)
fallback_ymin_log = -3
fallback_ymax_log =  0

# Reset y-related state between plots (prevents tick carry-over)
# Use as a "macro block": just paste these lines after each plot.
#   set yrange [*:*]
#   unset logscale y
#   unset format y
#   set ytics autofreq
#   unset mytics

# ============================================================
# Function-like blocks: GLOBAL stats over all files
# (gnuplot does not have true functions with side effects,
#  so we use repeated patterns per plot for clarity/reliability.)
# ============================================================


# ============================================================
# Plot 1: Sigma
# ============================================================

# --- GLOBAL min/max excluding non-positive values
ymin_all = BIG
ymax_all = SMALL
do for [i=1:N] {
    f = word(files,i)
    stats f using (column(c_sigma) > 0 ? column(c_sigma) : 1/0) name "S" nooutput
    if (S_records - S_invalid > 0) {
       ymin_all = (S_min < ymin_all) ? S_min : ymin_all
       ymax_all = (S_max > ymax_all) ? S_max : ymax_all
    }
}

# --- Decide y-range for log scale
if (ymin_all >= BIG/10 || ymax_all <= SMALL*10 || ymin_all <= 0 || ymax_all <= 0 || ymin_all >= ymax_all) {
    print "WARNING: no valid positive data for this plot; using fallback range"
    #ymin_log = fallback_ymin_log
    #ymax_log = fallback_ymax_log
    ymin_log = -8.0
    ymax_log = 3.0
} else {
    ymin_log = floor(log10(ymin_all))
    ymax_log = ceil(log10(ymax_all))
}
set yrange [10.0**ymin_log : 10.0**ymax_log]
set logscale y
set format y "10^{%L}"
set mytics 10
set ylabel "{/Symbol S} [g/cm^2]"

set key font ",10"
set key right top

set output "evol_sigma.pdf"
plot for [i=1:N:INC] word(files,i) using (column(c_r)/x0[i]):(column(c_sigma)) \
     with lines lw 1 notitle
unset output

# --- reset y state
set yrange [*:*]
unset logscale y
unset format y
#set ytics autofreq
#unset mytics


# ============================================================
# Plot 2: H/r
# ============================================================

# --- GLOBAL min/max of ratio, excluding invalid/non-positive values
ymin_all = BIG
ymax_all = SMALL
do for [i=1:N] {
    f = word(files,i)
    stats f using ( (column(c_r) != 0 && column(c_h)/column(c_r) > 0) \
                    ? (column(c_h)/column(c_r)) : 1/0 ) name "S" nooutput
    if (S_records - S_invalid > 0) {
       ymin_all = (S_min < ymin_all) ? S_min : ymin_all
       ymax_all = (S_max > ymax_all) ? S_max : ymax_all
    }
}

if (ymin_all >= BIG/10 || ymax_all <= SMALL*10 || ymin_all <= 0 || ymax_all <= 0 || ymin_all >= ymax_all) {
    print "WARNING: no valid positive data for this plot; using fallback range"
    ymin_log = fallback_ymin_log
    ymax_log = fallback_ymax_log
} else {
    ymin_log = floor(log10(ymin_all))
    ymax_log = ceil(log10(ymax_all))
}
set yrange [10.0**ymin_log : 10.0**ymax_log]
set logscale y
set format y "10^{%L}"
set ylabel "H/r"

# --- Force major ticks at 10^n (optional but robust), plus minor ticks
#unset ytics
#set ytics add
#do for [k=ymin_log:ymax_log] {
#    set ytics add sprintf("10^{%d}", k) 10**k
#}
#set mytics 10

set key font ",10"
set key right top

set output "evol_hor.pdf"
plot for [i=1:N:INC] word(files,i) using (column(c_r)/x0[i]):(column(c_h)/column(c_r)) \
     with lines lw 1 notitle
unset output

# --- reset y state (CRITICAL: avoid tick carry-over)
set yrange [*:*]
unset logscale y
unset format y
#set ytics autofreq
#unset mytics


# ============================================================
# Plot 3: rho
# ============================================================

ymin_all = BIG
ymax_all = SMALL
do for [i=1:N] {
    f = word(files,i)
    stats f using (column(c_rho) > 0 ? column(c_rho) : 1/0) name "S" nooutput
    if (S_records - S_invalid > 0) {
       ymin_all = (S_min < ymin_all) ? S_min : ymin_all
       ymax_all = (S_max > ymax_all) ? S_max : ymax_all
    }
}

if (ymin_all >= BIG/10 || ymax_all <= SMALL*10 || ymin_all <= 0 || ymax_all <= 0 || ymin_all >= ymax_all) {
    print "WARNING: no valid positive data for this plot; using fallback range"
    #ymin_log = fallback_ymin_log
    #ymax_log = fallback_ymax_log
    ymin_log = -16.0
    ymax_log = 3.0
} else {
    ymin_log = floor(log10(ymin_all))
    ymax_log = ceil(log10(ymax_all))
}
set yrange [10.0**ymin_log : 10.0**ymax_log]
set logscale y
set format y "10^{%L}"
set ylabel "{/Symbol r} [g/cm^3]"

set key font ",10"
set key right top

set output "evol_rho.pdf"
plot for [i=1:N:INC] word(files,i) using (column(c_r)/x0[i]):(column(c_rho)) \
     with lines lw 1 notitle
unset output

set yrange [*:*]
unset logscale y
unset format y
#set ytics autofreq
#unset mytics


# ============================================================
# Plot 4: nu
# ============================================================

ymin_all = BIG
ymax_all = SMALL
do for [i=1:N] {
    f = word(files,i)
    stats f using (column(c_nu) > 0 ? column(c_nu) : 1/0) name "S" nooutput
    if (S_records - S_invalid > 0) {
       ymin_all = (S_min < ymin_all) ? S_min : ymin_all
       ymax_all = (S_max > ymax_all) ? S_max : ymax_all
    }
}

if (ymin_all >= BIG/10 || ymax_all <= SMALL*10 || ymin_all <= 0 || ymax_all <= 0 || ymin_all >= ymax_all) {
    print "WARNING: no valid positive data for this plot; using fallback range"
    #ymin_log = fallback_ymin_log
    #ymax_log = fallback_ymax_log
    ymin_log = 12.0
    ymax_log = 19.0
} else {
    ymin_log = floor(log10(ymin_all))
    ymax_log = ceil(log10(ymax_all))
}
set yrange [10.0**ymin_log : 10.0**ymax_log]
set logscale y
set format y "10^{%L}"
set ylabel "{/Symbol n} [cm^2/s]"

set key font ",10"
set key right top

set output "evol_nu.pdf"
plot for [i=1:N:INC] word(files,i) using (column(c_r)/x0[i]):(column(c_nu)) \
     with lines lw 1 notitle
unset output

set yrange [*:*]
unset logscale y
unset format y
#set ytics autofreq
#unset mytics


# ============================================================
# Plot 5: Tmid and Tsurf (shared y-range)
# ============================================================

# --- global range across BOTH columns and ALL files
ymin_all = BIG
ymax_all = SMALL

do for [i=1:N] {
    f = word(files,i)

    stats f using (column(c_tmid)  > 0 ? column(c_tmid)  : 1/0) name "S" nooutput
    if (S_records - S_invalid > 0) {
       ymin_all = (S_min < ymin_all) ? S_min : ymin_all
       ymax_all = (S_max > ymax_all) ? S_max : ymax_all
    }

    stats f using (column(c_tsurf) > 0 ? column(c_tsurf) : 1/0) name "S" nooutput
    if (S_records - S_invalid > 0) {
       ymin_all = (S_min < ymin_all) ? S_min : ymin_all
       ymax_all = (S_max > ymax_all) ? S_max : ymax_all
    }
}

if (ymin_all >= BIG/10 || ymax_all <= SMALL*10 || ymin_all <= 0 || ymax_all <= 0 || ymin_all >= ymax_all) {
    print "WARNING: no valid positive data for this plot; using fallback range"
    ymin_log = fallback_ymin_log
    ymax_log = fallback_ymax_log
} else {
    ymin_log = floor(log10(ymin_all))
    ymax_log = ceil(log10(ymax_all))
}
    ymin_log = 1.0
    ymax_log = 7.0

set yrange [10.0**ymin_log : 10.0**ymax_log]
set logscale y
set format y "10^{%L}"
set key font ",10"
set key left bottom

set output "evol_tmid.pdf"
set ylabel "T_{mid} [K]"
plot for [i=1:N:INC] word(files,i) using (column(c_r)/x0[i]):(column(c_tmid)) \
     with lines lw 1 notitle

set output "evol_tsurf.pdf"
set ylabel "T_{surf} [K]"
plot for [i=1:N:INC] word(files,i) using (column(c_r)/x0[i]):(column(c_tsurf)) \
     with lines lw 1 notitle

unset output

set yrange [*:*]
unset logscale y
unset format y
#set ytics autofreq
#unset mytics


# ============================================================
# Plot 6: kappaR
# ============================================================

ymin_all = BIG
ymax_all = SMALL
do for [i=1:N] {
    f = word(files,i)
    stats f using (column(c_kappa) > 0 ? column(c_kappa) : 1/0) name "S" nooutput
    if (S_records - S_invalid > 0) {
       ymin_all = (S_min < ymin_all) ? S_min : ymin_all
       ymax_all = (S_max > ymax_all) ? S_max : ymax_all
    }
}

if (ymin_all >= BIG/10 || ymax_all <= SMALL*10 || ymin_all <= 0 || ymax_all <= 0 || ymin_all >= ymax_all) {
    print "WARNING: no valid positive data for this plot; using fallback range"
    ymin_log = fallback_ymin_log
    ymax_log = fallback_ymax_log
} else {
    ymin_log = floor(log10(ymin_all))
    ymax_log = ceil(log10(ymax_all))
}
    #ymin_log = -4.0
    #ymax_log = 3.0
set yrange [10.0**ymin_log : 10.0**ymax_log]
set logscale y
set format y "10^{%L}"
set ylabel "{/Symbol k}_{R} [cm^2/g]"

set key font ",10"
set key right top

set output "evol_kappaR.pdf"
plot for [i=1:N:INC] word(files,i) using (column(c_r)/x0[i]):(column(c_kappa)) \
     with lines lw 1 notitle
unset output

set yrange [*:*]
unset logscale y
unset format y
#set ytics autofreq
#unset mytics


# ============================================================
# Plot 7: tauR
# ============================================================

ymin_all = BIG
ymax_all = SMALL
do for [i=1:N] {
    f = word(files,i)
    stats f using (column(c_tau) > 0 ? column(c_tau) : 1/0) name "S" nooutput
    if (S_records - S_invalid > 0) {
       ymin_all = (S_min < ymin_all) ? S_min : ymin_all
       ymax_all = (S_max > ymax_all) ? S_max : ymax_all
    }
}

if (ymin_all >= BIG/10 || ymax_all <= SMALL*10 || ymin_all <= 0 || ymax_all <= 0 || ymin_all >= ymax_all) {
    print "WARNING: no valid positive data for this plot; using fallback range"
    ymin_log = fallback_ymin_log
    ymax_log = fallback_ymax_log
} else {
    ymin_log = floor(log10(ymin_all))
    ymax_log = ceil(log10(ymax_all))
}
    #ymin_log = -4.0
    #ymax_log = 3.0
set yrange [10.0**ymin_log : 10.0**ymax_log]
set logscale y
set format y "10^{%L}"
set ylabel "{/Symbol t}_{R}"

set key font ",10"
set key right top

set output "evol_tauR.pdf"
plot for [i=1:N:INC] word(files,i) using (column(c_r)/x0[i]):(column(c_tau)) \
     with lines lw 1 notitle
unset output

# Final reset (optional)
set yrange [*:*]
unset logscale y
unset format y
#set ytics autofreq
#unset mytics

print "Done: generated evol_*.pdf"
# ============================================================
