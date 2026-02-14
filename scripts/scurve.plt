# usage:
#        gnuplot scurve.plt
#
set encoding utf8
set term pdfcairo transparent enhanced font "Helvetica,12" fontscale 1.0 size 5.00in, 4.00in
INFILE = "scurve.dat"
OUTFILE1 = "scurve1.pdf"
OUTFILE2 = "scurve2.pdf"

# S-curve: Sigma vs T
set output OUTFILE1
plot "scurve.dat" u 2:1 with lines title "Sigma vs T"

unset output

# 不安定 branch のハイライト（stability 列が -1 の点）
set output OUTFILE2
plot "scurve.dat" u 2:1:8 with lines, \
     "scurve.dat" u ($8<0 ? $2:1/0):1 with points pt 7 title "unstable"
