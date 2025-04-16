reset
  
set grid
unset key
#set xrange [:]
#set yrange [:]
#set xtics -1.0, 0.2, 1.0
#set ytics -1.0, 0.2, 1.0
set size 0.6
#set size square
#set lmargin 5
#set rmargin 10
set xlabel "x label"
set ylabel "y label"
#set format x "%5.2f"
#set format y "%5.2f"
plot "test.dat" u 1:2

set terminal postscript color eps enhanced
set encoding iso_8859_1
set output "a.eps"
replot
set terminal X11
!eps2eps a.eps b.eps
!convert -density 300 b.eps Fig.tif
!rm a.eps b.eps

## MEMO
# lambda = "/Symbol l"
# angstrom = "{\305}"
