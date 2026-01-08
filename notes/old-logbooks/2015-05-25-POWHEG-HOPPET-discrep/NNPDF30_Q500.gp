set term postscript color enhanced
set output 'NNPDF30_Q500.ps'

set grid
set title 'NNPDF30\_nnlo\_as\_0118, Q = 500 GeV'
plot 'NNPDF30_Q500.res' u 1:($3/$5) t 'dbar/u',\
     0 w l lt 1 lc 0 t ''

set output