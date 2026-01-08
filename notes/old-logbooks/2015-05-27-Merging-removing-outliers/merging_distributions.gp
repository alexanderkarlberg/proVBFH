set term postscript enhanced color 
set datafile fortran
set output 'merging_distribution.ps'

set xrange [1:101]
set ylabel 'picobarn'
set xlabel '%'
file = 'index4.top'
set title 'index 4'

plot   file i 0 u (200.*$5):3:4 w e lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):3:4 w e lc rgb 'red' lw 2 ti 'Sigma method'


reset

set xrange [-1:101]
set ylabel 'picobarn'
set xlabel '%'
file = 'index7.top'
set title 'index 7'

plot   file i 0 u (200.*$5):3:4 w e lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):3:4 w e lc rgb 'red' lw 2 ti 'Sigma method'


reset

set xrange [-1:101]
set ylabel 'picobarn'
set xlabel '%'
file = 'index8.top'
set title 'index 8'

plot   file i 0 u (200.*$5):3:4 w e lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):3:4 w e lc rgb 'red' lw 2 ti 'Sigma method'


reset

set xrange [-1:101]
set ylabel 'picobarn'
set xlabel '%'
file = 'index9.top'
set title 'index 9'

plot   file i 0 u (200.*$5):3:4 w e lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):3:4 w e lc rgb 'red' lw 2 ti 'Sigma method'


reset

set xrange [-1:101]
set yrange [-0.0002:0.0006]
set ylabel 'picobarn'
set xlabel '%'
file = 'index17.top'
set title 'index 17'

plot   file i 0 u (200.*$5):3:4 w e lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):3:4 w e lc rgb 'red' lw 2 ti 'Sigma method'


reset

set xrange [1:101]
set ylabel 'picobarn'
set xlabel '%'
file = 'index4.top'
set title 'index 4 error'

plot   file i 0 u (200.*$5):4 lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):4 lc rgb 'red' lw 2 ti 'Sigma method'


reset

set xrange [-1:101]
set ylabel 'picobarn'
set xlabel '%'
file = 'index7.top'
set title 'index 7 error'

plot   file i 0 u (200.*$5):4 lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):4 lc rgb 'red' lw 2 ti 'Sigma method'


reset

set xrange [-1:101]
set ylabel 'picobarn'
set xlabel '%'
file = 'index8.top'
set title 'index 8 error'

plot   file i 0 u (200.*$5):4 lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):4 lc rgb 'red' lw 2 ti 'Sigma method'


reset

set xrange [-1:101]
set ylabel 'picobarn'
set xlabel '%'
file = 'index9.top'
set title 'index 9 error'

plot   file i 0 u (200.*$5):4 lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):4 lc rgb 'red' lw 2 ti 'Sigma method'


reset

set xrange [-1:101]
set yrange [-0.0002:0.0006]
set ylabel 'picobarn'
set xlabel '%'
file = 'index17.top'
set title 'index 17 error'

plot   file i 0 u (200.*$5):4 lc rgb 'green' lw 2 ti 'Central method',\
       file i 1 u (100.*$6):4 lc rgb 'red' lw 2 ti 'Sigma method'


reset

set output