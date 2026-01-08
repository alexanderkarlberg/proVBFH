# gnuplot file

set term postscript enhanced color
filename='apfel-hoppet-comp-mstw-new.ps'
set output filename

set log x
set xlabel 'x'
set xrange[1e-5:1]
set yrange [0.8:1.2]
unset multiplot
set lmargin at screen 0.1
set grid

#======================================================================
# F2
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NNLO(only) W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot 'apfel-F2Wpm.dat' u 1:(abs($6)) w l lw 2 t 'Apfel',\
     'hoppet-F2WpmZ.dat' u 1:(abs(($6))) w l lw 2 t 'hoppet',\
     'original-F2WpmZ.dat' u 1:(abs($6)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.9:1.1]
plot 1 w l lw 0.5 notitle,\
     '<paste apfel-F2Wpm.dat hoppet-F2WpmZ.dat' u 1:($13/($6)) w l t 'hoppet/Apfel',\
     '<paste apfel-F2Wpm.dat original-F2WpmZ.dat' u 1:($13/($6)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NNLO(only) W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot 'apfel-F2Wpm.dat' u 1:(abs($7)) w l lw 2 t 'Apfel',\
     'hoppet-F2WpmZ.dat' u 1:(abs(($7))) w l lw 2 t 'hoppet',\
     'original-F2WpmZ.dat' u 1:(abs($7)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.9:1.1]
plot 1 w l lw 0.5 notitle,\
     '<paste apfel-F2Wpm.dat hoppet-F2WpmZ.dat' u 1:($14/($7)) w l t 'hoppet/Apfel',\
     '<paste apfel-F2Wpm.dat original-F2WpmZ.dat' u 1:($14/($7)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#======================================================================
# F3
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, NNLO(only) W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot 'apfel-F3Wpm.dat' u 1:(abs($6)) w l lw 2 t 'Apfel',\
     'hoppet-F3WpmZ.dat' u 1:(abs(($6))) w l lw 2 t 'hoppet',\
     'original-F3WpmZ.dat' u 1:(abs($6)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.9:1.1]
plot 1 w l lw 0.5 notitle,\
     '<paste apfel-F3Wpm.dat hoppet-F3WpmZ.dat' u 1:($13/($6)) w l t 'hoppet/Apfel',\
     '<paste apfel-F3Wpm.dat original-F3WpmZ.dat' u 1:($13/($6)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, NNLO(only) W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot 'apfel-F3Wpm.dat' u 1:(abs($7)) w l lw 2 t 'Apfel',\
     'hoppet-F3WpmZ.dat' u 1:(abs(($7))) w l lw 2 t 'hoppet',\
     'original-F3WpmZ.dat' u 1:(abs($7)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.9:1.1]
plot 1 w l lw 0.5 notitle,\
     '<paste apfel-F3Wpm.dat hoppet-F3WpmZ.dat' u 1:($14/($7)) w l t 'hoppet/Apfel',\
     '<paste apfel-F3Wpm.dat original-F3WpmZ.dat' u 1:($14/($7)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

set output
