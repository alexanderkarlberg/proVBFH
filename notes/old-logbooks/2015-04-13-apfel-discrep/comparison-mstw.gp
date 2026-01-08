# gnuplot file

set term postscript enhanced color
filename='apfel-hoppet-comp-mstw.ps'
set output filename

set log x
set xlabel 'x'
set xrange[1e-5:1]

unset multiplot
set lmargin at screen 0.1
set grid

#======================================================================
# F2
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, LO W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot 'MSTW_F2e-p_LO_Q100.dat' u 2:(abs($3)) w l lw 2 t 'Apfel',\
     'hoppet-F2WpmZ.dat' u 1:(abs($2)) w l lw 2 t 'hoppet',\
     'original-F2WpmZ.dat' u 1:(abs($2)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.999:1.001]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F2e-p_LO_Q100.dat hoppet-F2WpmZ.dat' u 2:($5/$3) w l t 'hoppet/Apfel',\
     '<paste MSTW_F2e-p_LO_Q100.dat original-F2WpmZ.dat' u 2:($5/$3) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, LO W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot 'MSTW_F2e+p_LO_Q100.dat' u 2:(abs($3)) w l t 'Apfel',\
     'hoppet-F2WpmZ.dat' u 1:(abs($3)) w l lw 2 t 'hoppet',\
     'original-F2WpmZ.dat' u 1:(abs($3)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.999:1.001]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F2e+p_LO_Q100.dat hoppet-F2WpmZ.dat' u 2:($6/$3) w l t 'hoppet/Apfel',\
     '<paste MSTW_F2e+p_LO_Q100.dat original-F2WpmZ.dat' u 2:($6/$3) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NLO(only) W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste MSTW_F2e-p_LO_Q100.dat MSTW_F2e-p_NLO_Q100.dat' u 2:(abs($6-$3)) w l t 'Apfel',\
     'hoppet-F2WpmZ.dat' u 1:(abs(($4))) w l lw 2 t 'hoppet',\
     'original-F2WpmZ.dat' u 1:(abs($4)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.98:1.02]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F2e-p_LO_Q100.dat MSTW_F2e-p_NLO_Q100.dat hoppet-F2WpmZ.dat' u 2:(($10)/($6-$3)) w l t 'hoppet/Apfel',\
     '<paste MSTW_F2e-p_LO_Q100.dat MSTW_F2e-p_NLO_Q100.dat original-F2WpmZ.dat' u 2:(($10)/($6-$3)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NLO(only) W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste MSTW_F2e+p_LO_Q100.dat MSTW_F2e+p_NLO_Q100.dat' u 2:(abs(($6-$3))) w l t 'Apfel',\
     'hoppet-F2WpmZ.dat' u 1:(abs(($5))) w l lw 2 t 'hoppet',\
     'original-F2WpmZ.dat' u 1:(abs($5)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.98:1.02]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F2e+p_LO_Q100.dat MSTW_F2e+p_NLO_Q100.dat hoppet-F2WpmZ.dat' u 2:(($11)/($6-$3)) w l t 'hoppet/Apfel',\
     '<paste MSTW_F2e+p_LO_Q100.dat MSTW_F2e+p_NLO_Q100.dat original-F2WpmZ.dat' u 2:(($11)/($6-$3)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NNLO(only) W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste MSTW_F2e-p_NLO_Q100.dat MSTW_F2e-p_NNLO_Q100.dat' u 2:(abs(($6-$3))) w l t 'Apfel',\
     'hoppet-F2WpmZ.dat' u 1:(abs(($6))) w l lw 2 t 'hoppet',\
     'original-F2WpmZ.dat' u 1:(abs($6)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.9:1.1]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F2e-p_NLO_Q100.dat MSTW_F2e-p_NNLO_Q100.dat hoppet-F2WpmZ.dat' u 2:(($12)/($6-$3)) w l t 'hoppet/Apfel',\
     '<paste MSTW_F2e-p_NLO_Q100.dat MSTW_F2e-p_NNLO_Q100.dat original-F2WpmZ.dat' u 2:(($12)/($6-$3)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NNLO(only) W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste MSTW_F2e+p_NLO_Q100.dat MSTW_F2e+p_NNLO_Q100.dat' u 2:(abs(($6-$3))) w l t 'Apfel',\
     'hoppet-F2WpmZ.dat' u 1:(abs(($7))) w l lw 2 t 'hoppet',\
     'original-F2WpmZ.dat' u 1:(abs($7)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.:1.2]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F2e+p_NLO_Q100.dat MSTW_F2e+p_NNLO_Q100.dat hoppet-F2WpmZ.dat' u 2:(($13)/($6-$3)) w l t 'hoppet/Apfel',\
     '<paste MSTW_F2e+p_NLO_Q100.dat MSTW_F2e+p_NNLO_Q100.dat original-F2WpmZ.dat' u 2:(($13)/($6-$3)) w l dt 2 t 'Zaro/Apfel'
unset multiplot


#======================================================================
# F3
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, LO W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot 'MSTW_F3e-p_LO_Q100.dat' u 2:(abs($3)) w l lw 2 t 'Apfel',\
     'hoppet-F3WpmZ.dat' u 1:(abs($2)) w l lw 2 t 'hoppet',\
     'original-F3WpmZ.dat' u 1:(abs($2)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.998:1.002]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F3e-p_LO_Q100.dat hoppet-F3WpmZ.dat' u 2:($5/$3) w l t 'hoppet/Apfel',\
     '<paste MSTW_F3e-p_LO_Q100.dat original-F3WpmZ.dat' u 2:($5/$3) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, LO W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot 'MSTW_F3e+p_LO_Q100.dat' u 2:(abs($3)) w l t 'Apfel',\
     'hoppet-F3WpmZ.dat' u 1:(abs($3)) w l lw 2 t 'hoppet',\
     'original-F3WpmZ.dat' u 1:(abs($3)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.998:1.002]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F3e+p_LO_Q100.dat hoppet-F3WpmZ.dat' u 2:($6/$3) w l t 'hoppet/Apfel',\
     '<paste MSTW_F3e+p_LO_Q100.dat original-F3WpmZ.dat' u 2:($6/$3) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, NLO(only) W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste MSTW_F3e-p_LO_Q100.dat MSTW_F3e-p_NLO_Q100.dat' u 2:(abs(($6-$3))) w l t 'Apfel',\
     'hoppet-F3WpmZ.dat' u 1:(abs(($4))) w l lw 2 t 'hoppet',\
     'original-F3WpmZ.dat' u 1:(abs($4)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.98:1.02]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F3e-p_LO_Q100.dat MSTW_F3e-p_NLO_Q100.dat hoppet-F3WpmZ.dat' u 2:(($10)/($6-$3)) w l t 'hoppet/Apfel',\
     '<paste MSTW_F3e-p_LO_Q100.dat MSTW_F3e-p_NLO_Q100.dat original-F3WpmZ.dat' u 2:(($10)/($6-$3)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, NLO(only) W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste MSTW_F3e+p_LO_Q100.dat MSTW_F3e+p_NLO_Q100.dat' u 2:(abs(($6-$3))) w l t 'Apfel',\
     'hoppet-F3WpmZ.dat' u 1:(abs(($5))) w l lw 2 t 'hoppet',\
     'original-F3WpmZ.dat' u 1:(abs($5)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.98:1.02]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F3e+p_LO_Q100.dat MSTW_F3e+p_NLO_Q100.dat hoppet-F3WpmZ.dat' u 2:(($11)/($6-$3)) w l t 'hoppet/Apfel',\
     '<paste MSTW_F3e+p_LO_Q100.dat MSTW_F3e+p_NLO_Q100.dat original-F3WpmZ.dat' u 2:(($11)/($6-$3)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, NNLO(only) W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste MSTW_F3e-p_NLO_Q100.dat MSTW_F3e-p_NNLO_Q100.dat' u 2:(abs(($6-$3))) w l t 'Apfel',\
     'hoppet-F3WpmZ.dat' u 1:(abs(($6))) w l lw 2 t 'hoppet',\
     'original-F3WpmZ.dat' u 1:(abs($6)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.8:1.2]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F3e-p_NLO_Q100.dat MSTW_F3e-p_NNLO_Q100.dat hoppet-F3WpmZ.dat' u 2:(($12)/($6-$3)) w l t 'hoppet/Apfel',\
     '<paste MSTW_F3e-p_NLO_Q100.dat MSTW_F3e-p_NNLO_Q100.dat original-F3WpmZ.dat' u 2:(($12)/($6-$3)) w l dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, NNLO(only) W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste MSTW_F3e+p_NLO_Q100.dat MSTW_F3e+p_NNLO_Q100.dat' u 2:(abs(($6-$3))) w l t 'Apfel',\
     'hoppet-F3WpmZ.dat' u 1:(abs(($7))) w l lw 2 t 'hoppet',\
     'original-F3WpmZ.dat' u 1:(abs($7)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.8:1.2]
plot 1 w l lw 0.5 notitle,\
     '<paste MSTW_F3e+p_NLO_Q100.dat MSTW_F3e+p_NNLO_Q100.dat hoppet-F3WpmZ.dat' u 2:(($13)/($6-$3)) w l t 'hoppet/Apfel',\
     '<paste MSTW_F3e+p_NLO_Q100.dat MSTW_F3e+p_NNLO_Q100.dat original-F3WpmZ.dat' u 2:(($13)/($6-$3)) w l dt 2 t 'Zaro/Apfel'
unset multiplot


set output
