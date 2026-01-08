# gnuplot file

set term postscript enhanced color
filename='apfel-hoppet-comp-nnpdf.ps'
set output filename

set log x
set xlabel 'x'
set xrange[0.001:1]

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
plot 'NNPDF_F2e-p_LO_Q100.dat' u 2:(abs($3)) w l lw 2 t 'Apfel',\
     'hoppet-F2Wpm_nnpdf.dat' u 1:(abs($2)) w l lw 2 t 'hoppet',\
     'original-F2WpmZ_nnpdf.dat' u 1:(abs($2)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.999:1.001]
plot 1 w l lw 0.5 notitle,\
     '<paste NNPDF_F2e-p_LO_Q100.dat hoppet-F2Wpm_nnpdf.dat' u 2:($5/$3) w l lw 2 t 'hoppet/Apfel',\
     '<paste NNPDF_F2e-p_LO_Q100.dat original-F2WpmZ_nnpdf.dat' u 2:($5/$3) w l lw 2 dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, LO W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot 'NNPDF_F2e+p_LO_Q100.dat' u 2:(abs($3)) w l lw 2 t 'Apfel',\
     'hoppet-F2Wpm_nnpdf.dat' u 1:(abs($3)) w l lw 2 t 'hoppet',\
     'original-F2WpmZ_nnpdf.dat' u 1:(abs($3)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.999:1.001]
plot 1 w l lw 0.5 notitle,\
     '<paste NNPDF_F2e+p_LO_Q100.dat hoppet-F2Wpm_nnpdf.dat' u 2:($6/$3) w l lw 2 t 'hoppet/Apfel',\
     '<paste NNPDF_F2e+p_LO_Q100.dat original-F2WpmZ_nnpdf.dat' u 2:($6/$3) w l lw 2 dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NLO(only) W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste NNPDF_F2e-p_LO_Q100.dat NNPDF_F2e-p_NLO_Q100.dat' u 2:(abs($6-$3)) w l lw 2 t 'Apfel',\
     'hoppet-F2Wpm_nnpdf.dat' u 1:(abs($4)) w l lw 2 t 'hoppet',\
     'original-F2WpmZ_nnpdf.dat' u 1:(abs($4)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.99:1.01]
plot 1 w l lw 0.5 notitle,\
     '<paste NNPDF_F2e-p_LO_Q100.dat NNPDF_F2e-p_NLO_Q100.dat hoppet-F2Wpm_nnpdf.dat' u 2:(($10)/($6-$3)) w l lw 2 t 'hoppet/Apfel',\
     '<paste NNPDF_F2e-p_LO_Q100.dat NNPDF_F2e-p_NLO_Q100.dat original-F2WpmZ_nnpdf.dat' u 2:(($10)/($6-$3)) w l lw 2 dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NLO(only) W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste NNPDF_F2e+p_LO_Q100.dat NNPDF_F2e+p_NLO_Q100.dat' u 2:(abs($6-$3)) w l lw 2 t 'Apfel',\
     'hoppet-F2Wpm_nnpdf.dat' u 1:(abs($5)) w l lw 2 t 'hoppet',\
     'original-F2WpmZ_nnpdf.dat' u 1:(abs($5)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.99:1.01]
plot 1 w l lw 0.5 notitle,\
     '<paste NNPDF_F2e+p_LO_Q100.dat NNPDF_F2e+p_NLO_Q100.dat hoppet-F2Wpm_nnpdf.dat' u 2:(($11)/($6-$3)) w l lw 2 t 'hoppet/Apfel',\
     '<paste NNPDF_F2e+p_LO_Q100.dat NNPDF_F2e+p_NLO_Q100.dat original-F2WpmZ_nnpdf.dat' u 2:(($11)/($6-$3)) w l lw 2 dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NNLO(only) W+ (e-)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste NNPDF_F2e-p_NLO_Q100.dat NNPDF_F2e-p_NNLO_Q100.dat' u 2:(abs($6-$3)) w l lw 2 t 'Apfel',\
     'hoppet-F2Wpm_nnpdf.dat' u 1:(abs($6)) w l lw 2 t 'hoppet',\
     'original-F2WpmZ_nnpdf.dat' u 1:(abs($6)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.8:1.1]
plot 1 w l lw 0.5 notitle,\
     '<paste NNPDF_F2e-p_NLO_Q100.dat NNPDF_F2e-p_NNLO_Q100.dat hoppet-F2Wpm_nnpdf.dat' u 2:(($12)/($6-$3)) w l lw 2 t 'hoppet/Apfel',\
     '<paste NNPDF_F2e-p_NLO_Q100.dat NNPDF_F2e-p_NNLO_Q100.dat original-F2WpmZ_nnpdf.dat' u 2:(($12)/($6-$3)) w l lw 2 dt 2 t 'Zaro/Apfel'
unset multiplot

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, NNLO(only) W- (e+)'
set auto y; set log y; set format y '10^{%T}'
plot '<paste NNPDF_F2e+p_NLO_Q100.dat NNPDF_F2e+p_NNLO_Q100.dat' u 2:(abs($6-$3)) w l lw 2 t 'Apfel',\
     'hoppet-F2Wpm_nnpdf.dat' u 1:(abs($7)) w l lw 2 t 'hoppet',\
     'original-F2WpmZ_nnpdf.dat' u 1:(abs($7)) w l lw 2 dt 2 t 'Zaro'

unset log y; set format y '%g'; set yrange [0.:1.2]
plot 1 w l lw 0.5 notitle,\
     '<paste NNPDF_F2e+p_NLO_Q100.dat NNPDF_F2e+p_NNLO_Q100.dat hoppet-F2Wpm_nnpdf.dat' u 2:(($13)/($6-$3)) w l lw 2 t 'hoppet/Apfel',\
     '<paste NNPDF_F2e+p_NLO_Q100.dat NNPDF_F2e+p_NNLO_Q100.dat original-F2WpmZ_nnpdf.dat' u 2:(($13)/($6-$3)) w l lw 2 dt 2 t 'Zaro/Apfel'
unset multiplot

set output
