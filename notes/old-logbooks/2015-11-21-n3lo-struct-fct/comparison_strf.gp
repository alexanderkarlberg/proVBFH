# gnuplot file

set term postscript enhanced color dashed

set macros
resety='set auto y; set log y; set format y "10^{%T}"'
diffy='unset log y; set yrange [0.:1]; set format y "% g"'

# do for [sc in '_scQ _scMH _scMIX'] {
sc='_scQ'
filename='compare_struct_funct'.sc.'.ps'
set output filename

set log x
set xlabel 'x'
set xrange[1e-5:1]

unset multiplot
set lmargin at screen 0.1
set grid

set macros

# scalesch ='_mur1.0_muf1.0'.sc.'.dat _mur1.0_muf4.0'.sc.'.dat _mur4.0_muf1.0'.sc.'.dat _mur4.0_muf3.0'.sc.'.dat'
# do for [scale in scalesch] {
scale='_mur1.0_muf1.0'.sc.'.dat'

#======================================================================
# F1 Z
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F1, Z; xmur='.scale[5:5].', xmuf='.scale[12:12]
@resety
plot 'str_fct/hoppet-F1WpmZ'.scale u 1:(abs($13)) w l lw 2 t 'N3LO',\
     'str_fct/hoppet-F1WpmZ'.scale u 1:(abs($12)) w l lw 2 dt 2 t 'NNLO',\
     'str_fct/hoppet-F1WpmZ'.scale u 1:(abs($11)) w l lw 2 dt 2 t 'NLO'

@diffy
plot 1 w l lw 0.5 notitle,\
     '<paste str_fct/hoppet-F1WpmZ'.scale.' str_fct/hoppet-F1WpmZ'.scale u 1:(abs($26/$12)) w l t 'N3LO/NNLO',\
     '<paste str_fct/hoppet-F1WpmZ'.scale.' str_fct/hoppet-F1WpmZ'.scale u 1:(abs($25/$11)) w l t 'NNLO/NLO'

unset multiplot


#======================================================================
# F2 Z
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, Z; xmur='.scale[5:5].', xmuf='.scale[12:12]
@resety
plot 'str_fct/hoppet-F2WpmZ'.scale u 1:(abs($13)) w l lw 2 t 'N3LO',\
     'str_fct/hoppet-F2WpmZ'.scale u 1:(abs($12)) w l lw 2 dt 2 t 'NNLO',\
     'str_fct/hoppet-F2WpmZ'.scale u 1:(abs($11)) w l lw 2 dt 2 t 'NLO'

@diffy
plot 1 w l lw 0.5 notitle,\
     '<paste str_fct/hoppet-F2WpmZ'.scale.' str_fct/hoppet-F2WpmZ'.scale u 1:(abs($26/$12)) w l t 'N3LO/NNLO',\
     '<paste str_fct/hoppet-F2WpmZ'.scale.' str_fct/hoppet-F2WpmZ'.scale u 1:(abs($25/$11)) w l t 'NNLO/NLO'
unset multiplot

#======================================================================
# F3 Z
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, Z; xmur='.scale[5:5].', xmuf='.scale[12:12]
@resety
plot 'str_fct/hoppet-F3WpmZ'.scale u 1:(abs($13)) w l lw 2 t 'N3LO',\
     'str_fct/hoppet-F3WpmZ'.scale u 1:(abs($12)) w l lw 2 dt 2 t 'NNLO',\
     'str_fct/hoppet-F3WpmZ'.scale u 1:(abs($11)) w l lw 2 dt 2 t 'NLO'

@diffy
plot 1 w l lw 0.5 notitle,\
     '<paste str_fct/hoppet-F3WpmZ'.scale.' str_fct/hoppet-F3WpmZ'.scale u 1:(abs($26/$12)) w l t 'N3LO/NNLO',\
     '<paste str_fct/hoppet-F3WpmZ'.scale.' str_fct/hoppet-F3WpmZ'.scale u 1:(abs($25/$11)) w l t 'NNLO/NLO'
unset multiplot

#}


#======================================================================
#======================================================================
#======================================================================


# W+ and W-


# scalesch ='_mur1.0_muf1.0'.sc.'.dat _mur1.0_muf4.0'.sc.'.dat _mur4.0_muf1.0'.sc.'.dat _mur4.0_muf3.0'.sc.'.dat'
# do for [scale in scalesch] {

# scaleold = scale[1:strlen(scale)-4].'_old.dat'

#======================================================================
# F1 W+
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F1, W+ (e-); xmur='.scale[5:5].', xmuf='.scale[12:12]
@resety
plot 'str_fct/hoppet-F1WpmZ'.scale u 1:(abs($8)) w l lw 2 t 'N3LO',\
     'str_fct/hoppet-F1WpmZ'.scale u 1:(abs($6)) w l lw 2 dt 2 t 'NNLO',\
     'str_fct/hoppet-F1WpmZ'.scale u 1:(abs($4)) w l lw 2 dt 2 t 'NLO'

@diffy
plot 1 w l lw 0.5 notitle,\
     '<paste str_fct/hoppet-F1WpmZ'.scale.' str_fct/hoppet-F1WpmZ'.scale u 1:(abs($21/$6)) w l t 'N3LO/NNLO',\
     '<paste str_fct/hoppet-F1WpmZ'.scale.' str_fct/hoppet-F1WpmZ'.scale u 1:(abs($19/$4)) w l t 'NNLO/NLO'
unset multiplot

#======================================================================
# F2 W+
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, W+ (e-); xmur='.scale[5:5].', xmuf='.scale[12:12]
@resety
plot 'str_fct/hoppet-F2WpmZ'.scale u 1:(abs($8)) w l lw 2 t 'N3LO',\
     'str_fct/hoppet-F2WpmZ'.scale u 1:(abs($6)) w l lw 2 dt 2 t 'NNLO',\
     'str_fct/hoppet-F2WpmZ'.scale u 1:(abs($4)) w l lw 2 dt 2 t 'NLO'

@diffy
plot 1 w l lw 0.5 notitle,\
     '<paste str_fct/hoppet-F2WpmZ'.scale.' str_fct/hoppet-F2WpmZ'.scale u 1:(abs($21/$6)) w l t 'N3LO/NNLO',\
     '<paste str_fct/hoppet-F2WpmZ'.scale.' str_fct/hoppet-F2WpmZ'.scale u 1:(abs($19/$4)) w l t 'NNLO/NLO'
unset multiplot

#======================================================================
# F3 W+
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, W+ (e-); xmur='.scale[5:5].', xmuf='.scale[12:12]
@resety
plot 'str_fct/hoppet-F3WpmZ'.scale u 1:(abs($8)) w l lw 2 t 'N3LO',\
     'str_fct/hoppet-F3WpmZ'.scale u 1:(abs($6)) w l lw 2 dt 2 t 'NNLO',\
     'str_fct/hoppet-F3WpmZ'.scale u 1:(abs($4)) w l lw 2 dt 2 t 'NLO'

@diffy
plot 1 w l lw 0.5 notitle,\
      '<paste str_fct/hoppet-F3WpmZ'.scale.' str_fct/hoppet-F3WpmZ'.scale u 1:(abs($21/$6)) w l t 'N3LO/NNLO',\
      '<paste str_fct/hoppet-F3WpmZ'.scale.' str_fct/hoppet-F3WpmZ'.scale u 1:(abs($19/$4)) w l t 'NNLO/NLO'
unset multiplot

#======================================================================
# F1 W-
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F1, W- (e+); xmur='.scale[5:5].', xmuf='.scale[12:12]
@resety
plot 'str_fct/hoppet-F1WpmZ'.scale u 1:(abs($9)) w l lw 2 t 'N3LO',\
     'str_fct/hoppet-F1WpmZ'.scale u 1:(abs($7)) w l lw 2 dt 2 t 'NNLO',\
     'str_fct/hoppet-F1WpmZ'.scale u 1:(abs($5)) w l lw 2 dt 2 t 'NLO'

@diffy
plot 1 w l lw 0.5 notitle,\
     '<paste str_fct/hoppet-F1WpmZ'.scale.' str_fct/hoppet-F1WpmZ'.scale u 1:(abs($22/$7)) w l t 'N3LO/NNLO',\
     '<paste str_fct/hoppet-F1WpmZ'.scale.' str_fct/hoppet-F1WpmZ'.scale u 1:(abs($20/$5)) w l t 'NNLO/NLO'
unset multiplot

#======================================================================
# F2 W-
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F2, W- (e+); xmur='.scale[5:5].', xmuf='.scale[12:12]
@resety
plot 'str_fct/hoppet-F2WpmZ'.scale u 1:(abs($9)) w l lw 2 t 'N3LO',\
     'str_fct/hoppet-F2WpmZ'.scale u 1:(abs($7)) w l lw 2 dt 2 t 'NNLO',\
     'str_fct/hoppet-F2WpmZ'.scale u 1:(abs($5)) w l lw 2 dt 2 t 'NLO'

@diffy
plot 1 w l lw 0.5 notitle,\
     '<paste str_fct/hoppet-F2WpmZ'.scale.' str_fct/hoppet-F2WpmZ'.scale u 1:(abs($22/$7)) w l t 'N3LO/NNLO',\
     '<paste str_fct/hoppet-F2WpmZ'.scale.' str_fct/hoppet-F2WpmZ'.scale u 1:(abs($20/$5)) w l t 'NNLO/NLO'
unset multiplot

#======================================================================
# F3 W-
#======================================================================

#----------------------------------------------------------------------
set multiplot layout 2,1
set title 'F3, W- (e+); xmur='.scale[5:5].', xmuf='.scale[12:12]
@resety
plot 'str_fct/hoppet-F3WpmZ'.scale u 1:(abs($9)) w l lw 2 t 'N3LO',\
     'str_fct/hoppet-F3WpmZ'.scale u 1:(abs($7)) w l lw 2 dt 2 t 'NNLO',\
     'str_fct/hoppet-F3WpmZ'.scale u 1:(abs($5)) w l lw 2 dt 2 t 'NLO'

@diffy
plot 1 w l lw 0.5 notitle,\
     '<paste str_fct/hoppet-F3WpmZ'.scale.' str_fct/hoppet-F3WpmZ'.scale u 1:(abs($22/$7)) w l t 'N3LO/NNLO',\
     '<paste str_fct/hoppet-F3WpmZ'.scale.' str_fct/hoppet-F3WpmZ'.scale u 1:(abs($20/$5)) w l t 'NNLO/NLO'
unset multiplot

# }

set output
# }
