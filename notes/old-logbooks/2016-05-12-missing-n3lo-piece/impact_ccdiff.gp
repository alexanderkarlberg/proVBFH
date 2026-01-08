set term pdfcairo enhanced color transparent font "Times, 26" size 22cm,20cm

fnpth0='ptH_n3lo-nnlo_imod0.dat '
fnpth1='ptH_n3lo-nnlo_imod1.dat '
fnpth2='ptH_n3lo-nnlo_imod2.dat '
fnyh0='yH_n3lo-nnlo_imod0.dat '
fnyh1='yH_n3lo-nnlo_imod1.dat '
fnyh2='yH_n3lo-nnlo_imod2.dat '

set output 'impact_ccdiff.pdf'
set title 'Uncertainty due to approx coef. fct.'
set xlabel 'ptH [GeV]'
set yrange [0:0.01]
plot '< paste '.fnpth2.fnpth1.fnpth0 u (($1+$2)/2.0):(abs(($9-$6)/$3)) w l lw 2

set xlabel 'yH'
set yrange [0:0.002]
plot '< paste '.fnyh2.fnyh1.fnyh0 u (($1+$2)/2.0):(abs(($9-$6)/$3)) w l lw 2
