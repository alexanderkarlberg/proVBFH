set term pdfcairo enhanced color transparent font "Times, 26" size 16cm,12cm

fn='n3lo_pdfuncert_13tev.dat'
fnBQ5='n3lo_pdfuncert_13tev_schemeB_Q5.dat'
fnBQ8='n3lo_pdfuncert_13tev_schemeB_Q8.dat'
fnBQ10='n3lo_pdfuncert_13tev_schemeB_Q10.dat'
fnBQ20='n3lo_pdfuncert_13tev_schemeB_Q20.dat'
fnBQ40='n3lo_pdfuncert_13tev_schemeB_Q40.dat'

colLO   = '#FFCC11'
colNLO  = '#0099CC'
colNNLO = '#37BC61'
colN3LO = '#FF4500'

set output 'pdf-uncert-hists.pdf'
set grid front
set label 2000 'PDF4LHC15\_nnlo\_mc' at graph 0.03,0.09
set label 2001 'LHC 13 TeV' at graph 0.03,0.20
set xlabel 'p_{t,H} [GeV]'
set ylabel '{/Symbol d}^{PDF}'

set ytics -0.02,0.005,0.02
set yrange [0.0:0.02]
set xrange [0:300]
plot '<mergeidx.pl -f '.fn    .' NNLO.ptH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_A',\
     '<mergeidx.pl -f '.fnBQ5 .' N3LO.ptH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_B(5 GeV)',\
     '<mergeidx.pl -f '.fnBQ8.' N3LO.ptH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_B(8 GeV )',\
     '<mergeidx.pl -f '.fnBQ10.' N3LO.ptH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_B(10 GeV )'

     # '<mergeidx.pl -f '.fnBQ40.' N3LO.ptH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_B(40 GeV )'

# YH plot
set grid front
set label 2000 'PDF4LHC15\_nnlo\_mc' at graph 0.03,0.09
set label 2001 'LHC 13 TeV' at graph 0.03,0.20
set xlabel 'y_{H} [GeV]'
set ylabel '{/Symbol d}^{PDF}'
set ytics -0.02,0.005,0.02
set yrange [0.0:0.03]
set xrange [0:4.5]
plot '<mergeidx.pl -f '.fn    .' NNLO.yH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_A',\
     '<mergeidx.pl -f '.fnBQ5 .' N3LO.yH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_B(5 GeV)',\
     '<mergeidx.pl -f '.fnBQ8.' N3LO.yH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_B(8 GeV )',\
     '<mergeidx.pl -f '.fnBQ10.' N3LO.yH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_B(10 GeV )'

#     '<mergeidx.pl -f '.fnBQ40.' N3LO.yH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_B(40 GeV )'

set output 'pdf-uncert-energy.pdf'

unset grid
set label 2000 'PDF4LHC15\_nnlo\_mc' at graph 0.96,0.10 right
set label 2001 '' at graph 0.97,0.20 right
fn='pdf-uncert_energy-scan.dat'
fnBQ5='pdf-uncert_energy-scan_schemeB_Q5.dat'
fnBQ8='pdf-uncert_energy-scan_schemeB_Q8.dat'
fnBQ10='pdf-uncert_energy-scan_schemeB_Q10.dat'
fnBQ20='pdf-uncert_energy-scan_schemeB_Q20.dat'
fnBQ40='pdf-uncert_energy-scan_schemeB_Q40.dat'
set logscale x
set xtics add (7,13,20,30,50)
set ytics 0.0,0.01,0.03
set mytics 5
set title '{/Symbol d}^{PDF}' offset 0,-0.5
set ylabel ''
set xlabel '{/Symbol \326} s [TeV]' offset 0,0.3
#set key top left Left reverse invert vertical
set key top left reverse Left
set key samplen 1 maxrow 2 
set yrange [0.0:0.03]
set xrange [7:100]
plot '<mergeidx.pl -f '.fnBQ10 .' N3LO.cross' u 1:2 w lines lw 3 lt 4 t ' {/Symbol d}_B(10 GeV)',\
     '<mergeidx.pl -f '.fnBQ8.' N3LO.cross' u 1:2 w lines lw 3 lt 3 t ' {/Symbol d}_B(8 GeV)',\
     '<mergeidx.pl -f '.fnBQ5.' N3LO.cross' u 1:2 w lines lw 3 lt 2 t ' {/Symbol d}_B(5 GeV)',\
     '<mergeidx.pl -f '.fn    .' NNLO.cross' u 1:2 w lines lw 3 lt 1 t ' {/Symbol d}_A'


#     '<mergeidx.pl -f '.fnBQ40.' N3LO.cross' u 1:2 w lines lw 3 t ' {/Symbol d}_B(40 GeV)'
set key default
set output 'pdf-nnlo-uncert.pdf'
set term pdfcairo size 20cm,15cm
unset grid
set label 2000 'PDF4LHC15\_nlo\_mc' at graph 0.97,0.05 right
set label 2001 '' at graph 0.97,0.20 right
fn='pdf-uncert_energy-scan.dat'
fnBQ5='pdf-uncert_energy-scan_schemeB_Q5.dat'
fnBQ8='pdf-uncert_energy-scan_schemeB_Q8.dat'
fnBQ10='pdf-uncert_energy-scan_schemeB_Q10.dat'
fnBQ20='pdf-uncert_energy-scan_schemeB_Q20.dat'
fnBQ40='pdf-uncert_energy-scan_schemeB_Q40.dat'
set logscale x
set xtics add (7,13,20,30,50)
set ytics 0.0,0.01,0.8
set title '{/Symbol d}^{PDF} at NNLO'
set ylabel ''
set xlabel '{/Symbol \326} s [TeV]'
set key top left Left reverse invert
set yrange [0.0:0.08]
set xrange [7:100]
plot '<mergeidx.pl -f '.fn    .' NNLO.cross' u 1:(2.0*$2) w lines lw 3 t ' 2{/Symbol d}_A',\
     '<mergeidx.pl -f '.fnBQ5 .' NNLO.cross' u 1:2 w lines lw 3 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (5 GeV)',\
     '<mergeidx.pl -f '.fnBQ8.' NNLO.cross' u 1:2 w lines lw 3 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (8 GeV)',\
     '<mergeidx.pl -f '.fnBQ10.' NNLO.cross' u 1:2 w lines lw 3 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (10 GeV)'

#     '<mergeidx.pl -f '.fnBQ40.' NNLO.cross' u 1:2 w lines lw 3 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (40 GeV)'


fn='n3lo_pdfuncert_13tev.dat'
fnBQ5='nnlo_pdfuncert_13tev_schemeB_Q5.dat'
fnBQ8='nnlo_pdfuncert_13tev_schemeB_Q8.dat'
fnBQ10='nnlo_pdfuncert_13tev_schemeB_Q10.dat'
fnBQ20='nnlo_pdfuncert_13tev_schemeB_Q20.dat'
fnBQ40='nnlo_pdfuncert_13tev_schemeB_Q40.dat'

set grid front
unset xtics
unset logscale x
set xtics
set label 2000 'PDF4LHC15\_nlo\_mc' at graph 0.97,0.9 right
set label 2001 'LHC 13 TeV' at graph 0.97,0.8 right
set xlabel 'p_{t,H} [GeV]'
print '<mergeidx.pl -f '.fnBQ10.' NNLO.ptH'
set ytics -0.02,0.01,0.10
set yrange [0.0:0.05]
set xrange [0:300]
plot '<mergeidx.pl -f '.fn    .' NNLO.ptH' u (($1+$2)/2):(2.0*$3) w lines lw 2.0 t ' 2{/Symbol d}_A',\
     '<mergeidx.pl -f '.fnBQ5 .' NNLO.ptH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (5 GeV)',\
     '<mergeidx.pl -f '.fnBQ8.' NNLO.ptH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (8 GeV )',\
     '<mergeidx.pl -f '.fnBQ10.' NNLO.ptH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (10 GeV )'

#     '<mergeidx.pl -f '.fnBQ40.' NNLO.ptH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (40 GeV )'

# YH plot
set grid front
set label 2000 'PDF4LHC15\_nlo\_mc' at graph 0.97,0.9 right
set label 2001 'LHC 13 TeV' at graph 0.97,0.8
set xlabel 'y_{H} [GeV]'
# set ytics -0.02,0.005,0.02
set yrange [0.0:0.05]
set xrange [0:4.5]
plot '<mergeidx.pl -f '.fn    .' NNLO.yH' u (($1+$2)/2):(2.0*$3) w lines lw 2.0 t ' 2{/Symbol d}_A',\
     '<mergeidx.pl -f '.fnBQ5 .' NNLO.yH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (5 GeV)',\
     '<mergeidx.pl -f '.fnBQ8.' NNLO.yH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (8 GeV )',\
     '<mergeidx.pl -f '.fnBQ10.' NNLO.yH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (10 GeV )'

#     '<mergeidx.pl -f '.fnBQ40.' NNLO.yH' u (($1+$2)/2):3 w lines lw 2.0 t ' {/Symbol d}_@{/*0.9B}^{/*0.8(NNLO)} (40 GeV )'

set output
